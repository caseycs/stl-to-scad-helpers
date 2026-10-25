#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
#   "pillow",
# ]
# ///
"""Render an STL (or OpenSCAD) model from many angles, one color per surface.

For each rotation axis (elevation, azimuth, roll) sweeps in --step
increments over 0-360 degrees, producing one PNG per angle. Each planar
face is colored by its normal direction, so coplanar triangles (i.e. the
same real-world surface) share a color across all views. Hidden surfaces
are removed with a real per-pixel z-buffer (not matplotlib's painter's
algorithm, which can show occluded faces through nearer ones). Meant to
be fed to a vision model (e.g. Haiku) to inspect the model from many sides.
"""
import argparse
import concurrent.futures
import itertools
import os
import pathlib
import sys

import numpy as np
from PIL import Image

from stl_tools import add_define_option, load_mesh


def parse_box(text):
    values = [float(v) for v in text.split(",")]
    if len(values) != 6:
        raise argparse.ArgumentTypeError("--focus needs 6 values: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX")
    return values


def focus_triangles(triangles, box):
    # keep triangles whose centroid lies inside the box (STL coordinates)
    xmin, xmax, ymin, ymax, zmin, zmax = box
    c = triangles.mean(axis=1)
    inside = (
        (c[:, 0] >= xmin) & (c[:, 0] <= xmax)
        & (c[:, 1] >= ymin) & (c[:, 1] <= ymax)
        & (c[:, 2] >= zmin) & (c[:, 2] <= zmax)
    )
    if not inside.any():
        sys.exit(f"no triangles inside focus box {box}")
    return triangles[inside]


def face_colors(triangles):
    # Derive normals from vertex geometry (not the file's stored normals) and
    # round them before mapping to color, so coplanar triangles always land
    # on the exact same RGB value even if the stored normals have noise.
    v1 = triangles[:, 1] - triangles[:, 0]
    v2 = triangles[:, 2] - triangles[:, 0]
    normals = np.cross(v1, v2)
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    lengths[lengths == 0] = 1
    unit = np.round(normals / lengths, 4)
    rgb = unit * 0.5 + 0.5  # map [-1,1] -> [0,1] per axis, distinct per orientation
    return (rgb * 255).astype(np.uint8)


def camera_basis(elev, azim, roll):
    el, az, ro = np.radians([elev, azim, roll])
    forward = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    world_up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(forward, world_up)) > 0.999:
        world_up = np.array([0.0, 1.0, 0.0])
    # right = cross(world_up, forward), not cross(forward, world_up): forward here points
    # from the scene toward the camera (see depth below), and swapping the cross-product
    # order mirrors every image left-right (verified: a camera on +Y facing the origin
    # must see -X, i.e. world West, on its right; the swapped order gave +X instead).
    right = np.cross(world_up, forward)
    right /= np.linalg.norm(right)
    up = np.cross(forward, right)

    # rotate right/up around the view axis (forward) by roll, via Rodrigues' formula
    axis = forward
    for v_name in ("right", "up"):
        v = right if v_name == "right" else up
        rotated = (
            v * np.cos(ro)
            + np.cross(axis, v) * np.sin(ro)
            + axis * np.dot(axis, v) * (1 - np.cos(ro))
        )
        if v_name == "right":
            right = rotated
        else:
            up = rotated
    return forward, right, up


def rasterize(triangles, colors, elev, azim, roll, imgsize, frame_points=None):
    forward, right, up = camera_basis(elev, azim, roll)
    # frame_points (e.g. a --focus box's corners) sets the camera framing independently of
    # which triangles get drawn, so one oversized triangle pulled in by --focus (its centroid
    # falls in the box, but a vertex may sit far outside it) can't force the camera to zoom out.
    frame_points = triangles.reshape(-1, 3) if frame_points is None else frame_points
    center = frame_points.mean(axis=0)
    rel = triangles - center
    frame_rel = frame_points - center

    x = rel @ right
    y = rel @ up
    depth = rel @ forward  # larger = nearer camera

    margin = 1.1
    half_extent = max(np.abs(frame_rel @ right).max(), np.abs(frame_rel @ up).max()) * margin
    scale = (imgsize / 2) / half_extent
    tx = x * scale + imgsize / 2
    ty = imgsize / 2 - y * scale

    depth_buf = np.full(imgsize * imgsize, -np.inf)
    tri_buf = np.full(imgsize * imgsize, -1, dtype=np.int64)  # winning triangle per pixel

    # per-triangle range of pixels whose centres (k + 0.5) fall inside its bounding box,
    # clipped to the image; this drops sub-pixel triangles that cover no pixel centre at all,
    # which in a dense mesh is most of them. Off-screen and degenerate triangles go too.
    xmin = np.maximum(np.ceil(tx.min(axis=1) - 0.5), 0).astype(np.int32)
    xmax = np.minimum(np.floor(tx.max(axis=1) - 0.5), imgsize - 1).astype(np.int32)
    ymin = np.maximum(np.ceil(ty.min(axis=1) - 0.5), 0).astype(np.int32)
    ymax = np.minimum(np.floor(ty.max(axis=1) - 0.5), imgsize - 1).astype(np.int32)
    denom = (ty[:, 1] - ty[:, 2]) * (tx[:, 0] - tx[:, 2]) + (tx[:, 2] - tx[:, 1]) * (ty[:, 0] - ty[:, 2])
    keep = np.flatnonzero((xmin <= xmax) & (ymin <= ymax) & (denom != 0))
    if not len(keep):
        return np.full((imgsize, imgsize, 3), 255, dtype=np.uint8)

    # barycentric weights as affine functions of the pixel offset from vertex 2:
    # w0 = a0*dx + b0*dy, w1 = a1*dx + b1*dy, depth = z2 + w0*dz0 + w1*dz1
    d = denom[keep]
    a0 = (ty[keep, 1] - ty[keep, 2]) / d
    b0 = (tx[keep, 2] - tx[keep, 1]) / d
    a1 = (ty[keep, 2] - ty[keep, 0]) / d
    b1 = (tx[keep, 0] - tx[keep, 2]) / d
    ox, oy = tx[keep, 2], ty[keep, 2]
    z2 = depth[keep, 2]
    dz0 = depth[keep, 0] - z2
    dz1 = depth[keep, 1] - z2
    bx, by = xmin[keep], ymin[keep]
    bw = xmax[keep] - bx + 1
    area = bw * (ymax[keep] - ymin[keep] + 1)

    # Rasterize triangles in chunks of ~FRAGMENTS candidate pixels at once: every
    # triangle expands to all pixels of its bounding box, the inside test and the
    # z-buffer then run as flat array operations instead of a Python loop per triangle.
    FRAGMENTS = 1 << 21
    ends = np.cumsum(area)
    lo = 0
    while lo < len(keep):
        hi = max(int(np.searchsorted(ends, (ends[lo - 1] if lo else 0) + FRAGMENTS, "right")), lo + 1)
        sl = slice(lo, hi)
        n = area[sl]
        local = np.repeat(np.arange(hi - lo, dtype=np.int32), n)
        off = np.arange(n.sum(), dtype=np.int32) - np.repeat((np.cumsum(n) - n).astype(np.int32), n)
        w = bw[sl][local]
        gx = bx[sl][local] + off % w
        gy = by[sl][local] + off // w
        dx = gx + 0.5 - ox[sl][local]
        dy = gy + 0.5 - oy[sl][local]
        w0 = a0[sl][local] * dx + b0[sl][local] * dy
        w1 = a1[sl][local] * dx + b1[sl][local] * dy
        inside = (w0 >= 0) & (w1 >= 0) & (w0 + w1 <= 1)
        local, w0, w1 = local[inside], w0[inside], w1[inside]
        pix = (gy * imgsize + gx)[inside]
        pix_depth = z2[sl][local] + w0 * dz0[sl][local] + w1 * dz1[sl][local]

        # z-buffer merge: only fragments strictly nearer than what earlier chunks drew
        # compete, the nearest one per pixel wins, and an exact depth tie goes to the
        # earlier triangle (same result as drawing the triangles one by one with ">")
        nearer = pix_depth > depth_buf[pix]
        pix, pix_depth, tri = pix[nearer], pix_depth[nearer], keep[lo + local[nearer]]
        np.maximum.at(depth_buf, pix, pix_depth)
        won = pix_depth == depth_buf[pix]
        pix, tri = pix[won], tri[won]
        tri_buf[pix] = len(triangles)
        np.minimum.at(tri_buf, pix, tri)
        lo = hi

    color_buf = np.full((imgsize * imgsize, 3), 255, dtype=np.uint8)
    hit = tri_buf >= 0
    color_buf[hit] = colors[tri_buf[hit]]
    return color_buf.reshape(imgsize, imgsize, 3)


_scene = None  # per-worker-process copy of the mesh and render settings, set by _init_worker


def _init_worker(*scene):
    global _scene
    _scene = scene


def _render_job(job):
    triangles, colors, imgsize, frame_points, out_dir = _scene
    name, elev, azim, roll = job
    out_path = out_dir / f"{name}.png"
    buf = rasterize(triangles, colors, elev, azim, roll, imgsize, frame_points)
    Image.fromarray(buf, "RGB").save(out_path)
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stl", type=pathlib.Path, help="STL file (binary or ASCII), or an .scad file to export with openscad first")
    parser.add_argument("--step", type=int, default=15, help="degrees between shots per axis")
    parser.add_argument("--imgsize", type=int, default=1024, help="output image size (square, px)")
    parser.add_argument("--out", type=pathlib.Path, default=None, help="output directory")
    parser.add_argument(
        "--focus", type=parse_box, default=None, metavar="XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX",
        help="render only triangles whose centre is inside this box (STL coordinates); "
             "the camera frames the box itself",
    )
    add_define_option(parser)
    args = parser.parse_args()

    if not args.stl.exists():
        sys.exit(f"missing {args.stl}")

    if args.out is None:
        args.out = args.stl.parent / f".stl-renders-{args.stl.stem}"
    args.out.mkdir(parents=True, exist_ok=True)
    for old in args.out.glob("*.png"):  # drop stale images from a previous run (e.g. different --step/--focus)
        old.unlink()

    triangles = load_mesh(args.stl, args.out, args.define)
    frame_points = None
    if args.focus is not None:
        triangles = focus_triangles(triangles, args.focus)
        print(f"focused on {len(triangles)} triangles")
        xmin, xmax, ymin, ymax, zmin, zmax = args.focus
        frame_points = np.array(list(itertools.product((xmin, xmax), (ymin, ymax), (zmin, zmax))))
    colors = face_colors(triangles)

    angles = range(0, 360, args.step)
    base_elev, base_azim, base_roll = 25, -60, 0  # default isometric-ish view

    jobs = []
    for a in angles:
        jobs.append((f"elev{a:03d}", a, base_azim, base_roll))
        jobs.append((f"azim{a:03d}", base_elev, a, base_roll))
        jobs.append((f"roll{a:03d}", base_elev, base_azim, a))

    # one view per task across all cores; the mesh is sent to each worker once, not per view
    scene = (triangles, colors, args.imgsize, frame_points, args.out)
    workers = min(os.cpu_count() or 1, len(jobs))
    with concurrent.futures.ProcessPoolExecutor(workers, initializer=_init_worker, initargs=scene) as pool:
        for out_path in pool.map(_render_job, jobs):
            print(out_path)

    print(f"{len(jobs)} images written to {args.out}/")


if __name__ == "__main__":
    main()
