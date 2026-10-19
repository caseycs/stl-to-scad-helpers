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
import itertools
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
    px = (x * scale + imgsize / 2)
    py = (imgsize / 2 - y * scale)

    depth_buf = np.full((imgsize, imgsize), -np.inf)
    color_buf = np.full((imgsize, imgsize, 3), 255, dtype=np.uint8)

    for i in range(len(triangles)):
        tx, ty, tz = px[i], py[i], depth[i]
        xmin = max(int(np.floor(tx.min())), 0)
        xmax = min(int(np.ceil(tx.max())), imgsize - 1)
        ymin = max(int(np.floor(ty.min())), 0)
        ymax = min(int(np.ceil(ty.max())), imgsize - 1)
        if xmin > xmax or ymin > ymax:
            continue

        gx, gy = np.meshgrid(
            np.arange(xmin, xmax + 1) + 0.5, np.arange(ymin, ymax + 1) + 0.5
        )

        denom = (ty[1] - ty[2]) * (tx[0] - tx[2]) + (tx[2] - tx[1]) * (ty[0] - ty[2])
        if denom == 0:
            continue
        w0 = ((ty[1] - ty[2]) * (gx - tx[2]) + (tx[2] - tx[1]) * (gy - ty[2])) / denom
        w1 = ((ty[2] - ty[0]) * (gx - tx[2]) + (tx[0] - tx[2]) * (gy - ty[2])) / denom
        w2 = 1 - w0 - w1

        inside = (w0 >= 0) & (w1 >= 0) & (w2 >= 0)
        if not inside.any():
            continue

        pix_depth = w0 * tz[0] + w1 * tz[1] + w2 * tz[2]
        sub_depth = depth_buf[ymin : ymax + 1, xmin : xmax + 1]
        nearer = inside & (pix_depth > sub_depth)
        if not nearer.any():
            continue

        sub_depth[nearer] = pix_depth[nearer]
        sub_color = color_buf[ymin : ymax + 1, xmin : xmax + 1]
        sub_color[nearer] = colors[i]

    return color_buf


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stl", type=pathlib.Path, help="binary STL file, or an .scad file to export with openscad first")
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

    for name, elev, azim, roll in jobs:
        out_path = args.out / f"{name}.png"
        buf = rasterize(triangles, colors, elev, azim, roll, args.imgsize, frame_points)
        Image.fromarray(buf, "RGB").save(out_path)
        print(out_path)

    print(f"{len(jobs)} images written to {args.out}/")


if __name__ == "__main__":
    main()
