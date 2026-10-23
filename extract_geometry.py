#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
# ]
# ///
"""Extract structured planar-face geometry from a binary STL file as JSON.

Pixel renders are good for a sanity check but a vision model can't read
exact dimensions off them. This instead groups triangles into flat
surfaces (same plane) and, for each one, reports its normal, plane
equation, bounding box, and boundary loops (outer outlines plus any holes)
as exact 2D/3D coordinates - numbers an LLM can transcribe directly into
OpenSCAD primitives (cube, polygon+linear_extrude, etc). The top level also
reports volume and a watertight check.

Subcommands (see `extract_geometry.py section --help`):
  section  cut the mesh with axis-aligned planes: material intervals along
           scan lines and boundary peaks/troughs (profile)
  compare  volume, bounding box, per-normal area and section differences
           between two meshes
"""
import argparse
import json
import pathlib
import sys

import numpy as np

from stl_tools import add_define_option, load_mesh, parse_define

AXES = "xyz"
NUDGE = 1e-6


def rnd(values, decimals):
    # + 0.0 turns -0.0 into 0.0
    if np.ndim(values) == 0:
        return round(float(values), decimals) + 0.0
    return (np.round(np.asarray(values, dtype=float), decimals) + 0.0).tolist()


def bbox(pts, decimals):
    return {a: [rnd(pts[:, i].min(), decimals), rnd(pts[:, i].max(), decimals)] for i, a in enumerate(AXES)}


def dedupe_vertices(triangles, decimals=4):
    pts = triangles.reshape(-1, 3)
    rounded = np.round(pts, decimals)
    uniq, inverse = np.unique(rounded, axis=0, return_inverse=True)
    return uniq, inverse.reshape(-1, 3)


def face_normals_and_areas(triangles):
    v1 = triangles[:, 1] - triangles[:, 0]
    v2 = triangles[:, 2] - triangles[:, 0]
    cross = np.cross(v1, v2)
    lengths = np.linalg.norm(cross, axis=1)
    areas = lengths / 2
    safe = np.where(lengths == 0, 1, lengths)
    unit = cross / safe[:, None]
    return unit, areas


def group_coplanar(triangles, normals, areas, normal_tol=2e-3, offset_tol=0.02):
    # Greedy clustering, largest triangle first: a triangle joins a plane when
    # its normal is within normal_tol (per component) and all its vertices lie
    # within offset_tol of that plane (the earliest such plane wins).
    # offset_tol=inf groups by normal only.
    # Planes are bucketed by normal into cells of size normal_tol, so a
    # candidate within normal_tol is always in one of the 27 neighbouring cells.
    check_offset = np.isfinite(offset_tol)
    cells = np.floor(normals / normal_tol).astype(int).tolist()
    normals_l, tris_l = normals.tolist(), triangles.tolist()
    rep_n, rep_d, buckets = [], [], {}
    neighbours = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    group = np.empty(len(triangles), dtype=int)
    for i in np.argsort(-areas, kind="stable").tolist():
        (nx, ny, nz), (cx, cy, cz) = normals_l[i], cells[i]
        best = None
        for dx, dy, dz in neighbours:
            for r in buckets.get((cx + dx, cy + dy, cz + dz), ()):
                if best is not None and r > best:
                    continue
                rx, ry, rz = rep_n[r]
                if abs(rx - nx) >= normal_tol or abs(ry - ny) >= normal_tol or abs(rz - nz) >= normal_tol:
                    continue
                if check_offset and any(abs(rx * x + ry * y + rz * z - rep_d[r]) >= offset_tol
                                        for x, y, z in tris_l[i]):
                    continue
                best = r
        if best is None:
            best = len(rep_n)
            x, y, z = tris_l[i][0]
            rep_n.append((nx, ny, nz))
            rep_d.append(nx * x + ny * y + nz * z)
            buckets.setdefault((cx, cy, cz), []).append(best)
        group[i] = best
    return group


def split_groups(group_ids):
    """(group id, triangle indices) per group, in increasing group id order."""
    order = np.argsort(group_ids, kind="stable")
    ids, starts = np.unique(group_ids[order], return_index=True)
    return zip(ids.tolist(), np.split(order, starts[1:]))


# (u, v) for the six axis-aligned normals, world-axis-aligned so e.g. a +Z face's
# points_2d reads directly as (x, y). v always equals cross(normal, u) (see in_plane_basis).
_AXIS_FRAMES = {
    (0, 0, 1): (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
    (0, 0, -1): (np.array([1.0, 0.0, 0.0]), np.array([0.0, -1.0, 0.0])),
    (1, 0, 0): (np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    (-1, 0, 0): (np.array([0.0, -1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    (0, 1, 0): (np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])),
    (0, -1, 0): (np.array([0.0, 0.0, 1.0]), np.array([-1.0, 0.0, 0.0])),
}


def in_plane_basis(normal):
    # right-handed (u, v, normal): loops wound CCW around the normal have positive 2D area.
    key = (int(round(normal[0])), int(round(normal[1])), int(round(normal[2])))
    if key in _AXIS_FRAMES:
        return _AXIS_FRAMES[key]
    helper = np.array([0.0, 0.0, 1.0]) if abs(normal[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(normal, helper)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    return u, v


def boundary_loops(tri_vert_ids):
    # Directed edges keep the triangles' winding: an edge whose reverse is not
    # in the group is boundary, and following them traces outlines
    # counter-clockwise and holes clockwise (seen from the normal side).
    edges = {(p, q) for a, b, c in tri_vert_ids for p, q in ((a, b), (b, c), (c, a))}
    outgoing = {}
    for p, q in edges:
        if (q, p) not in edges:
            outgoing.setdefault(p, []).append(q)

    loops = []
    for start in list(outgoing):
        while outgoing.get(start):
            loop = [start]
            current = outgoing[start].pop()
            while current != start and outgoing.get(current):
                loop.append(current)
                current = outgoing[current].pop()
            if current != start:
                print(f"warning: open boundary loop ({len(loop)} points) - mesh is not "
                      f"watertight here", file=sys.stderr)
            loops.append(loop)
    return loops


def polygon_area_2d(points):
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def extract_faces(triangles, decimals=4):
    normals, areas = face_normals_and_areas(triangles)
    degenerate = areas <= 1e-9
    if degenerate.any():
        # a zero-area sliver would otherwise form its own group with a zero/NaN normal
        # and crash the weighted average below (weights summing to zero)
        print(f"warning: dropping {int(degenerate.sum())} degenerate (zero-area) "
              f"triangle(s)", file=sys.stderr)
        triangles, normals, areas = triangles[~degenerate], normals[~degenerate], areas[~degenerate]

    vertices, tri_vert_ids = dedupe_vertices(triangles, decimals)
    group_ids = group_coplanar(triangles, normals, areas)

    faces = []
    for g, mask in split_groups(group_ids):
        group_tris = tri_vert_ids[mask]
        normal = (normals[mask] * areas[mask, None]).sum(axis=0)
        normal /= np.linalg.norm(normal)
        offset = np.average(triangles[mask].mean(axis=1) @ normal, weights=areas[mask])
        pts = vertices[np.unique(group_tris)]

        loops = []
        u, v = in_plane_basis(normal)
        origin = vertices[group_tris[0, 0]]
        for loop_ids in boundary_loops(group_tris):
            loop_pts = vertices[loop_ids]
            rel = loop_pts - origin
            uv = np.column_stack([rel @ u, rel @ v])
            loops.append({
                "kind": "outer" if polygon_area_2d(uv) > 0 else "hole",
                "points_3d": rnd(loop_pts, decimals),
                "points_2d": rnd(uv, decimals),
            })
        loops.sort(key=lambda loop: loop["kind"] != "outer")

        faces.append({
            "id": int(g),
            "normal": rnd(normal, decimals),
            "plane_offset": rnd(offset, decimals),
            "area": rnd(areas[mask].sum(), decimals),
            "triangle_count": len(mask),
            "bounding_box": bbox(pts, decimals),
            # points_2d = origin + u*p[0] + v*p[1] in model coordinates
            "frame": {"origin": rnd(origin, decimals), "u": rnd(u, decimals), "v": rnd(v, decimals)},
            "loops": loops,
        })

    faces.sort(key=lambda f: -f["area"])
    return faces


def mesh_stats(triangles, decimals=4):
    # signed volume: sum of tetrahedra to the origin, positive for outward normals
    signed = np.einsum("ij,ij->i", triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum() / 6
    _, ids = dedupe_vertices(triangles, decimals)
    edges = np.sort(np.concatenate([ids[:, [0, 1]], ids[:, [1, 2]], ids[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    bad = int((counts != 2).sum())
    return {
        "volume": rnd(abs(signed), decimals),
        "volume_sign": 1 if signed >= 0 else -1,
        "watertight": bad == 0,
        "bad_edge_count": bad,
        "triangle_count": len(triangles),
    }


# ---- section ---------------------------------------------------------------

def parse_axis_values(text):
    axis, sep, values = text.partition("=")
    axis = axis.strip().lower()
    try:
        nums = [float(v) for v in values.split(",")]
    except ValueError:
        nums = []
    if not sep or axis not in AXES or not nums:
        raise argparse.ArgumentTypeError(f"expected AXIS=V1[,V2...] with AXIS x, y or z, got {text!r}")
    return AXES.index(axis), nums


def parse_xyz(text):
    values = [float(v) for v in text.split(",")]
    if len(values) != 3:
        raise argparse.ArgumentTypeError("--shift needs 3 values: X,Y,Z")
    return np.array(values)


def parse_profile_axis(text):
    text = text.strip().lower()
    if text not in AXES:
        raise argparse.ArgumentTypeError(f"--profile expects x, y or z, got {text!r}")
    return AXES.index(text)


def nudge(value, coords, label):
    # move a cut/scan value off any vertex lying on it, so crossings are never degenerate
    original = value
    while len(coords) and np.abs(coords - value).min() < NUDGE / 2:
        value += NUDGE
    if value != original:
        print(f"warning: {label}={original:g} hits a vertex, using {value:.7f}", file=sys.stderr)
    return value


def cut_segments(triangles, normals, axis, value):
    """Segments where the plane coord[axis] == value crosses the mesh, with each triangle's normal."""
    d = triangles[:, :, axis] - value
    side = d > 0
    crossing = side.any(axis=1) & ~side.all(axis=1)
    segs = []
    for tri, dist, s in zip(triangles[crossing], d[crossing], side[crossing]):
        pts = [tri[i] + dist[i] / (dist[i] - dist[j]) * (tri[j] - tri[i])
               for i, j in ((0, 1), (1, 2), (2, 0)) if s[i] != s[j]]
        segs.append(pts)
    return np.array(segs).reshape(-1, 2, 3), normals[crossing]


def scan_intervals(segs, seg_normals, scan_axis, value, along):
    """Material intervals along `along` on the line coord[scan_axis] == value.

    A crossing enters material when the surface normal points against the
    scan direction (normal[along] < 0) and leaves when it points with it.
    """
    a, b = segs[:, 0], segs[:, 1]
    hit = (a[:, scan_axis] - value) * (b[:, scan_axis] - value) < 0
    t = (value - a[hit, scan_axis]) / (b[hit, scan_axis] - a[hit, scan_axis])
    pos = a[hit, along] + t * (b[hit, along] - a[hit, along])
    sign = np.sign(seg_normals[hit, along])
    intervals, depth, start = [], 0, None
    for p, s in sorted(zip(pos, sign)):
        if s < 0:
            depth += 1
            if depth == 1:
                start = p
        elif s > 0:
            depth -= 1
            if depth == 0:
                intervals.append([start, p])
            elif depth < 0:
                print(f"warning: unpaired exit at {p:.4f}", file=sys.stderr)
                depth = 0
    if depth:
        print(f"warning: {depth} unpaired entry crossing(s) on scan line", file=sys.stderr)
    return intervals


def turning_points(pos, val, prominence):
    """Local maxima/minima of val(pos) whose rise/fall to the next turn exceeds prominence (zigzag filter)."""
    out, trend, lo, hi, cand = [], 0, 0, 0, 0

    def emit(kind, i):
        j, k = i, i
        while j > 0 and abs(val[j - 1] - val[i]) < 1e-6:
            j -= 1
        while k < len(val) - 1 and abs(val[k + 1] - val[i]) < 1e-6:
            k += 1
        if j > 0 and k < len(val) - 1:  # skip extremes that touch the ends
            out.append({"kind": kind, "value": val[i], "at": [pos[j], pos[k]] if j != k else pos[i]})

    for i in range(1, len(val)):
        if trend == 0:
            lo = i if val[i] < val[lo] else lo
            hi = i if val[i] > val[hi] else hi
            if val[i] - val[lo] >= prominence:
                emit("min", lo)
                trend, cand = 1, i
            elif val[hi] - val[i] >= prominence:
                emit("max", hi)
                trend, cand = -1, i
        elif trend == 1:
            if val[i] > val[cand]:
                cand = i
            elif val[cand] - val[i] >= prominence:
                emit("max", cand)
                trend, cand = -1, i
        else:
            if val[i] < val[cand]:
                cand = i
            elif val[i] - val[cand] >= prominence:
                emit("min", cand)
                trend, cand = 1, i
    return out


def section_profile(segs, along, by, prominence, decimals):
    """Low and high boundary (min/max of `along`) as a function of the other in-plane axis `by`.

    The boundary is piecewise linear, so sampling it at every section vertex is exact.
    """
    a, c = segs[:, 0], segs[:, 1]
    breaks = np.unique(segs[:, :, by])
    lo_hi = []
    for b in breaks:
        on = (np.minimum(a[:, by], c[:, by]) <= b) & (np.maximum(a[:, by], c[:, by]) >= b)
        p, q = a[on], c[on]
        span = q[:, by] - p[:, by]
        t = (b - p[:, by]) / np.where(span == 0, 1, span)
        vals = np.concatenate([p[:, along] + t * (q[:, along] - p[:, along]), q[span == 0, along]])
        lo_hi.append((vals.min(), vals.max()))
    lo_hi = np.array(lo_hi)
    result = {"along": AXES[along], "by": AXES[by], "prominence": prominence}
    for name, col in (("low", 0), ("high", 1)):
        result[name] = [{k: v if k == "kind" else rnd(v, decimals) for k, v in tp.items()}
                        for tp in turning_points(breaks, lo_hi[:, col], prominence)]
    return result


def section_cut(triangles, normals, axis, value, scan, profile, prominence, decimals):
    value = nudge(value, np.unique(triangles[:, :, axis]), AXES[axis])
    segs, seg_normals = cut_segments(triangles, normals, axis, value)
    result = {"axis": AXES[axis], "value": rnd(value, 7), "segment_count": len(segs)}
    if not len(segs):
        return result
    result["bounding_box"] = bbox(segs.reshape(-1, 3), decimals)

    if scan is not None:
        scan_axis, scan_values = scan
        along = 3 - axis - scan_axis
        result["scans"] = []
        for sv in scan_values:
            sv = nudge(sv, np.unique(segs[:, :, scan_axis]), AXES[scan_axis])
            intervals = scan_intervals(segs, seg_normals, scan_axis, sv, along)
            result["scans"].append({
                "line": {AXES[scan_axis]: rnd(sv, 7)},
                "along": AXES[along],
                "intervals": rnd(intervals, decimals),
            })

    if profile is not None:
        result["profile"] = section_profile(segs, profile, 3 - axis - profile, prominence, decimals)
    return result


def transformed(triangles, shift=None, mirror_x=False):
    if mirror_x:
        triangles = triangles[:, [0, 2, 1]] * np.array([-1.0, 1.0, 1.0])  # mirror and keep normals outward
    if shift is not None:
        triangles = triangles + shift
    return triangles


def run_sections(triangles, args):
    normals, _ = face_normals_and_areas(triangles)
    axis, values = args.cut
    if args.scan is not None and args.scan[0] == axis:
        sys.exit("--scan axis must differ from --cut axis")
    if args.profile is not None and args.profile == axis:
        sys.exit("--profile axis must differ from --cut axis")
    return [section_cut(triangles, normals, axis, v, args.scan, args.profile, args.prominence, args.decimals)
            for v in values]


def add_section_options(parser):
    parser.add_argument("--cut", type=parse_axis_values, required=True, metavar="AXIS=V1[,V2...]",
                        help="cutting plane(s), e.g. Y=40,100,160")
    parser.add_argument("--scan", type=parse_axis_values, default=None, metavar="AXIS=V1[,V2...]",
                        help="scan lines inside each cut, e.g. Z=-13,-11; reports material intervals "
                             "along the remaining axis")
    parser.add_argument("--profile", type=parse_profile_axis, default=None, metavar="AXIS",
                        help="report peaks/troughs of the low and high boundary along AXIS "
                             "(e.g. X) as a function of the remaining in-plane axis")
    parser.add_argument("--prominence", type=float, default=0.2,
                        help="minimum rise/fall (mm) for a profile peak or trough (default 0.2)")
    parser.add_argument("--shift", type=parse_xyz, default=None, metavar="X,Y,Z",
                        help="translate the mesh by this before cutting")
    parser.add_argument("--decimals", type=int, default=4, help="output rounding")
    add_define_option(parser)


def write_json(result, out):
    text = json.dumps(result, indent=2)
    if out is None:
        print(text)
    else:
        out.write_text(text)
        print(f"written to {out}", file=sys.stderr)


def section_main(argv):
    parser = argparse.ArgumentParser(prog="extract_geometry.py section", description=(
        "Cut a mesh with axis-aligned planes. Values that hit a vertex are nudged by 1e-6 (warning on stderr)."))
    parser.add_argument("stl", type=pathlib.Path, help="binary STL file, or an .scad file to export first")
    add_section_options(parser)
    parser.add_argument("--out", type=pathlib.Path, default=None, help="output JSON path (default: stdout)")
    args = parser.parse_args(argv)
    triangles = transformed(load_mesh(args.stl, pathlib.Path.cwd(), args.define), args.shift)
    write_json({"source": str(args.stl), "sections": run_sections(triangles, args)}, args.out)


# ---- compare ---------------------------------------------------------------

def compare_main(argv):
    parser = argparse.ArgumentParser(prog="extract_geometry.py compare", description=(
        "Compare two meshes: volume, bounding box, area per normal direction, and optionally "
        "sections. --mirror-x and --shift apply to B (mirror first); there is no automatic alignment."))
    parser.add_argument("a", type=pathlib.Path)
    parser.add_argument("b", type=pathlib.Path)
    parser.add_argument("--mirror-x", action="store_true", help="mirror B in X (x -> -x) before shifting")
    parser.add_argument("--normal-tol", type=float, default=2e-3,
                        help="per-component normal grouping tolerance for area_by_normal (default 0.002)")
    parser.add_argument("--cut", type=parse_axis_values, default=None, metavar="AXIS=V1[,V2...]",
                        help="also compare sections at these planes (see section --help)")
    parser.add_argument("--scan", type=parse_axis_values, default=None, metavar="AXIS=V1[,V2...]")
    parser.add_argument("--profile", type=parse_profile_axis, default=None, metavar="AXIS")
    parser.add_argument("--prominence", type=float, default=0.2)
    parser.add_argument("--shift", type=parse_xyz, default=None, metavar="X,Y,Z", help="translate B by this")
    parser.add_argument("--decimals", type=int, default=4, help="output rounding")
    parser.add_argument("--out", type=pathlib.Path, default=None, help="output JSON path (default: stdout)")
    add_define_option(parser)
    parser.add_argument(
        "--define-b", "-Db", action="append", default=[], type=parse_define, metavar="NAME=VALUE",
        help="extra OpenSCAD variable override for B only, applied after --define "
             "(e.g. compare one .scad against itself under a different -D)",
    )
    args = parser.parse_args(argv)
    dec = args.decimals

    meshes = {
        "a": load_mesh(args.a, pathlib.Path.cwd(), args.define),
        "b": transformed(load_mesh(args.b, pathlib.Path.cwd(), args.define + args.define_b),
                          args.shift, args.mirror_x),
    }
    result = {"a": str(args.a), "b": str(args.b)}
    for key, tris in meshes.items():
        result[f"{key}_stats"] = mesh_stats(tris, dec) | {"bounding_box": bbox(tris.reshape(-1, 3), dec)}
    result["volume_diff"] = rnd(result["b_stats"]["volume"] - result["a_stats"]["volume"], dec)

    both = np.concatenate([meshes["a"], meshes["b"]])
    normals, areas = face_normals_and_areas(both)
    group = group_coplanar(both, normals, areas, args.normal_tol, np.inf)
    from_b = np.arange(len(both)) >= len(meshes["a"])
    rows = []
    for _, m in split_groups(group):
        n = (normals[m] * areas[m, None]).sum(axis=0)
        area_a, area_b = areas[m[~from_b[m]]].sum(), areas[m[from_b[m]]].sum()
        rows.append({"normal": rnd(n / np.linalg.norm(n), dec), "area_a": rnd(area_a, dec),
                     "area_b": rnd(area_b, dec), "diff": rnd(area_b - area_a, dec)})
    result["area_by_normal"] = sorted(rows, key=lambda r: -max(r["area_a"], r["area_b"]))

    if args.cut is not None:
        sections = {key: run_sections(tris, args) for key, tris in meshes.items()}
        result["sections"] = []
        for sa, sb in zip(sections["a"], sections["b"]):
            entry = {"axis": sa["axis"], "value": sa["value"], "a": sa, "b": sb}
            diffs = []
            for la, lb in zip(sa.get("scans", []), sb.get("scans", [])):
                ia, ib = np.array(la["intervals"]).reshape(-1, 2), np.array(lb["intervals"]).reshape(-1, 2)
                diffs.append({"line": la["line"], "same_count": len(ia) == len(ib),
                              "b_minus_a": rnd(ib - ia, dec) if len(ia) == len(ib) else None})
            entry["scan_diffs"] = diffs
            result["sections"].append(entry)
    write_json(result, args.out)


# ---- faces (default) -------------------------------------------------------

def faces_main(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stl", type=pathlib.Path, help="binary STL file, or an .scad file to export with openscad first")
    parser.add_argument("--decimals", type=int, default=4, help="rounding used to merge vertices and in output")
    parser.add_argument("--out", type=pathlib.Path, default=None, help="output JSON path (default: stdout)")
    parser.add_argument("--min-area", type=float, default=None,
                        help="drop faces smaller than this area (mm^2), e.g. to skip small holes "
                             "when the full output is too large to read in one go")
    add_define_option(parser)
    args = parser.parse_args(argv)

    triangles = load_mesh(args.stl, pathlib.Path.cwd(), args.define)
    faces = extract_faces(triangles, args.decimals)
    if args.min_area is not None:
        faces = [f for f in faces if f["area"] >= args.min_area]

    result = {
        "source": str(args.stl),
        "bounding_box": bbox(triangles.reshape(-1, 3), args.decimals),
        **mesh_stats(triangles, args.decimals),
        "face_count": len(faces),
        "faces": faces,
    }
    write_json(result, args.out)
    if args.out is not None:
        print(f"{len(faces)} planar faces", file=sys.stderr)


def main():
    commands = {"section": section_main, "compare": compare_main}
    if len(sys.argv) > 1 and sys.argv[1] in commands:
        commands[sys.argv[1]](sys.argv[2:])
    else:
        faces_main(sys.argv[1:])


if __name__ == "__main__":
    main()
