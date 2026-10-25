"""Shared helpers for render_views.py and extract_geometry.py: STL loading
(binary or ASCII) and .scad -> STL export through the openscad CLI."""
import argparse
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

STL_RECORD = np.dtype([
    ("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2"),
])


def load_stl(path):
    """Triangles (n, 3, 3) from a binary or ASCII STL."""
    data = pathlib.Path(path).read_bytes()
    if len(data) >= 84:
        count = int(np.frombuffer(data, "<u4", 1, 80)[0])
        if len(data) == 84 + 50 * count:
            return np.frombuffer(data, STL_RECORD, count, 84)["v"].astype(np.float64)
    # a binary header may also start with "solid", so the size check above comes first
    if data.lstrip()[:5] == b"solid":
        coords = re.findall(rb"vertex\s+(\S+)\s+(\S+)\s+(\S+)", data)
        if coords and len(coords) % 3 == 0:
            return np.array(coords).astype(np.float64).reshape(-1, 3, 3)
    sys.exit(f"{path}: neither a binary nor an ASCII STL ({len(data)} bytes)")


def parse_define(text):
    name, sep, value = text.partition("=")
    if not sep or not name.strip():
        raise argparse.ArgumentTypeError(f"--define needs NAME=VALUE, got {text!r}")
    return text


def add_define_option(parser):
    parser.add_argument(
        "--define", "-D", action="append", default=[], type=parse_define, metavar="NAME=VALUE",
        help="OpenSCAD variable override for .scad input (repeatable), e.g. -D 'part=\"left\"'",
    )


def export_scad_to_stl(scad, out, defines=()):
    # arguments go to openscad as a list (no shell), so quotes in values survive
    if shutil.which("openscad") is None:
        sys.exit("openscad not found on PATH")
    cmd = ["openscad", "--export-format", "binstl", "-o", str(out)]
    for d in defines:
        cmd += ["-D", d]
    subprocess.run(cmd + [str(scad)], check=True)
    print(f"exported {scad} -> {out}", file=sys.stderr)
    return pathlib.Path(out)


def load_mesh(path, scad_out_dir, defines=()):
    """Triangles from an STL, or from a .scad exported into scad_out_dir first."""
    path = pathlib.Path(path)
    if not path.exists():
        sys.exit(f"missing {path}")
    if path.suffix == ".scad":
        # distinct suffix so this never silently overwrites a real same-named .stl
        out = pathlib.Path(scad_out_dir) / f"{path.stem}.scad-export.stl"
        path = export_scad_to_stl(path, out, defines)
    return load_stl(path)
