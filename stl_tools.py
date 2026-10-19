"""Shared helpers for render_views.py and extract_geometry.py: binary STL
loading and .scad -> STL export through the openscad CLI."""
import argparse
import pathlib
import shutil
import subprocess
import sys

import numpy as np

STL_RECORD = np.dtype([
    ("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2"),
])


def load_binary_stl(path):
    data = pathlib.Path(path).read_bytes()
    if len(data) < 84:
        sys.exit(f"{path}: too small to be a binary STL ({len(data)} bytes)")
    count = int(np.frombuffer(data, "<u4", 1, 80)[0])
    expected = 84 + 50 * count
    if len(data) != expected:
        hint = (" (looks like an ASCII STL - export/convert to binary first, "
                 "e.g. `openscad --export-format binstl`)" if data[:5] == b"solid" else "")
        sys.exit(f"{path}: not a binary STL (expected {expected} bytes for {count} "
                  f"triangles, got {len(data)}){hint}")
    records = np.frombuffer(data, STL_RECORD, count, 84)
    return records["v"].astype(np.float64)


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
    return load_binary_stl(path)
