---
name: stl-to-scad-helpers
description: Inspect an STL or OpenSCAD model's exact geometry and shape when writing, replicating, or editing OpenSCAD code. Use when the user wants dimensions, faces, or holes read off a 3D model, wants to sanity-check a model from multiple angles, or wants to compare two STL/SCAD files.
---

Reading dimensions off a rendered picture is unreliable — an LLM can't
transcribe exact numbers from a PNG. Use `extract_geometry.py` for the
numbers and `render_views.py` only as a visual sanity check.

Both scripts are `uv run`-able standalone (dependencies are declared
inline in each file via PEP 723) — no install step. They accept a `.stl`
or a `.scad` file (the latter is exported through the `openscad` CLI,
which must be on `PATH`).

## Steps

1. **Get exact geometry first.**

   ```
   uv run extract_geometry.py <model.stl|model.scad> --out geometry.json
   ```

   Read `geometry.json`: each entry in `faces` is one flat surface with its
   `normal`, `bounding_box`, and `loops` (outer boundary + holes, as both
   3D and 2D polygon points). Use these to write or check OpenSCAD
   primitives (`cube`, `polygon()` + `linear_extrude`, hole positions,
   etc). The top-level `volume`, `volume_sign`, and `watertight` fields
   flag a broken mesh before you trust anything else in the file.

2. **Render views only to sanity-check a shape or orientation visually.**

   ```
   uv run render_views.py <model.stl|model.scad> --step 30
   ```

   Each flat surface gets a distinct color (grouped by normal direction),
   so you can tell surfaces apart across angles. Use `--focus
   XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX` to zoom the camera on one feature
   instead of scanning every rendered PNG. Decoding which angle each PNG
   filename shows is in [`docs/usage.md`](docs/usage.md).

3. **For localized measurements or diffing two models**, use the
   `section` and `compare` subcommands of `extract_geometry.py` instead of
   re-deriving them from the full face list — see
   [`docs/usage.md`](docs/usage.md) for their flags.

Full flag reference, the `--define`/`-D` passthrough for `.scad` variables,
and the view-filename decoding table are in [`docs/usage.md`](docs/usage.md).

## Installing this skill

Clone this repo into a skills directory Claude Code reads from (e.g.
`~/.claude/skills/stl-to-scad-helpers` or a project's `.claude/skills/`),
or add it as a plugin marketplace/skill source. No build step — the
scripts run as-is via `uv run`.
