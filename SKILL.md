---
name: stl-to-scad-helpers
description: Inspect an STL or OpenSCAD model's exact geometry and shape when writing, replicating, or editing OpenSCAD code. Use when the user wants dimensions, faces, or holes read off a 3D model, wants to sanity-check a model from multiple angles, or wants to compare two STL/SCAD files.
---

Reading dimensions off a rendered picture is unreliable — an LLM can't
transcribe exact numbers from a PNG. Use `extract_geometry.py` for the
numbers and `render_views.py` only as a visual sanity check. Both scripts
run standalone via `uv run` (see [`README.md`](README.md) for how) and
accept a `.stl` or `.scad` file.

Run every command below with this skill's own absolute directory in the
script path (e.g. `uv run /abs/path/to/this/skill/extract_geometry.py
<model> ...`), not a bare filename — the model file can be anywhere, but a
relative script path only resolves if your shell's cwd happens to be this
directory.

## Steps

1. **Get exact geometry first.**

   ```
   uv run extract_geometry.py <model.stl|model.scad> --out geometry.json
   ```

   `geometry.json` can be large (every hexagon in a honeycomb grid is its
   own hole loop) — for a big model, don't read it whole: list faces by
   area/normal first (e.g. with `jq` or Python), then pull loop data for
   one face of interest, or re-run with `--min-area` to drop small holes
   from the output entirely.

   Each entry in `faces` is one flat surface with its `normal`,
   `bounding_box`, and `loops` (outer boundary + holes, as both 3D and 2D
   polygon points — `frame` gives the 2D coordinate system). Use these to
   write or check OpenSCAD primitives (`cube`, `polygon()` +
   `linear_extrude`, hole positions, etc). The top-level `volume`,
   `volume_sign`, and `watertight` fields flag a broken mesh before you
   trust anything else in the file.

2. **Render views only to sanity-check a shape or orientation visually.**

   ```
   uv run render_views.py <model.stl|model.scad> --step 30
   ```

   Each flat surface gets a distinct color (grouped by normal direction),
   so you can tell surfaces apart across angles. Look at 2-3 named views
   (e.g. `azim000`, `elev090`, `roll090`) rather than scanning all ~72 —
   see [`docs/usage.md`](docs/usage.md) for what each filename means.

   Use `--focus XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX` to zoom the camera on one
   feature instead. Pad the box a little past the feature's own bounding
   box: `--focus` keeps whole triangles by centroid, not a true geometric
   clip, so a too-tight box can drop a real edge of the feature or, more
   rarely, pull in an oversized unrelated triangle whose centroid just
   barely lands inside.

3. **For localized measurements or diffing two models**, use the
   `section`, `slice` and `compare` subcommands of `extract_geometry.py`
   instead of re-deriving them from the full face list — see
   [`docs/usage.md`](docs/usage.md) for their flags. `slice` is the one
   to reach for when you need a feature's exact outline: it cuts with any
   plane (not just X/Y/Z) and returns lines and arcs you can paste into
   `polygon()`; clip it to the feature with `--window`.

Full flag reference, the `--define`/`-D` passthrough for `.scad` variables,
and the view-filename decoding table are in [`docs/usage.md`](docs/usage.md).

## Installing this skill

Clone this repo into a skills directory Claude Code reads from (e.g.
`~/.claude/skills/stl-to-scad-helpers` or a project's `.claude/skills/`),
or add it as a plugin marketplace/skill source. No build step — the
scripts run as-is via `uv run`.
