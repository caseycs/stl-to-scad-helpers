# stl-to-scad-heplers

Tools for an LLM to inspect an STL (or OpenSCAD) model when replicating or
editing it in OpenSCAD:

- **`render_views.py`** renders the model from many camera angles into PNGs,
  one color per surface, for a vision model to eyeball from all sides.
- **`extract_geometry.py`** extracts exact planar-face geometry (normals,
  dimensions, boundary polygons) as JSON. A vision model can't read exact
  dimensions off a picture, so use this for replicating precise measurements
  in OpenSCAD, and keep the renders as a visual sanity check only.

Both accept a `.scad` file too: it's exported to binary STL first via the
`openscad` CLI (must be on `PATH`).

Each script is a [PEP 723](https://peps.python.org/pep-0723/) single-file
script with its dependencies declared inline — run with `uv run`, no
virtualenv setup needed.

```
uv run render_views.py model.stl --step 30
uv run extract_geometry.py model.scad --out geometry.json
```

See [`docs/usage.md`](docs/usage.md) for full options (focus boxes,
`section`/`compare` subcommands, `--define` for `.scad` variables) and how
view-angle filenames decode. See [`SKILL.md`](SKILL.md) to use this as a
Claude Code skill.

## Requirements

- Python ≥ 3.10, [`uv`](https://docs.astral.sh/uv/) on `PATH`.
- [`openscad`](https://openscad.org/) on `PATH`, only if you pass a `.scad`
  file instead of an `.stl`.
