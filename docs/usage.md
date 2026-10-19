# Full usage

`render_views.py` writes the exported STL (from a `.scad` input) into its
`--out` directory; `extract_geometry.py` writes it as `<scad filename>.stl`
next to where you ran the command.

## render_views.py

```
uv run render_views.py <path-to.stl|.scad> [--step DEGREES] [--imgsize PX] [--out DIR] [--focus XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX]
```

- `--step` (default `15`): degrees between shots on each rotation axis.
- `--imgsize` (default `1024`): output image size in pixels (square).
- `--out` (default `.stl-renders-<stl filename>` next to the STL): output directory.
- `--focus` (optional): only render triangles whose centre lies inside this
  box (STL coordinates), and frame the camera on just those triangles —
  useful to zoom in on one feature (e.g. a handle or a stop block).

## extract_geometry.py

```
uv run extract_geometry.py <path-to.stl|.scad> [--decimals N] [--out file.json]
```

- `--decimals` (default `4`): rounding used both to merge shared vertices
  into planes/loops and in the output numbers.
- `--out` (default: prints to stdout): write the JSON to a file instead.

Output is a JSON object with the overall bounding box and a `faces` list
(largest area first). Each face groups every triangle sharing the same
plane (normal + offset) and reports:

- `normal`, `plane_offset`: the plane equation.
- `area`, `triangle_count`, `bounding_box`.
- `loops`: boundary polygons on that plane — one `"outer"` loop (the
  face's outline) plus any `"hole"` loops (e.g. each hexagon in a
  honeycomb grid, or a handle cutout), each as both `points_3d` (model
  coordinates) and `points_2d` (coordinates in the face's own 2D plane,
  ready to drop into an OpenSCAD `polygon()`).

Planes are grouped with a tolerance (normal 2e-3, offset 0.02 mm), so one
flat face is one entry; a plane with several separate regions has several
`"outer"` loops (holes are decided by loop winding, not size). The top
level also reports `volume`, `volume_sign`, `watertight`, `bad_edge_count`.

### `--define` for `.scad` input

Both scripts take `--define NAME=VALUE` / `-D` (repeatable) for `.scad`
input, passed to `openscad -D` without a shell, e.g. `-D 'part="left"'`.

### `section` subcommand

```
extract_geometry.py section FILE --cut Y=40,160 [--scan Z=-13] [--profile X] [--prominence 0.2] [--shift X,Y,Z]
```

Material intervals along each scan line, and peaks/troughs of the low and
high boundary (exact, at section vertices). Values that hit a vertex are
nudged by `1e-6` with a warning.

### `compare` subcommand

```
extract_geometry.py compare A B [--mirror-x] [--shift X,Y,Z] [--cut ... --scan ...]
```

Volume, bounding box, area per normal, and section differences; mirror and
shift apply to `B`, there is no automatic alignment.

## How surfaces are colored

Each triangle is colored by its normal direction (mapped from `[-1,1]` to
`[0,1]` per axis), so every triangle belonging to the same flat surface gets
the exact same color, and differently-oriented surfaces get different
colors. Hidden surfaces are removed with a real per-pixel z-buffer, so
occluded faces never show through nearer ones.

## Decoding view settings from filenames

Each image sweeps one of three camera rotations while holding the other two
at their base value (`elev=25`, `azim=-60`, `roll=0`):

| Filename pattern | Meaning |
| --- | --- |
| `elevNNN.png` | elevation = `NNN` degrees, azim = -60, roll = 0 |
| `azimNNN.png` | azimuth = `NNN` degrees, elev = 25, roll = 0 |
| `rollNNN.png` | roll = `NNN` degrees, elev = 25, azim = -60 |

`NNN` is the swept angle in degrees (0, step, 2*step, ... up to but
excluding 360). For example `azim090.png` is the view with elevation 25°,
azimuth 90°, roll 0°.

- **elevation**: angle of the camera above the horizontal plane.
- **azimuth**: angle of the camera around the vertical axis.
- **roll**: rotation of the camera around its own viewing direction (tilts
  the horizon).
