# Full usage

## render_views.py

```
uv run render_views.py <path-to.stl|.scad> [--step DEGREES] [--imgsize PX] [--out DIR] [--focus XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX]
```

- `--step` (default `15`): degrees between shots on each rotation axis. The
  default sweeps 72 full-resolution images; for a quick pass on a large or
  fine mesh, start with `--step 30 --imgsize 512` and zoom in with `--focus`
  instead of cranking resolution.
- `--imgsize` (default `1024`): output image size in pixels (square).
- `--out` (default `.stl-renders-<stl filename>` next to the STL): output
  directory. Cleared of old `*.png` files at the start of each run, so a
  later run with different `--step`/`--focus` doesn't leave stale images
  mixed in with the new set.
- `--focus` (optional): only render triangles whose centre lies inside this
  box (STL coordinates), and frame the camera on exactly that box — useful
  to zoom in on one feature (e.g. a handle or a stop block). A negative
  value needs the `=` form, e.g. `--focus=-40,40,-10,10,0,20` (otherwise
  argparse reads the leading `-` as a new flag). Note this selects *whole
  triangles* by centroid, not a true geometric clip: a triangle whose
  centroid falls just inside the box but has a vertex far outside it is
  still drawn in full (and a triangle whose centroid falls just outside is
  dropped even if most of it is inside) — pad the box a bit and expect the
  occasional oversized or missing neighbour at the edges.

If `render_views.py` is run on a `.scad` file, the exported STL is written
into `--out` as `<scad filename>.scad-export.stl`.

## Camera convention

Azimuth 0°/elevation 0° is a camera sitting on the model's +X axis, looking
back toward the origin — so you're looking at the model's -X side. Azimuth
rotates the camera around the vertical (Z) axis from there (90° moves it to
+Y, 180° to -X, and so on). Elevation 0° is level with the XY plane; 90° is
straight overhead (top-down); past 90° the camera continues over the top and
looks up at the model from the far side; negative elevation looks up from
below. Roll tilts the horizon around the camera's own viewing direction.

Colors encode surface *orientation* only (from each triangle's normal) —
there's no lighting model, so two parallel faces at different positions or
heights are indistinguishable by color alone; use the bounding box in each
image, or `extract_geometry.py`, to tell them apart.

### Decoding view settings from filenames

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

### How surfaces are colored

Each triangle is colored by its normal direction (mapped from `[-1,1]` to
`[0,1]` per axis), so every triangle belonging to the same flat surface gets
the exact same color, and differently-oriented surfaces get different
colors. Hidden surfaces are removed with a real per-pixel z-buffer, so
occluded faces never show through nearer ones.

## extract_geometry.py

```
uv run extract_geometry.py <path-to.stl|.scad> [--decimals N] [--out file.json] [--min-area MM2]
```

- `--decimals` (default `4`): rounding used both to merge shared vertices
  into planes/loops and in the output numbers.
- `--out` (default: prints to stdout): write the JSON to a file instead.
- `--min-area` (optional): drop faces smaller than this area (mm²) from the
  output. The full output can be large (every hexagon in a honeycomb grid
  is its own hole loop) — for a quick look, list faces by area/normal first
  (e.g. with `jq` or Python) before reading loop data for one face, or use
  `--min-area` to cut small holes out of the output entirely.

If run on a `.scad` file, the exported STL is written as
`<scad filename>.scad-export.stl` in the current directory.

Output is a JSON object with the overall bounding box and a `faces` list
(largest area first). Each face groups every triangle sharing the same
plane (normal + offset) and reports:

- `normal`, `plane_offset`: the plane equation.
- `area`, `triangle_count`, `bounding_box`.
- `frame`: `{origin, u, v}` — the face's own 2D coordinate system, in model
  coordinates: a `points_2d` point `[p, q]` is `origin + u*p + v*q`. For an
  axis-aligned normal, `u`/`v` are world axes (e.g. a `+Z` face has
  `u=[1,0,0], v=[0,1,0]`, so `points_2d` is just `(x, y)`); otherwise they're
  an arbitrary orthonormal basis in the plane.
- `loops`: boundary polygons on that plane — one `"outer"` loop (the
  face's outline) plus any `"hole"` loops (e.g. each hexagon in a
  honeycomb grid, or a handle cutout), each as both `points_3d` (model
  coordinates) and `points_2d` (coordinates in `frame`, ready to drop into
  an OpenSCAD `polygon()`).

Planes are grouped with a tolerance (normal 2e-3, offset 0.02 mm), so one
flat face is one entry; a plane with several separate regions has several
`"outer"` loops (holes are decided by loop winding, not size). Zero-area
(degenerate) triangles are dropped before grouping, with a warning on
stderr. The top level also reports `volume`, `volume_sign` (-1 for an
inside-out mesh, which also swaps which loops are `"outer"` vs `"hole"`),
`watertight`, `bad_edge_count` (a non-zero count, or an "open boundary
loop" warning on stderr, means the mesh isn't watertight there).

### `--define` for `.scad` input

Both scripts take `--define NAME=VALUE` / `-D` (repeatable) for `.scad`
input, passed to `openscad -D` without a shell, e.g. `-D 'part="left"'`.

### `section` subcommand

```
extract_geometry.py section FILE --cut Y=40,160 [--scan Z=-13] [--profile X] [--prominence 0.2] [--shift X,Y,Z] [--out file.json]
```

For each `--cut` value, cuts the mesh with the plane `AXIS=value` and
reports:

- `axis`, `value` (nudged off a vertex by up to `1e-6` if needed, with a
  warning), `segment_count`, `bounding_box` of the cut.
- `scans` (if `--scan AXIS=V1,V2,...` given): for each scan line, `line`
  (the scan axis and value), `along` (the remaining axis), and `intervals`
  — `[start, end]` pairs where the mesh has material along that line.
- `profile` (if `--profile AXIS` given): `along`, `by` (the other in-plane
  axis), `prominence`, and `low`/`high` — the peaks and troughs of the
  low/high boundary (min/max of `along`) as a function of `by`. Each entry
  is `{"kind": "min"|"max", "value": ..., "at": pos}`, or
  `"at": [pos_from, pos_to]` when the peak/trough is a flat plateau rather
  than a single point.

### `slice` subcommand

```
extract_geometry.py slice FILE --origin X,Y,Z --normal NX,NY,NZ [--u UX,UY,UZ] [--window UMIN,UMAX,VMIN,VMAX [--within]] [--scad] [--svg cut.svg] [--out file.json]
```

Cuts the mesh with any plane and returns the cross-section as vector
geometry: closed loops (or, with a clipping window, open chains) of `line`,
`arc` and `circle` items, in the plane's own 2D coordinates.

- Frame: `--origin` becomes `(0, 0)`; `u` is `--u` projected onto the plane
  (default: the same frame `extract_geometry.py` uses for faces) and
  `v = normal × u`, so `(u, v, normal)` is right-handed. The output `frame`
  maps a point `(p, q)` back to `origin + u*p + v*q`. For a vertical cut,
  pick the normal's sign so `v` comes out as `+Z` (e.g. `--normal -1,-1,0
  --u 1,-1,0`); flipping the normal mirrors `v`.
- Winding: outer loops run counter-clockwise, holes clockwise, seen from the
  normal side — the same convention OpenSCAD's `polygon()` paths use.
- `--window`: clip the cut to a `(u, v)` box. Plan cuts of a plate are one
  huge outer loop, so clipping is how you zoom on a single feature; loops
  crossing the box edge become `"kind": "open"` chains. `--within` instead
  keeps only whole loops lying entirely inside the box, unclipped.
- Simplification: Douglas-Peucker with `--tol` (default 0.002 mm) removes
  collinear mesh vertices, then runs of at least `--min-arc-points` points
  that fit one circle within `--arc-tol` (default 0.005 mm) become an `arc`
  (`center`, `radius`, `start_deg`, `end_deg`, `ccw`); a whole loop on one
  circle is a single `circle`. Curves that aren't circular (e.g. a thread
  flank cut at an angle) come out as several short arcs or lines.
- `--scad` adds a ready `polygon(points, paths)` snippet (closed loops only,
  arcs expanded at `--fn` segments per turn); `--svg` writes a drawing
  (outer blue, holes red, open chains grey) — convert it with
  `rsvg-convert` if you need a PNG to look at.
- Vector options (`--origin`, `--normal`, `--u`, `--window`, `--shift`)
  accept negative values without the `=` form.

### `compare` subcommand

```
extract_geometry.py compare A B [--mirror-x] [--shift X,Y,Z] [--define-b NAME=VALUE] [--normal-tol T] [--cut ... --scan ...] [--out file.json]
```

`--mirror-x` and `--shift` apply to `B` only (mirror first) — there is no
automatic alignment. `--define` applies to both `A` and `B`; `--define-b`
adds an extra override for `B` only, so you can compare one `.scad` file
against itself under two different `-D` values (e.g.
`--define 'part="left"' --define-b 'part="right"'`). `--normal-tol`
(default `0.002`) is the per-component normal tolerance used to group `A`'s
and `B`'s faces together in `area_by_normal`.

Output:

- `a_stats`/`b_stats`: each mesh's `volume`, `volume_sign`, `watertight`,
  `bad_edge_count`, `bounding_box`.
- `volume_diff`: `b - a`.
- `area_by_normal`: one row per matched normal direction across both
  meshes, with `normal`, `area_a`, `area_b`, `diff` (`b - a`).
- `sections` (if `--cut` given): one entry per cut value, with `a`/`b`
  (each the full `section` output above) and `scan_diffs` — per scan line,
  `same_count` (whether `A` and `B` have the same number of intervals) and
  `b_minus_a` (interval-by-interval difference, `null` if counts differ).
