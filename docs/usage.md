# Output reference

Flags and defaults: `uv run <script>.py [subcommand] --help`. This file covers
what the outputs mean.

## render_views.py

Writes one PNG per angle into `--out` (default `.stl-renders-<name>` next to
the model), clearing old `*.png` there first. A `.scad` input's exported STL
lands in `--out` too.

Each image sweeps one camera rotation while the other two stay at their base
value (`elev=25`, `azim=-60`, `roll=0`):

| Filename | Camera |
| --- | --- |
| `elevNNN.png` | elevation `NNN`, azim -60, roll 0 |
| `azimNNN.png` | azimuth `NNN`, elev 25, roll 0 |
| `rollNNN.png` | roll `NNN`, elev 25, azim -60 |

- Azimuth 0 / elevation 0: camera on +X looking at the origin (you see the
  model's +X side, with +Y to the right). Azimuth turns the camera about Z
  (90 = on +Y, 180 = on -X).
- Elevation 90: straight down, +X right, +Y up. Past 90 the camera goes over
  the top; negative elevation looks up from below.
- Roll tilts the horizon around the viewing direction.
- Color = triangle normal mapped from `[-1,1]` to `[0,1]` per RGB channel; no
  lighting. Hidden surfaces are removed with a z-buffer.

## extract_geometry.py (faces)

Top level: `bounding_box`, `volume`, `volume_sign` (-1 = inside-out mesh,
which also swaps outer/hole), `watertight`, `bad_edge_count`,
`triangle_count`, `face_count`, and `faces`, largest area first.

A face is every triangle on one plane (normal within 2e-3, offset within
0.02 mm), so one plane with separate regions has several `outer` loops.
Zero-area triangles are dropped with a stderr warning. Each face has:

- `id`, `normal`, `plane_offset` (plane: `normal · p = offset`), `area`,
  `triangle_count`, `bounding_box`.
- `frame`: `{origin, u, v}`; a 2D point `[p, q]` is `origin + u*p + v*q`.
  `origin` is a vertex of the face, not the world origin. Axis-aligned normals
  use world axes for `u`/`v` (`+Z`: X, Y; `-Z`: X, -Y; `+X`: Y, Z; `+Y`: Z, X);
  other normals get an arbitrary orthonormal pair.
- `loops`: `outer` outlines and `hole`s (decided by winding), each as
  `points_3d` (world) and `points_2d` (in `frame`).

## section

One entry per `--cut` in `sections`, with the cut value nudged by up to 1e-6
off a vertex (stderr warning):

- `axis`, `value`, `segment_count`, `bounding_box`.
- `scans` (with `--scan AXIS=…`): per line, `line`, `along` (the remaining
  axis) and `intervals`: `[start, end]` stretches of material.
- `profile` (with `--profile AXIS`): `along`, `by`, `prominence`, and
  `low`/`high`: peaks and troughs of the cut's min/max boundary along `along`
  as a function of `by`, each `{"kind": "min"|"max", "value", "at"}`, with
  `"at": [from, to]` for a flat plateau.

## slice

`frame`, `loop_count`, `loops` (and `scad` with `--scad`).

- Frame: `--origin` is `(0, 0)`; `u` is `--u` projected onto the plane
  (default: the face frame above); `v = normal × u`. For a vertical cut, pick
  the normal's sign so `v` comes out `+Z` (e.g. `--normal -1,-1,0 --u 1,-1,0`);
  flipping the normal mirrors `v`.
- Each loop has a `kind` (`outer` CCW, `hole` CW seen from the normal side,
  matching `polygon()`; `open` for a chain clipped by `--window`) and `items`:
  `line` (`from`, `to`), `arc` (`from`, `to`, `center`, `radius`, `start_deg`,
  `end_deg`, `ccw`) or `circle` (`center`, `radius`); each item has a `type`.
- Simplification: collinear points within `--tol` are dropped; runs of at
  least `--min-arc-points` points fitting one circle within `--arc-tol` become
  an arc. Non-circular curves (a thread flank, an oblique cut) come out as
  several short arcs or lines.
- `--window` clips to a `(u, v)` box; `--within` keeps whole loops inside it
  instead. `--svg` draws outer blue, holes red, open chains grey.

## compare

`--mirror-x` then `--shift` apply to B only; there is no auto-alignment.
`--define` applies to both models, `--define-b` adds an override for B only
(compare one `.scad` under two `-D` values).

- `a_stats`, `b_stats`: `volume`, `volume_sign`, `watertight`,
  `bad_edge_count`, `triangle_count`, `bounding_box`.
- `volume_diff`: B − A.
- `area_by_normal`: rows of `normal`, `area_a`, `area_b`, `diff` (B − A),
  normals grouped within `--normal-tol`.
- `sections` (with `--cut`): per cut, `a`/`b` as in `section`, plus
  `scan_diffs`: per scan line `same_count` and `b_minus_a` (interval-wise,
  `null` when counts differ).
