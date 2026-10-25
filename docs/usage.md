# Output reference

Flags and defaults: `uv run <script>.py [subcommand] --help`. This file covers
what the outputs mean. STL input may be binary or ASCII.

`--brief` (faces, slice, compare) prints the same result as compact text, one
line per face, loop, item or diff, instead of JSON: `L from -> to` for a line,
`A from -> to c center r radius ccw|cw` for an arc, `C center r radius` for a
circle. Read it directly; reach for the JSON only for fields `--brief` leaves
out (face loop points, section profiles).

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
`collapsed_triangle_count`, `triangle_count`, `face_count`, and `faces`,
largest area first.

- Slivers thinner than the vertex rounding (`--decimals`) collapse to a line
  or point; they are counted in `collapsed_triangle_count` and left out of
  the edge check, which is common in CSG output and harmless.
- `bad_edges` (only when `bad_edge_count` > 0): where they are, grouped per
  kind on a 1 mm grid, biggest first: `kind` (`open`: one face, a hole in the
  surface; `over-shared`: three or more faces, e.g. two solids touching
  along an edge), `count`, `center`, `bounding_box`.
- `--box X0,X1,Y0,Y1,Z0,Z1` keeps only faces whose bounding box lies inside
  the box. `--min-area` and `--box` filter before the per-face loops are
  built, so a tight filter is also the fast path; many surviving faces are
  built on all cores.
- `--stats` prints only the top level above, without `faces`, and skips the
  plane grouping entirely.

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

## thread

Radial rays from the axis (`--center`, `--axis`, default `+Z`) at `--angles`
evenly spaced angles, every `--step` mm along it; the first surface each ray
meets gives `r(t)`. Angles turn counter-clockwise seen from the axis tip,
starting at `u` (`+X` for a `+Z` axis); `t` is measured along the axis from
`--center`.

- `minor_radius` (crest, the thread's inner tips) and `major_radius` (root),
  `depth`.
- `pitch` (axial distance between crests), `lead` (axial advance per turn),
  `starts` (`lead / pitch`), `handedness` (`right`: the crest climbs as the
  angle grows).
- `crest_flat`, `root_flat`, `flank_axial` (axial length of one flank),
  `flank_angle_deg` (flank against the radial direction: 45 for a square
  45° flank, 60 for a 30° half-angle).
- Phase: `crest_center_at_0deg` (the crest centre on the angle-0 ray),
  `crest_angle_at_ref` (the angle where a crest centre sits at `t = ref`),
  `repeats_every_deg`.
- `profile_0deg`: the simplified `[t, r]` corners of the angle-0 ray;
  `crest_centers_by_angle`: whole crests seen on each ray.
- Entry chamfers are recognised (outside the root, or sloped unlike the
  thread's flanks) and skipped, so the default `--range` is the whole hole.
  Flats are only measured between two thread runs, so a thread with less
  than about two pitches of length reports `null` pitch or phase: lengthen
  `--range` or add `--angles`.
- The bore need not be round (an octagon, a slot through it): crest and root
  radii are found per ray, rays that see no thread (a constant radius, or
  nothing) are left out, and the flats come from the rays whose crest sits at
  `minor_radius`. A crest split by another feature (a notch) is dropped as off
  the helix rather than bending `lead`.

## compare

`--scale SX,SY,SZ` then `--shift` apply to B only; there is no auto-alignment.
A negative scale mirrors (`-1,1,1`: x -> -x) and two make a 180° turn
(`-1,-1,1`: about Z), so a part's mirrored or turned variant lines up with it.
`--define` applies to both models, `--define-b` adds an override for B only
(compare one `.scad` under two `-D` values).

- `a_stats`, `b_stats`: as the faces command's top level (`volume`,
  `volume_sign`, `watertight`, `bad_edge_count`, `collapsed_triangle_count`,
  `bad_edges`, `triangle_count`) plus `bounding_box`.
- `volume_diff`: B − A.
- `area_by_normal`: rows of `normal`, `area_a`, `area_b`, `diff` (B − A),
  normals grouped within `--normal-tol`.
- `sections` (with `--cut`): per cut, `a`/`b` as in `section`, plus
  `scan_diffs`: per scan line `same_count` and `b_minus_a` (interval-wise,
  `null` when counts differ).
- `slices` (with `--slice AXIS=…`): both meshes cut as in `slice`, framed so
  `(u, v)` is `(y, z)` for X, `(x, z)` for Y and `(x, y)` for Z; `--window`
  clips in that frame. Per cut: `net_area_a`/`_b`, `area_diff`,
  `loop_count_a`/`_b`; `unmatched_a`/`_b` (loops with no same-kind partner
  overlapping by more than half its bounding box); and `loop_diffs` for
  matched loops whose areas differ by more than `--area-tol` or whose items
  differ, largest first: `area_a`, `area_b`, `diff`, and `only_a`/`only_b`,
  the lines/arcs with no partner within `--match-tol` (ends, centre and
  radius). Arcs fitted along a spiral edge (a thread) differ run to run, so
  expect a few matched-area loops listing arcs there.
