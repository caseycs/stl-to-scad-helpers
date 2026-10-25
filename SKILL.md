---
name: stl-to-scad-helpers
description: Replicate an STL as parametric OpenSCAD, measure a mesh, or check a .scad against a reference. Use when recreating a 3D model in OpenSCAD, reading exact dimensions, holes or profiles off an STL/.scad, comparing two models, or viewing a model from many angles.
---

Exact numbers come from `extract_geometry.py`; pictures from `render_views.py`
only show structure (a vision model can't read dimensions off a PNG). Every
number you write into OpenSCAD traces back to tool output.

Run the scripts as `uv run <this skill's absolute dir>/<script>.py`. Work from
a scratch directory: a `.scad` input is exported to
`<name>.scad-export.stl` in the current directory. Pass `--brief` to `faces`,
`slice` and `compare` and read their text output directly, instead of
writing a script to summarise the JSON.

## Tools

| Command | Gives you |
| --- | --- |
| `render_views.py MODEL --step 90 [--focus=X0,X1,Y0,Y1,Z0,Z1]` | PNGs from many angles, one color per surface orientation |
| `extract_geometry.py MODEL --min-area A [--box X0,X1,Y0,Y1,Z0,Z1] --brief` | health (`watertight`, `volume`, `bounding_box`, where any `bad_edges` are) + every flat face: normal, offset, bbox, hole centres; `--box` keeps one feature's faces; JSON (no `--brief`) adds the loops' points |
| `extract_geometry.py MODEL --stats` | health only, no faces: seconds even on a million-triangle mesh |
| `extract_geometry.py slice MODEL --origin P --normal N [--window U0,U1,V0,V1 [--within]] [--scad]` | cross-section on any plane as lines/arcs/circles (centers, radii), optional `polygon()` |
| `extract_geometry.py thread MODEL --center P [--axis A] [--rmax R]` | a threaded hole: minor/major radius, pitch, lead, starts, hand, crest/root flats, flank angle, phase |
| `extract_geometry.py section MODEL --cut Y=a,b --scan Z=c,d [--profile X]` | material intervals along scan lines; peaks/troughs of a cut's boundary |
| `extract_geometry.py compare A B [--scale SX,SY,SZ] [--cut … --scan …] [--slice Z=a,b [--window …]]` | volume, bbox and area-per-normal diffs (B − A), section diffs; `--slice` diffs cuts loop by loop down to the lines/arcs only one side has; `--scale` mirrors or turns B first |

Flags: `--help` on each command. Output fields, frames and camera angles:
[`docs/usage.md`](docs/usage.md).

## Replicating an STL

1. **Survey.** `render_views.py MODEL --step 90`; view `elev090` (top, +X right,
   +Y up), `azim000`, `azim090`, `elev000`. Done when you can name the base
   body, every feature type (through-hole, pocket, boss, chamfer, thread,
   slot) and every repeated pattern.

2. **Measure the envelope.** `extract_geometry.py MODEL --min-area A --brief`.
   Pick `A` so it lists tens to a few hundred faces: curved surfaces split
   into thousands of tiny faces (a 224k-triangle plate gave 146k faces at the
   default). Done when you have the bounding box, a watertight mesh with
   `sign 1` (otherwise note where it's broken), and the distinct plane offsets
   per axis (the step heights).

3. **Find patterns.** The `hole centres` under each large face give counts
   and pitch. Done when each repeated feature has a count, a pitch, and an
   origin.

4. **Plan the CSG.** List base body + features, each as an OpenSCAD
   primitive and the tool that will measure it:
   - axis-aligned outline or straight-edged hole → face loops (`points_3d`
     in the JSON, `--box` around the feature)
   - round or sloped feature → `slice` at mid-height of the
     feature, `--window` around one instance (`--within` for whole loops);
     arcs give center and radius, `--scad` gives a `polygon()`
   - threaded hole → `thread` on its axis; check its phase on a second
     instance (`crest_angle_at_ref` at the same `--ref`); cut it with
     [`docs/helix_groove.scad`](docs/helix_groove.scad), whose arguments map
     onto `thread`'s output
   - heights, wall thickness, chamfers, depths → `section --scan` intervals
     and `--profile`

   Done when every feature from step 1 has a primitive and a measurement source.

5. **Measure.** One instance of each repeated feature in full, then its
   center on two more to confirm the pitch. Done when every dimension in the
   plan holds a number from tool output.

6. **Write the SCAD** in the STL's own coordinates (no recentering, so
   `compare` lines up without `--shift`), measured values as named variables,
   patterns as `for` loops.

7. **Verify.** `compare MODEL model.scad --cut … --scan …` with cuts through
   the features. Fix the largest `area_by_normal` |diff| first: its normal
   names the wrong surface family. Repeat until the bounding boxes match,
   `volume_diff` is within ~0.5 % of volume, every flat-face normal matches
   in area, and each scan has `same_count: true` with `b_minus_a` near 0
   (leftover diffs on curved normals come from `$fn` tessellation). Then
   `compare … --slice Z=…` at each feature's heights: a feature with the right
   volume in the wrong place (a transposed pattern, a mirrored mark) passes
   every check above and shows up only here, as `unmatched_*` loops or
   `only_a`/`only_b` items. Done when every loop matches and the remaining
   `loop_diffs` are arc-fit noise (same ends, slightly different arcs along a
   thread or other spiral edge). Finish with `render_views.py model.scad
   --step 90` against step 1's views.

For a single measurement or a diff, run just the matching step's tool.

## A variant of a replicated part

When a second STL looks like a mirrored, turned or slightly changed copy of
one already replicated, diff it against the first before measuring anything:

1. Find the transform: `compare FIRST VARIANT --scale … --shift … --brief`
   with the bounding boxes lined up (`--scale -1,1,1` mirrors in X,
   `-1,-1,1` turns 180° about Z). Done when the `area_by_normal` diffs are
   down to the changed features.
2. Locate the changes: `--slice` on each axis through the part; the
   `only_a`/`only_b` loops and items are the differences, the rest is
   shared.
3. Measure only those, add them to the first SCAD behind a `part = N`
   switch, and verify each part against its own STL as in step 7.

## Gotchas

- `points_2d` are relative to the face's `frame.origin` (a vertex of the face),
  and a `-Z` face's `v` is world `-Y`. Use `points_3d` for world positions.
- `slice` with `--origin 0,0,Z --normal 0,0,1` reports `(u, v)` = world `(x, y)`.
- `section`/`slice` print vertex-nudge warnings on stderr: redirect `2>/dev/null`
  before parsing stdout as JSON.
- `--focus` keeps whole triangles by centroid: pad the box past the feature.
  Negative values need the `=` form (`--focus=-40,40,…`).
- Render colors encode orientation only: parallel faces at different heights
  share a color.
- `over-shared` bad edges in your own SCAD's export: two solids touch along an
  edge or at a vertex. Overlap them (extend one into the other by ~1 mm).
- A hand-written `polyhedron()` wound the wrong way is inside out:
  `volume_diff` jumps by about the feature's own volume and the mesh is not
  watertight. Reverse each face's vertex order.
