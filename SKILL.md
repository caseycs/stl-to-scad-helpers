---
name: stl-to-scad-helpers
description: Replicate an STL as parametric OpenSCAD, measure a mesh, or check a .scad against a reference. Use when recreating a 3D model in OpenSCAD, reading exact dimensions, holes or profiles off an STL/.scad, comparing two models, or viewing a model from many angles.
---

Exact numbers come from `extract_geometry.py`; pictures from `render_views.py`
only show structure (a vision model can't read dimensions off a PNG). Every
number you write into OpenSCAD traces back to tool output.

Run the scripts as `uv run <this skill's absolute dir>/<script>.py`. Work from
a scratch directory: a `.scad` input is exported to
`<name>.scad-export.stl` in the current directory.

## Tools

| Command | Gives you |
| --- | --- |
| `render_views.py MODEL --step 90 [--focus=X0,X1,Y0,Y1,Z0,Z1]` | PNGs from many angles, one color per surface orientation |
| `extract_geometry.py MODEL --min-area A --out geo.json` | health (`watertight`, `volume`, `bounding_box`) + every flat face: normal, offset, outline and hole loops |
| `extract_geometry.py slice MODEL --origin P --normal N [--window U0,U1,V0,V1 [--within]] [--scad]` | cross-section on any plane as lines/arcs/circles (centers, radii), optional `polygon()` |
| `extract_geometry.py section MODEL --cut Y=a,b --scan Z=c,d [--profile X]` | material intervals along scan lines; peaks/troughs of a cut's boundary |
| `extract_geometry.py compare A B [--cut … --scan …]` | volume, bbox and area-per-normal diffs (B − A), section diffs |

Flags: `--help` on each command. Output fields, frames and camera angles:
[`docs/usage.md`](docs/usage.md).

## Replicating an STL

1. **Survey.** `render_views.py MODEL --step 90`; view `elev090` (top, +X right,
   +Y up), `azim000`, `azim090`, `elev000`. Done when you can name the base
   body, every feature type (through-hole, pocket, boss, chamfer, thread,
   slot) and every repeated pattern.

2. **Measure the envelope.** `extract_geometry.py MODEL --min-area A --out geo.json`.
   Pick `A` so the file holds tens to a few hundred faces: curved surfaces
   split into thousands of tiny faces (a 224k-triangle plate gave 146k faces
   and 210 MB at the default). Read it through a summary, never whole:

   ```
   python3 -c '
   import json, sys
   g = json.load(open(sys.argv[1]))
   print({k: v for k, v in g.items() if k != "faces"})
   for f in g["faces"][:40]:
       k = [l["kind"] for l in f["loops"]]
       print(f["id"], f["area"], f["normal"], f["plane_offset"], f["bounding_box"], "outer", k.count("outer"), "holes", k.count("hole"))
   ' geo.json
   ```

   Done when you have the bounding box, a watertight mesh with
   `volume_sign` 1 (otherwise note where it's broken), and the distinct plane
   offsets per axis (the step heights).

3. **Find patterns.** Hole-loop centers of the large faces give counts and
   pitch:

   ```
   python3 -c '
   import json, sys
   f = next(f for f in json.load(open(sys.argv[1]))["faces"] if f["id"] == int(sys.argv[2]))
   for l in f["loops"]:
       p = l["points_3d"]; n = len(p)
       print(l["kind"], n, [round(sum(c[i] for c in p) / n, 3) for i in range(3)])
   ' geo.json FACE_ID
   ```

   Done when each repeated feature has a count, a pitch, and an origin.

4. **Plan the CSG.** List base body + features, each as an OpenSCAD
   primitive and the tool that will measure it:
   - axis-aligned outline or straight-edged hole → face loops (`points_3d`)
   - round, threaded or sloped feature → `slice` at mid-height of the
     feature, `--window` around one instance (`--within` for whole loops);
     arcs give center and radius, `--scad` gives a `polygon()`
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
   (leftover diffs on curved normals come from `$fn` tessellation). Finish
   with `render_views.py model.scad --step 90` against step 1's views.

For a single measurement or a diff, run just the matching step's tool.

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
