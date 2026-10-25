// A thread groove cut into a bore wall, along +Z: subtract it from the part.
// Single start; right hand climbs as the angle grows (counter-clockwise seen
// from +Z). The arguments map onto `extract_geometry.py thread` output:
//   pitch   lead
//   r_in    a little below minor_radius, so the flanks start in the bore's air
//   r_out   major_radius
//   root    root_flat: the groove's width at r_out
//   flank   flank_axial / depth: axial run per mm of radius (1 = 45 deg flank)
//   z0      the groove's centre on the +X ray: crest_center_at_0deg + pitch / 2
//   z_from, z_to  where the groove is cut
// Rotate or multmatrix the result onto the bore's axis; a pure rotation keeps
// the polyhedron's winding, a mirror is fine too (OpenSCAD fixes it).
module helix_groove(pitch, r_in, r_out, root, flank, z0, z_from, z_to, left = false, n = 120) {
    if (left) mirror([0, 1, 0]) helix_groove(pitch, r_in, r_out, root, flank, z0, z_from, z_to, false, n);
    else {
        hw_out = root / 2;
        hw_in = hw_out + flank * (r_out - r_in);
        turns0 = floor((z_from - z0) / pitch) - 1;
        turns1 = ceil((z_to - z0) / pitch) + 1;
        steps = (turns1 - turns0) * n;
        prof = [[r_in, -hw_in], [r_out, -hw_out], [r_out, hw_out], [r_in, hw_in]];
        pts = [for (i = [0 : steps]) let(a = 360 * i / n, z = z0 + pitch * (turns0 + i / n))
                  for (p = prof) [p[0] * cos(a), p[0] * sin(a), z + p[1]]];
        // polyhedron() wants every face clockwise seen from outside; the other way
        // round renders inside out (volume off by the groove's, not watertight)
        faces = concat(
            [[3, 2, 1, 0]],
            [for (i = [0 : steps - 1]) for (j = [0 : 3])
                let(a = 4 * i + j, b = 4 * i + (j + 1) % 4) [b, b + 4, a + 4, a]],
            [[4 * steps, 4 * steps + 1, 4 * steps + 2, 4 * steps + 3]]);
        intersection() {
            polyhedron(pts, faces);
            translate([-r_out - 1, -r_out - 1, z_from]) cube([2 * r_out + 2, 2 * r_out + 2, z_to - z_from]);
        }
    }
}
