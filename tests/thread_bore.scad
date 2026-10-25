// A block with a threaded bore along +Z through (0, 0): known answers for
// `thread` (pitch 3.125, minor 6.75, major 7.15, groove centre at z 5.44375 on +X).
use <../docs/helix_groove.scad>

oct = false;     // octagonal bore with a slot out to +Y, instead of a round one
left = false;    // left-hand thread
threaded = true; // false: a plain bore

$fn = 180;

difference() {
    translate([-10, -10, 0]) cube([20, 20, 12]);
    translate([0, 0, -1]) {
        if (oct) linear_extrude(14) intersection() {
            square(13.5, center = true);
            rotate(45) square(13.5, center = true);
        }
        else cylinder(r = 6.75, h = 14);
        if (oct) translate([-4, 0, 0]) cube([8, 11, 14]);
    }
    if (threaded)
        helix_groove(pitch = 3.125, r_in = 6.4, r_out = 7.15, root = 0.625, flank = 0.41667 / 0.4,
                     z0 = 5.44375, z_from = -1, z_to = 13, left = left);
}
