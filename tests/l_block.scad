// An L-shaped block with a round hole: no symmetry, so a wrong transform shows.
turn = false;    // turned 180 deg about Z
mirrored = false; // mirrored in X

module body() difference() {
    union() {
        cube([20, 10, 5]);
        cube([5, 20, 5]);
    }
    translate([12, 5, -1]) cylinder(r = 2, h = 7, $fn = 64);
}

rotate([0, 0, turn ? 180 : 0])
    if (mirrored) mirror([1, 0, 0]) body();
    else body();
