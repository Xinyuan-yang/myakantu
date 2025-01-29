L = 1;
s = L / 4;

Point(1) = {-L / 2., -L / 2., 0, s};

line1[] = Extrude{L / 2, 0, 0} { Point{1}; };
line2[] = Extrude{L / 2, 0, 0} { Point{2}; };

Physical Point("corner") = {1};

Physical Point("x-minus") = {1};
Physical Point("x-plus") = {3};

Physical Point("traction") = {2};
