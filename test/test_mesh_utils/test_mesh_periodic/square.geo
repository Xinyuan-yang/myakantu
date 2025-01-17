L = 1;
s = 2;

Point(1) = {-L / 2., -L / 2., 0, s};

line[] = Extrude{L, 0, 0} { Point{1}; };

surface[] = Extrude{0, L, 0} { Line{line[1]}; };

Physical Point("corner", 6) = {1};
Physical Curve("bottom", 7) = {1};
Physical Curve("top", 8) = {2};

Transfinite Surface "*";
Recombine Surface "*";
//+
Physical Surface("bulk", 9) = {5};
