L = 1;
s = .3;

Point(1) = {-L / 2., -L / 2., 0, s};

line[] = Extrude{L, 0, 0} { Point{1}; };

surface[] = Extrude{0, L, 0} { Line{line[1]}; };

Transfinite Surface "*";
Recombine Surface "*";
//+

Physical Point("corner", 6) = {1};
Physical Curve("x-minus", 7) = {3};
Physical Curve("x-plus", 8) = {4};

Physical Curve("y-minus", 9) = {1};
Physical Curve("y-plus", 10) = {2};

Physical Surface("bulk", 9) = {5};
