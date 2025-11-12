L = 1;
s = 0.5;

Point(1) = {0, 0, 0, s};

line[] = Extrude{L, 0, 0} { Point{1}; };

surface[] = Extrude{0, L, 0} { Line{line[1]}; };

volume[] = Extrude{0, 0, L} { Surface{surface[1]}; };

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Point("corner") = {1};

Physical Curve("line") = {1};

Physical Surface("x-minus") = {26};
Physical Surface("x-plus") = {18};

Physical Surface("y-minus") = {14};
Physical Surface("y-plus") = {22};

Physical Surface("z-minus") = {5};
Physical Surface("z-plus") = {27};

Physical Volume("bulk") = {1};
