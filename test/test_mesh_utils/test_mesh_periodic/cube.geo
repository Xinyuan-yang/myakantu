L = 1;
s = 0.5;

Point(1) = {0, 0, 0, s};

line[] = Extrude{L, 0, 0} {
  Point{1};
};

surface[] = Extrude{0, L, 0} {
  Line{line[1]};
};

volume[] = Extrude{0, 0, L} {
  Surface{surface[1]};
};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Point("corner") = {1};
Physical Line("line") = {1};
Physical Surface("bottom") = {5};
Physical Surface("top") = {27};
Physical Volume("body") = {1};
