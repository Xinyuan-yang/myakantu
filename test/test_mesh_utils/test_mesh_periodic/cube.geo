L = 1;
s = 0.5;

Point(1) = {-L/2., -L/2., 0, s};

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

Physical Volume("body") = {volume[1]};
