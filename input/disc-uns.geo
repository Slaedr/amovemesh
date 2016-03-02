srad = 1.0;
ffs = 10.0;
num_points = 9;

hs = 0.25;
hf = 1.0;
Point(1) = {0, 0, 0, hs};
Point(2) = {srad, 0, 0, hs};
Point(3) = {0, srad, 0, hs};
Point(4) = {-srad, 0, 0, hs};
Point(5) = {0, -srad, 0, hs};

Point(8) = {ffs, ffs, 0, hf};
Point(9) = {-ffs, ffs, 0, hf};
Point(10) = {-ffs, -ffs, 0, hf};
Point(11) = {ffs, -ffs, 0, hf};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2}; 

Line Loop(15) = {1,2,3,4};

Line(32) = {8, 9};
Line(33) = {9, 10};
Line(34) = {10, 11};
Line(35) = {11, 8};

Line Loop(44) = {32,33,34,35};
Plane Surface(45) = {44,15};
//Recombine Surface(45);
Physical Surface(46) = {45};

// specify boundary layer size field
Field[1] = BoundaryLayer;
Field[1].EdgesList = {1,2,3,4};
//Field[1].NodesList = {2, 3, 4, 5};
//Field[1].FacesList = {45};
Field[1].Quads = 0;
Field[1].hwall_n = 0.01;
Field[1].hwall_t = 0.25;
Field[1].thickness = 0.5;
Field[1].IntersectMetrics = 1;
Background Field = 1;
