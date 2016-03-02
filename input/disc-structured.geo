srad = 1.0;
partrad = 3.0;
ffs = 15.0;
num_points = 9;

hs = 0.25;
hf = 1.5;
Point(1) = {0, 0, 0, hs};
Point(2) = {srad, 0, 0, hs};
Point(3) = {0, srad, 0, hs};
Point(4) = {-srad, 0, 0, hs};
Point(5) = {0, -srad, 0, hs};

Point(8) = {ffs, ffs, 0, hf};
Point(9) = {-ffs, ffs, 0, hf};
Point(10) = {-ffs, -ffs, 0, hf};
Point(11) = {ffs, -ffs, 0, hf};

Circle(1) = {2, 1, 3}; Transfinite Line {1} = num_points;
Circle(2) = {3, 1, 4}; Transfinite Line {2} = num_points;
Circle(3) = {4, 1, 5}; Transfinite Line {3} = num_points;
Circle(4) = {5, 1, 2}; Transfinite Line {4} = num_points;

Line Loop(15) = {1,2,3,4};

Line(32) = {8, 9};
Line(33) = {9, 10};
Line(34) = {10, 11};
Line(35) = {11, 8};

Line Loop(44) = {32,33,34,35};
Ruled Surface(45) = {44,15};
//Transfinite Surface {45};
Recombine Surface(45);
Physical Surface(46) = {45};