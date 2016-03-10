// flow past sphere

// actual boundary
srad = 1.0;
// for partition
prad = 1.25;
// far field
ffs = 15.0;
// mesh size at far field
hf = 2.0;
// mesh size at inner circle
hi = 0.2;
// mesh size at outer circle
ho = 0.5;

// number of points along the inner circle
n1 = 13;
n2 = 50;
// ratio (r) for geometric progression of transfinite points
rc = 1.2;

//center
Point(1) = {   0,   0,  0, 0.1};

//sphere-1
Point(2) = { srad,   0,  0, 0.1};
Point(3) = {   0, srad,  0, 0.1};
Point(4) = {-srad,   0,  0, 0.1};
Point(5) = {   0,-srad,  0, 0.1};
Point(6) = {0, 0, srad, hi};
Point(7) = {0, 0, -srad, hi};

// Generate arcs of sphere 1
Circle(1) = {2,1,3}; 
Circle(2) = {3,1,4}; 
Circle(3) = {4,1,5}; 
Circle(4) = {5,1,2}; 
Circle(5) = {3,1,6}; 
Circle(6) = {6,1,5}; 
Circle(7) = {5,1,7}; 
Circle(8) = {7,1,3}; 
Circle(9) = {2,1,6}; 
Circle(10) = {6,1,4};
Circle(11) = {4,1,7};
Circle(12) = {7,1,2};

// generate surfaces of sphere 1 (Line loops need to be ordered!)
Line Loop(50) = {1,5,-9}; Ruled Surface(1) = {50};
Line Loop(51) = {9,6,4}; Ruled Surface(2) = {51};
Line Loop(52) = {6,-3,-10}; Ruled Surface(3) = {52};
Line Loop(53) = {10,5,-2}; Ruled Surface(4) = {53}; 
Line Loop(54) = {1,-8,12}; Ruled Surface(5) = {54}; 
Line Loop(55) = {12,7,-4}; Ruled Surface(6) = {55}; 
Line Loop(56) = {7,-11,3}; Ruled Surface(7) = {56}; 
Line Loop(57) = {11,8,2}; Ruled Surface(8) = {57};
// surface of sphere 1
Surface Loop(80) = {1,2,3,4,5,6,7,8};

//cube- farfield
Point(31) = { -ffs, -ffs, -ffs, hf};
Point(32) = {  ffs, -ffs, -ffs, hf};
Point(33) = {  ffs,  ffs, -ffs, hf};
Point(34) = { -ffs,  ffs, -ffs, hf};
Point(35) = { -ffs, -ffs, ffs, hf};
Point(36) = {  ffs, -ffs, ffs, hf};
Point(37) = {  ffs,  ffs, ffs, hf};
Point(38) = { -ffs,  ffs, ffs, hf};
Line(41)  = {31,32};
Line(42)  = {32,33};
Line(43)  = {33,34};
Line(44)  = {34,31};
Line(45) = {35,36};
Line(46) = {36,37};
Line(47) = {37,38};
Line(48) = {38,35};
Line(49) = {31,35};
Line(50) = {32,36};
Line(51) = {33,37};
Line(52) = {34,38};

// Far-field volume
Line Loop(131) = {44, 49, -48, -52};
Plane Surface(132) = {131};
Line Loop(133) = {47, 48, 45, 46};
Plane Surface(134) = {133};
Line Loop(135) = {46, -51, -42, 50};
Plane Surface(136) = {135};
Line Loop(137) = {43, 44, 41, 42};
Plane Surface(138) = {137};
Line Loop(139) = {52, -47, -51, 43};
Plane Surface(140) = {139};
Line Loop(141) = {45, -50, -41, 49};
Plane Surface(142) = {141};
Surface Loop(143) = {132, 138, 140, 134, 142, 136};
Volume(144) = {80, 143};

// Interior sphere boundary
Physical Surface(2) = {4, 1, 5, 8, 7, 3, 2, 6};
Physical Volume(1) = {144};
