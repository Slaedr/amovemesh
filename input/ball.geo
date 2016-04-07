// actual boundary
srad = 1.0;

// mesh size at inner circle
//hi = 0.5;
//hi = 0.25;
//hi = 0.125;
hi = 0.0625;

//center
Point(1) = {   0,   0,  0, 0.1};

//sphere-1
Point(2) = { srad,   0,  0, hi};
Point(3) = {   0, srad,  0, hi};
Point(4) = {-srad,   0,  0, hi};
Point(5) = {   0,-srad,  0, hi};
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
Volume(100) = {80};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};

Physical Volume(10) = {100};
