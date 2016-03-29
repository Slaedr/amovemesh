// actual boundary
srad = 1.0;

// mesh size at inner circle
hi = 0.5;
// number of points along curves
np = 5;

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
Circle(1) = {2,1,3}; Transfinite Line {1} = np;
Circle(2) = {3,1,4}; Transfinite Line {2} = np;
Circle(3) = {4,1,5}; Transfinite Line {3} = np;
Circle(4) = {5,1,2}; Transfinite Line {4} = np;
Circle(5) = {3,1,6}; Transfinite Line {5} = np;
Circle(6) = {6,1,5}; Transfinite Line {6} = np;
Circle(7) = {5,1,7}; Transfinite Line {7} = np;
Circle(8) = {7,1,3}; Transfinite Line {8} = np;
Circle(9) = {2,1,6}; Transfinite Line {9} = np;
Circle(10) = {6,1,4}; Transfinite Line {10} = np;
Circle(11) = {4,1,7}; Transfinite Line {11} = np;
Circle(12) = {7,1,2}; Transfinite Line {12} = np;

// generate surfaces of sphere 1 (Line loops need to be ordered!)
Line Loop(50) = {1,5,-9}; Ruled Surface(1) = {50}; Transfinite Surface {1};
Line Loop(51) = {9,6,4}; Ruled Surface(2) = {51}; Transfinite Surface {2};
Line Loop(52) = {6,-3,-10}; Ruled Surface(3) = {52}; Transfinite Surface {3};
Line Loop(53) = {-10,-5,2}; Ruled Surface(4) = {53}; Transfinite Surface {4};
Line Loop(54) = {1,-8,12}; Ruled Surface(5) = {54}; Transfinite Surface {5};
Line Loop(55) = {12,7,-4}; Ruled Surface(6) = {55}; Transfinite Surface {6};
Line Loop(56) = {7,-11,3}; Ruled Surface(7) = {56}; Transfinite Surface {7};
Line Loop(57) = {11,8,2}; Ruled Surface(8) = {57}; Transfinite Surface {8};
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
