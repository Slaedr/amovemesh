// flow past sphere

// actual boundary
srad = 1.0;
// for partition
prad = 1.25;
// far field
ffs = 15.0;
// mesh size at far field
hf = 2.0;
// mesh size at inner circle (not used)
hi = 0.2;
// mesh size at outer circle (not used)
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
Circle(1) = {2,1,3}; Transfinite Line {1} = n1;
Circle(2) = {3,1,4}; Transfinite Line {2} = n1;
Circle(3) = {4,1,5}; Transfinite Line {3} = n1;
Circle(4) = {5,1,2}; Transfinite Line {4} = n1;
Circle(5) = {3,1,6}; Transfinite Line {5} = n1;
Circle(6) = {6,1,5}; Transfinite Line {6} = n1;
Circle(7) = {5,1,7}; Transfinite Line {7} = n1;
Circle(8) = {7,1,3}; Transfinite Line {8} = n1;
Circle(9) = {2,1,6}; Transfinite Line {9} = n1;
Circle(10) = {6,1,4}; Transfinite Line {10} = n1;
Circle(11) = {4,1,7}; Transfinite Line {11} = n1;
Circle(12) = {7,1,2}; Transfinite Line {12} = n1;

// generate surfaces of sphere 1 (Line loops need to be ordered!)
Line Loop(50) = {1,5,-9}; Ruled Surface(1) = {50};
Line Loop(51) = {9,6,4}; Ruled Surface(2) = {51};
Line Loop(52) = {6,-3,-10}; Ruled Surface(3) = {52};
Line Loop(53) = {10,5,-2}; Ruled Surface(4) = {53};
Line Loop(54) = {1,-8,12}; Ruled Surface(5) = {54};
Line Loop(55) = {12,7,-4}; Ruled Surface(6) = {55};
Line Loop(56) = {7,-11,3}; Ruled Surface(7) = {56};
Line Loop(57) = {11,8,2}; Ruled Surface(8) = {57};
// outer surface of sphere 1
Surface Loop(80) = {1,2,3,4,5,6,7,8};

//sphere-2
Point(12) = { prad,   0,  0, 0.1};
Point(13) = {   0, prad,  0, 0.1};
Point(14) = {-prad,   0,  0, 0.1};
Point(15) = {   0,-prad,  0, 0.1};
Point(16) = {0, 0, prad, ho};
Point(17) = {0, 0, -prad, ho};
// arcs of sphere 2
Circle(21) = {12,1,13}; Transfinite Line {21} = n1;
Circle(22) = {13,1,14}; Transfinite Line {22} = n1;
Circle(23) = {14,1,15}; Transfinite Line {23} = n1;
Circle(24) = {15,1,12}; Transfinite Line {24} = n1;
Circle(25) = {13,1,16}; Transfinite Line {25} = n1;
Circle(26) = {16,1,15}; Transfinite Line {26} = n1;
Circle(27) = {15,1,17}; Transfinite Line {27} = n1;
Circle(28) = {17,1,13}; Transfinite Line {28} = n1;
Circle(29) = {12,1,16}; Transfinite Line {29} = n1;
Circle(30) = {16,1,14}; Transfinite Line {30} = n1;
Circle(31) = {14,1,17}; Transfinite Line {31} = n1;
Circle(32) = {17,1,12}; Transfinite Line {32} = n1;

// surfaces of sphere 2
Line Loop(60) = {21,25,-29}; Ruled Surface(11) = {60};
Line Loop(61) = {29,26,24}; Ruled Surface(12) = {61};
Line Loop(62) = {26,-23,-30}; Ruled Surface(13) = {62};
Line Loop(63) = {30,25,-22}; Ruled Surface(14) = {63};
Line Loop(64) = {21,-28,32}; Ruled Surface(15) = {64};
Line Loop(65) = {32,27,-24}; Ruled Surface(16) = {65};
Line Loop(66) = {27,-31,23}; Ruled Surface(17) = {66};
Line Loop(67) = {31,28,22}; Ruled Surface(18) = {67};
// outer surface of sphere 2
Surface Loop(90) = {11,12,13,14,15,16,17,18};

// lines joining inner circle (1) points to outer circle (2) points
Line(13) = {2,12}; Transfinite Line {13} = n2 Using Progression rc;
Line(14) = {3,13}; Transfinite Line {14} = n2 Using Progression rc;
Line(15) = {4,14}; Transfinite Line {15} = n2 Using Progression rc;
Line(16) = {5,15}; Transfinite Line {16} = n2 Using Progression rc;
Line(17) = {6,16}; Transfinite Line {17} = n2 Using Progression rc;
Line(18) = {7,17}; Transfinite Line {18} = n2 Using Progression rc;

//vol-1 (1st octant)
/*Line Loop(70)={14,25,-17,-5}; Plane Surface(20) = {70};
Transfinite Surface {20}; //Recombine Surface {20};
Line Loop(71) = {17,-29,-13,9}; Plane Surface(21) = {71};
Transfinite Surface {21}; //Recombine Surface {21};
Line Loop(72) = {1,14,-21,-13}; Plane Surface(22) = {72};
Transfinite Surface {22}; //Recombine Surface {22};
Surface Loop(100) = {20,21,22,11,1};
Volume(1000) = {100};

//vol-2 (2nd octant..)
Line Loop(73) = {14,22,-15,-2}; Ruled Surface(23) = {73};
Transfinite Surface {23}; //Recombine Surface{2};
Line Loop(74) = {15,-30,-17,10}; Ruled Surface(24) = {74};
Transfinite Surface {24};
Surface Loop(101) = {20,24,23,4,14};
Volume(1001) = {101};

//vol-3
Line Loop(75) = {3,16,-23,-15}; Plane Surface(25) = {75};
Transfinite Surface {25}; //
Line Loop(76) = {16,-26,-17,6}; Plane Surface(26) = {76};
Transfinite Surface {26};
Surface Loop(102) = {25,26,24,3,13};

//vol-4
Line Loop(77) = {4,13,-24,-16}; Plane Surface(27) = {77};
Transfinite Surface {27}; //Recombine Surface {27};*/

//cube-1
Point(31) = { -ffs, -ffs, -ffs, hf};
Point(32) = {  ffs, -ffs, -ffs, hf};
Point(33) = {  ffs,  ffs, -ffs, hf};
Point(34) = { -ffs,  ffs, -ffs, hf};
Point(35) = { -ffs, -ffs, ffs, hf};
Point(36) = {  ffs, -ffs, ffs, hf};
Point(37) = {  ffs,  ffs, ffs, hf};
Point(38) = { -ffs,  ffs, ffs, hf};
Line(31)  = {31,32};
Line(32)  = {32,33};
Line(33)  = {33,34};
Line(34)  = {34,31};
Line(35) = {35,36};
Line(36) = {36,37};
Line(37) = {37,38};
Line(38) = {38,35};
Line(39) = {31,35};
Line(40) = {32,36};
Line(41) = {33,37};
Line(42) = {34,38};

