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
n1 = 20;
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
Line Loop(50) = {1,5,-9}; Ruled Surface(1) = {50}; Transfinite Surface {1};  //Recombine Surface {1};
Line Loop(51) = {9,6,4}; Ruled Surface(2) = {51}; Transfinite Surface {2}; //Recombine Surface {2};
Line Loop(52) = {6,-3,-10}; Ruled Surface(3) = {52}; Transfinite Surface {3}; //Recombine Surface {3};
Line Loop(53) = {10,5,-2}; Ruled Surface(4) = {53}; Transfinite Surface {4}; //Recombine Surface {4};
Line Loop(54) = {1,-8,12}; Ruled Surface(5) = {54}; Transfinite Surface {5}; //Recombine Surface {5};
Line Loop(55) = {12,7,-4}; Ruled Surface(6) = {55}; Transfinite Surface {6}; //Recombine Surface {6};
Line Loop(56) = {7,-11,3}; Ruled Surface(7) = {56}; Transfinite Surface {7}; //Recombine Surface {7};
Line Loop(57) = {11,8,2}; Ruled Surface(8) = {57}; Transfinite Surface {8}; //Recombine Surface {8};
// surface of sphere 1
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
Line Loop(60) = {21,25,-29}; Ruled Surface(11) = {60};  Transfinite Surface {11}; //Recombine Surface {11};
Line Loop(61) = {29,26,24}; Ruled Surface(12) = {61}; Transfinite Surface {12}; //Recombine Surface {12};
Line Loop(62) = {26,-23,-30}; Ruled Surface(13) = {62}; Transfinite Surface {13}; //Recombine Surface {13};
Line Loop(63) = {30,25,-22}; Ruled Surface(14) = {63}; Transfinite Surface {14}; //Recombine Surface {14};
Line Loop(64) = {21,-28,32}; Ruled Surface(15) = {64}; Transfinite Surface {15}; //Recombine Surface {15};
Line Loop(65) = {32,27,-24}; Ruled Surface(16) = {65}; Transfinite Surface {16}; //Recombine Surface {16};
Line Loop(66) = {27,-31,23}; Ruled Surface(17) = {66}; Transfinite Surface {17}; //Recombine Surface {17};
Line Loop(67) = {31,28,22}; Ruled Surface(18) = {67}; Transfinite Surface {18}; //Recombine Surface {18};
// outer surface of sphere 2
Surface Loop(90) = {11,12,13,14,15,16,17,18};

// lines joining inner circle (1) points to outer circle (2) points
Line(13) = {2,12}; Transfinite Line {13} = n2 Using Progression rc;
Line(14) = {3,13}; Transfinite Line {14} = n2 Using Progression rc;
Line(15) = {4,14}; Transfinite Line {15} = n2 Using Progression rc;
Line(16) = {5,15}; Transfinite Line {16} = n2 Using Progression rc;
Line(17) = {6,16}; Transfinite Line {17} = n2 Using Progression rc;
Line(18) = {7,17}; Transfinite Line {18} = n2 Using Progression rc;

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

// Plane surfaces in between the two spheres
Line Loop(91) = {22, -15, -2, 14};
Plane Surface(92) = {91}; Transfinite Surface {92}; //Recombine Surface {92};
Line Loop(93) = {21, -14, -1, 13};
Plane Surface(94) = {93}; Transfinite Surface {94}; //Recombine Surface {94};
Line Loop(95) = {13, 29, -17, -9};
Plane Surface(96) = {95}; Transfinite Surface {96}; //Recombine Surface {96};
Line Loop(97) = {15, -30, -17, 10};
Plane Surface(98) = {97}; Transfinite Surface {98}; //Recombine Surface {96};
Line Loop(99) = {31, -18, -11, 15};
Plane Surface(100) = {99}; Transfinite Surface {100}; //Recombine Surface {100};
Line Loop(101) = {13, -32, -18, 12};
Plane Surface(102) = {101}; Transfinite Surface {102}; //Recombine Surface {102};
Line Loop(103) = {25, -17, -5, 14};
Plane Surface(104) = {103}; Transfinite Surface {104}; //Recombine Surface {104};
Line Loop(105) = {14, -28, -18, 8};
Plane Surface(106) = {105}; Transfinite Surface {106}; //Recombine Surface {106};
Line Loop(107) = {7, 18, -27, -16};
Plane Surface(108) = {107}; Transfinite Surface {108}; //Recombine Surface {108};
Line Loop(109) = {26, -16, -6, 17};
Plane Surface(110) = {109}; Transfinite Surface {110}; //Recombine Surface {110};
Line Loop(111) = {23, -16, -3, 15};
Plane Surface(112) = {111}; Transfinite Surface {112}; //Recombine Surface {112};
Line Loop(113) = {16, 24, -13, -4};
Plane Surface(114) = {113}; Transfinite Surface {114}; //Recombine Surface {114};

// Volumes corresponding to the boundary layer
Surface Loop(115) = {14, 98, 4, 92, 104};
Volume(116) = {115}; Transfinite Volume {116}; Recombine Volume {116};
Surface Loop(117) = {11, 1, 96, 104, 94};
Volume(118) = {117}; Transfinite Volume {118}; Recombine Volume {118};
Surface Loop(119) = {15, 5, 94, 106, 102};
Volume(120) = {119}; Transfinite Volume {120}; Recombine Volume {120};
Surface Loop(121) = {18, 8, 100, 92, 106};
Volume(122) = {121}; Transfinite Volume {122}; Recombine Volume {122};
Surface Loop(123) = {16, 6, 102, 108, 114};
Volume(124) = {123}; Transfinite Volume {124}; Recombine Volume {124};
Surface Loop(125) = {17, 7, 108, 100, 112};
Volume(126) = {125}; Transfinite Volume {126}; Recombine Volume {126};
Surface Loop(127) = {3, 13, 112, 110, 98};
Volume(128) = {127}; Transfinite Volume {128}; Recombine Volume {128};
Surface Loop(129) = {2, 12, 110, 114, 96};
Volume(130) = {129}; Transfinite Volume {130}; Recombine Volume {130};

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
Volume(144) = {90, 143};

// Interior sphere boundary
Physical Surface(2) = {4, 1, 5, 8, 7, 3, 2, 6};
Physical Surface(3) = {11,12,13,14,15,16,17,18};
Physical Volume(1) = {144, 116,118,120,122,124,126,128,130};