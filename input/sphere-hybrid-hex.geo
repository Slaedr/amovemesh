srad = 1.0;
ffs = 10.0;
num_points = 9;

hs = 0.2;
hf = 1.0;
Point(1) = {0, 0, 0, hs};
Point(2) = {srad, 0, 0, hs};
Point(3) = {0, srad, 0, hs};
Point(4) = {-srad, 0, 0, hs};
Point(5) = {0, -srad, 0, hs};
Point(6) = {0, 0, srad, hs};
Point(7) = {0, 0, -srad, hs};

Point(8) = {ffs, ffs, -ffs, hf};
Point(9) = {-ffs, ffs, -ffs, hf};
Point(10) = {-ffs, -ffs, -ffs, hf};
Point(11) = {ffs, -ffs, -ffs, hf};
Point(12) = {ffs, ffs, ffs, hf};
Point(13) = {-ffs, ffs, ffs, hf};
Point(14) = {-ffs, -ffs, ffs, hf};
Point(15) = {ffs, -ffs, ffs, hf};

Circle(1) = {6, 1, 2};
Circle(2) = {2, 1, 7};
Circle(3) = {7, 1, 4};
Circle(4) = {4, 1, 6}; 
Circle(5) = {3, 1, 6};
Circle(6) = {6, 1, 5};
Circle(7) = {5, 1, 7}; 
Circle(8) = {7, 1, 3}; 
Circle(9) = {4, 1, 3}; 
Circle(10) = {3, 1, 2};
Circle(11) = {2, 1, 5};
Circle(12) = {5, 1, 4};

Line Loop(15) = {9, 5, -4};
Ruled Surface(15) = {15}; Recombine Surface(15);
Line Loop(17) = {5, 1, -10};
Ruled Surface(17) = {17}; Recombine Surface(17);
Line Loop(19) = {10, 2, 8};
Ruled Surface(19) = {19}; Recombine Surface(19);
Line Loop(21) = {8, -9, -3};
Ruled Surface(21) = {21}; Recombine Surface(21);
Line Loop(23) = {6, 12, 4};
Ruled Surface(23) = {23}; Recombine Surface(23);
Line Loop(25) = {6, -11, -1};
Ruled Surface(25) = {25}; Recombine Surface(25);
Line Loop(27) = {11, 7, -2};
Ruled Surface(27) = {27}; Recombine Surface(27);
Line Loop(29) = {7, 3, -12};
Ruled Surface(29) = {29}; Recombine Surface(29);

Surface Loop(31) = {19, 17, 15, 21, 29, 27, 25, 23};
//Volume(31) = {31}; Transfinite Volume{31};

Line(32) = {13, 12};
Line(33) = {12, 15};
Line(34) = {15, 14};
Line(35) = {14, 13};
Line(36) = {9, 8};
Line(37) = {8, 11};
Line(38) = {11, 10};
Line(39) = {10, 9};
Line(40) = {9, 13};
Line(41) = {8, 12};
Line(42) = {10, 14};
Line(43) = {11, 15};

Line Loop(44) = {40, -35, -42, 39};
Plane Surface(45) = {44};
Line Loop(46) = {35, 32, 33, 34};
Plane Surface(47) = {46};
Line Loop(48) = {37, 43, -33, -41};
Plane Surface(49) = {48};
Line Loop(50) = {42, -34, -43, 38};
Plane Surface(51) = {50};
Line Loop(52) = {40, 32, -41, -36};
Plane Surface(53) = {52};
Line Loop(54) = {36, 37, 38, 39};
Plane Surface(55) = {54};

Surface Loop(56) = {45, 53, 47, 49, 55, 51};

Volume(57) = {31, 56};

Extrude {
	Surface {19, 17, 15, 21, 29, 27, 25, 23}; Layers{ {10,10}, {1e-2,5e-2} };
	Recombine;
}

Mesh.Algorithm = 6; // Frontal

