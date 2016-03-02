srad = 1.0;
ffs = 10.0;
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

Extrude {
		Line {1,2,3,4}; Layers{ {10,10,10}, {1e-2,5e-2,2.5e-1} };
		Recombine;
		QuadTriAddVerts RecombLaterals;
	}