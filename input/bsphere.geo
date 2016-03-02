/* GMSH geometry file for boundary layer mesh of sphere
 *
 * Aditya Kashi
 */
 
lc = 0.2;

//mesh size at far field
lf = 1.0;

// sphere radius
srad = 1.0;
// farfield coordinate
ffsize = 10.0;
 
// center
Point(1) = { 0.0, 0.0, 0.0, lc};

// points for arcs of circle in xy plane
Point(2) = { srad, 0.0, 0.0, lc};
Point(3) = { 0.0,  srad, 0.0, lc};
Point(4) = {-srad,  0.0, 0.0, lc};
Point(5) = {0.0, -srad, 0.0, lc};
Point(6) = {0.0, 0.0, srad, lc};
Point(7) = {0.0, 0.0, -srad, lc};

/*Point(10) = {ffsize, ffsize, -ffsize, lf};
Point(11) = {-ffsize, ffsize, -ffsize, lf};
Point(12) = {-ffsize, -ffsize, -ffsize, lf};
Point(13) = {ffsize, -ffsize, -ffsize, lf};
Point(14) = {ffsize, ffsize, ffsize, lf};
Point(15) = {-ffsize, ffsize, ffsize, lf};
Point(16) = {-ffsize, -ffsize, ffsize, lf};
Point(17) = {ffsize, -ffsize, ffsize, lf};*/
 
Circle(30) = {2,1,3}; Transfinite Line {30} = 9;
Circle(31) = {3,1,4}; Transfinite Line {31} = 9;
Circle(32) = {4,1,5}; Transfinite Line {32} = 9;
Circle(33) = {5,1,2}; Transfinite Line {33} = 9;
 
Line Loop(50) = {30,31,32,33};
Ruled Surface(50) = {50};
Transfinite Surface{50};
//Recombine Surface(50);

// rotate-extrude the disc to get a sphere
Extrude{ {1,0,0}, {1,0,0}, 2*Pi } {Surface{50};}

/*Surface Loop(50) = {extrusion[0]};

Volume(100) = {50};
Transfinite Volume{100};
//Recombine Volume(100);

Physical Volume(200) = {100};*/

 //Mesh.Color.Triangles{160} 
/*Physical Line(01) = {10};
Physical Line(04) = {11};
Physical Line(02) = {12};
Physical Line(03) = {13};*/

//Physical Surface(100) = {50}; 
