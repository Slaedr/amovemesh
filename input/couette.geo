/*********************************************************************
 * GMSH geometry file for couette flow
 *
 * cjbuaa Dec, 29 2014
 *********************************************************************/
//define geometry
 
 
lc = 0.2;
 
Point(1) = { 0.0, 0.0, 0.0, lc};
Point(2) = { 1.0, 0.0, 0.0, lc};
Point(3) = { 1.0,  1.0, 0.0, lc};
Point(4) = {0.0,  1.0, 0.0, lc};
 
Line(5) = {1,2}; Transfinite Line {5} = 9;
Line(6) = {2,3}; Transfinite Line {6} = 9;
Line(7) = {3,4}; Transfinite Line {7} = 9;
Line(8) = {4,1}; Transfinite Line {8} = 9;
 
Line Loop(1) = {5,6,7,8};
Ruled Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface(1);
 
 //Mesh.Color.Triangles{160} 
Physical Line(01) = {5};
Physical Line(04) = {6};
Physical Line(02) = {7};
Physical Line(03) = {8};
Physical Surface(100) = {1}; 
