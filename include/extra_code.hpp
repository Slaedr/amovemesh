/** Unused, mostly useless, code. */

/// Computes circumcentre and circumradius of a [tetrahedron](@ref Tet).
/** NOTE: The tetrahedron's jacobian (2*volume) should be stored in [elem.D](&ref D) beforehand.
* The center and radius of the circumsphere are calculated as follows. The circumcenter is
* \f[ \mathbf{O} = \frac{|\mathbf{a}^2(\mathbf{b}\times \mathbf{c}) + \mathbf{b}^2(\mathbf{c}\times \mathbf{a}) + \mathbf{c}^2(\mathbf{a}\times \mathbf{b})|}{12 V} \f]
* and the radius \f$ R = |\mathbf{O}|\f$. 
*/
void Delaunay3d::compute_circumsphere2(Tet& elem)
{
	vector<double> a(ndim), b(ndim), c(ndim), n1(ndim), n2(ndim), n3(ndim), fin(ndim);
	for(int idim = 0; idim < ndim; idim++)
	{
		a[idim] = nodes[elem.p[0]][idim]-nodes[elem.p[3]][idim];
		b[idim] = nodes[elem.p[1]][idim]-nodes[elem.p[3]][idim];
		c[idim] = nodes[elem.p[2]][idim]-nodes[elem.p[3]][idim];
	}

	cross_product3(n1,b,c);
	cross_product3(n2,c,a);
	cross_product3(n3,a,b);
	for(int idim = 0; idim < ndim; idim++)
	{
		fin[idim] = dot(a,a)*n1[idim] + dot(b,b)*n2[idim] + dot(c,c)*n3[idim];
		fin[idim] /= 2.0*elem.D;
		elem.centre[idim] = fin[idim];
	}
	elem.radius = l2norm(fin);
	cout << "Circumsphere data: centre " << elem.centre[0] << "," << elem.centre[1] << "," << elem.centre[2] << ", radius " << elem.radius << ", elem jacobian " << elem.D << endl;
}

/// Computes circumcentre and circumradius of a [tetrahedron](@ref Tet).
/** NOTE: The tetrahedron's jacobian (6*volume) should be stored in [elem.D](@ref elem::D) beforehand.
* Reference: [Weisstein, Eric W. "Circumsphere." From MathWorld--A Wolfram Web Resource.](@ref http://mathworld.wolfram.com/Circumsphere.html)
* \note Probably wrong!
*/
void Delaunay3d::compute_circumsphere3(Tet& elem)
{
	Matrix<double> dx(nnode,nnode), dy(nnode,nnode), dz(nnode,nnode), cc(nnode,nnode);
	double a = elem.D, c, d_x, d_y, d_z;
	#if DEBUGW==1
	if(fabs(a) < ZERO_TOL) cout << "Delaunay3d: compute_circumcircle(): ! Jacobian of the element is zero!!" << endl;
	#endif
	
	dx(0,0) = nodes[elem.p[0]][0]*nodes[elem.p[0]][0] + nodes[elem.p[0]][1]*nodes[elem.p[0]][1] + nodes[elem.p[0]][2]*nodes[elem.p[0]][2];
	dx(0,1) = nodes[elem.p[0]][1]; dx(0,2) = nodes[elem.p[0]][2]; dx(0,3) = 1;
	dx(1,0) = nodes[elem.p[1]][0]*nodes[elem.p[1]][0] + nodes[elem.p[1]][1]*nodes[elem.p[1]][1] + nodes[elem.p[1]][2]*nodes[elem.p[1]][2];
	dx(1,1) = nodes[elem.p[1]][1]; dx(1,2) = nodes[elem.p[1]][2]; dx(1,3) = 1;
	dx(2,0) = nodes[elem.p[2]][0]*nodes[elem.p[2]][0] + nodes[elem.p[2]][1]*nodes[elem.p[2]][1] + nodes[elem.p[2]][2]*nodes[elem.p[2]][2];
	dx(2,1) = nodes[elem.p[2]][1]; dx(2,2) = nodes[elem.p[2]][2]; dx(2,3) = 1;
	dx(3,0) = nodes[elem.p[3]][0]*nodes[elem.p[3]][0] + nodes[elem.p[3]][1]*nodes[elem.p[3]][1] + nodes[elem.p[3]][2]*nodes[elem.p[3]][2];
	dx(3,1) = nodes[elem.p[3]][1]; dx(3,2) = nodes[elem.p[3]][2]; dx(3,3) = 1;
	
	dy(0,0) = nodes[elem.p[0]][0]*nodes[elem.p[0]][0] + nodes[elem.p[0]][1]*nodes[elem.p[0]][1] + nodes[elem.p[0]][2]*nodes[elem.p[0]][2];
	dy(0,1) = nodes[elem.p[0]][0]; dy(0,2) = nodes[elem.p[0]][2]; dy(0,3) = 1;
	dy(1,0) = nodes[elem.p[1]][0]*nodes[elem.p[1]][0] + nodes[elem.p[1]][1]*nodes[elem.p[1]][1] + nodes[elem.p[1]][2]*nodes[elem.p[1]][2];
	dy(1,1) = nodes[elem.p[1]][0]; dy(1,2) = nodes[elem.p[1]][2]; dy(1,3) = 1;
	dy(2,0) = nodes[elem.p[2]][0]*nodes[elem.p[2]][0] + nodes[elem.p[2]][1]*nodes[elem.p[2]][1] + nodes[elem.p[2]][2]*nodes[elem.p[2]][2];
	dy(2,1) = nodes[elem.p[2]][0]; dy(2,2) = nodes[elem.p[2]][2]; dy(2,3) = 1;
	dy(3,0) = nodes[elem.p[3]][0]*nodes[elem.p[3]][0] + nodes[elem.p[3]][1]*nodes[elem.p[3]][1] + nodes[elem.p[3]][2]*nodes[elem.p[3]][2];
	dy(3,1) = nodes[elem.p[3]][0]; dy(3,2) = nodes[elem.p[3]][2]; dy(3,3) = 1;
	
	dz(0,0) = nodes[elem.p[0]][0]*nodes[elem.p[0]][0] + nodes[elem.p[0]][1]*nodes[elem.p[0]][1] + nodes[elem.p[0]][2]*nodes[elem.p[0]][2];
	dz(0,1) = nodes[elem.p[0]][0]; dz(0,2) = nodes[elem.p[0]][1]; dz(0,3) = 1;
	dz(1,0) = nodes[elem.p[1]][0]*nodes[elem.p[1]][0] + nodes[elem.p[1]][1]*nodes[elem.p[1]][1] + nodes[elem.p[1]][2]*nodes[elem.p[1]][2];
	dz(1,1) = nodes[elem.p[1]][0]; dz(1,2) = nodes[elem.p[1]][1]; dz(1,3) = 1;
	dz(2,0) = nodes[elem.p[2]][0]*nodes[elem.p[2]][0] + nodes[elem.p[2]][1]*nodes[elem.p[2]][1] + nodes[elem.p[2]][2]*nodes[elem.p[2]][2];
	dz(2,1) = nodes[elem.p[2]][0]; dz(2,2) = nodes[elem.p[2]][1]; dz(2,3) = 1;
	dz(3,0) = nodes[elem.p[3]][0]*nodes[elem.p[3]][0] + nodes[elem.p[3]][1]*nodes[elem.p[3]][1] + nodes[elem.p[3]][2]*nodes[elem.p[3]][2];
	dz(3,1) = nodes[elem.p[3]][0]; dz(3,2) = nodes[elem.p[3]][1]; dz(3,3) = 1;
	
	cc(0,0) = nodes[elem.p[0]][0]*nodes[elem.p[0]][0] + nodes[elem.p[0]][1]*nodes[elem.p[0]][1] + nodes[elem.p[0]][2]*nodes[elem.p[0]][2];
	cc(0,1) = nodes[elem.p[0]][0]; cc(0,2) = nodes[elem.p[0]][1]; cc(0,3) = nodes[elem.p[0]][2];
	cc(1,0) = nodes[elem.p[1]][0]*nodes[elem.p[1]][0] + nodes[elem.p[1]][1]*nodes[elem.p[1]][1] + nodes[elem.p[1]][2]*nodes[elem.p[1]][2];
	cc(1,1) = nodes[elem.p[1]][0]; cc(1,2) = nodes[elem.p[1]][1]; cc(1,3) = nodes[elem.p[1]][2];
	cc(2,0) = nodes[elem.p[2]][0]*nodes[elem.p[2]][0] + nodes[elem.p[2]][1]*nodes[elem.p[2]][1] + nodes[elem.p[2]][2]*nodes[elem.p[2]][2];
	cc(2,1) = nodes[elem.p[2]][0]; cc(2,2) = nodes[elem.p[2]][1]; cc(2,3) = nodes[elem.p[2]][2];
	cc(3,0) = nodes[elem.p[3]][0]*nodes[elem.p[3]][0] + nodes[elem.p[3]][1]*nodes[elem.p[3]][1] + nodes[elem.p[3]][2]*nodes[elem.p[3]][2];
	cc(3,1) = nodes[elem.p[3]][0]; cc(3,2) = nodes[elem.p[3]][1]; cc(3,3) = nodes[elem.p[3]][2];

	d_x = determinant(dx);
	d_y = determinant(dy);
	d_z = determinant(dz);
	c = determinant(cc);

	elem.centre[0] = d_x/a*0.5;
	elem.centre[1] = d_y/a*0.5;
	elem.centre[2] = d_z/a*0.5;
	elem.radius = sqrt(d_x*d_x + d_y*d_y + d_z*d_z - 4.0*a*c)/2*fabs(a);

	cout << "Circumcircle data: centre " << elem.centre[0] << "," << elem.centre[1] << "," << elem.centre[2] << ", radius " << elem.radius << ", elem jacobian " << a << endl;
}
