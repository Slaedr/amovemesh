#include <arotation2d.hpp>
#include "alinelast_p1_sparse2.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

const int pi = 3.14159265;

int main()
{
	ifstream confile("elast2.control");
	string dum;

	double E, nu, ydisp, tol; string inmesh, outmesh, finmesh; int maxiter;

	confile >> dum; confile >> inmesh;
	confile >> dum; confile >> finmesh;
	confile >> dum; confile >> E;
	confile >> dum; confile >> nu;
	confile >> dum; confile >> ydisp;		// note that ydisp is rotation angle
	confile >> dum; confile >> tol;
	confile >> dum; confile >> maxiter;
	confile.close();

	double lambda = nu*E/((1+nu)*(1-2*nu));
	double mu = E/(2*(1+nu));

	cout << "Linear Elasticity: lambda = " << lambda << ", mu = " << mu << ".\n";

	string meshcheck = "smal.msh";

	ifstream meshfile(inmesh);
	UTriMesh m(meshfile);
	meshfile.close();
	//m.writeGmsh2(meshcheck);

	LinElastP1 prob(&m, mu, lambda);

	//set Dirichlet BCs
	Matrix<int> extra(m.gnpoin(),1);
	extra.zeros();
	//MRotation2d r(&m, ydisp*pi/180, m.gcoords(420,0), m.gcoords(420,1), n_rot);
	/*MRotation2d r(&m, ydisp*pi/180, 1.0, 0.0, 2);
	Matrix<double> bc = r.rhsvect_rotate();*/

	// set bc
	Matrix<double> bc(m.gnpoin(),2);
	bc.zeros();
	bc(64,1) = 0.25; bc(79,1) = 0.25;
	bc(4,1) = bc(18,1) = 0.5;
	bc(65,1) = bc(78,1) = 0.75;
	bc(5,1) = bc(17,1) = 0.9;
	bc(66,1) = bc(77,1) = 0.7;
	bc(6,1) = bc(16,1) = 0.5;

	prob.assembleStiffnessMatrix();
	prob.assembleLoadVector();
	prob.dirichletBC_onAllBface(bc, extra, 2, 0);

	SpMatrix K = prob.stiffnessMatrix();
	Matrix<double> f = prob.loadVector();

	/*ofstream outf("check.dat");
	K.fprint(outf);
	outf.close();*/

	Matrix<double> xold(2*m.gnpoin(),1); xold.zeros();
	Matrix<double> u(2*m.gnpoin(), 1);		// final displacements

	u = sparseCG_d(&K, f, xold, tol, maxiter);
	//u = sparse_bicgstab(&K, f, xold, tol, maxiter);
	//u = sparsegaussseidel(&K, f, xold, tol, maxiter);
	//u = cholesky(K,f);
	m.movemesh(u);
	m.writeGmsh2(finmesh);

	string jacs = "jacobian_check.dat";
	m.compute_jacobians();
	cout << "Checking jacobians\n";
	ofstream ojac(jacs);
	m.detect_negative_jacobians(ojac);
	ojac.close();

	/*Matrix<double> disp(m.gnpoin()*2,1);
	for(int i = 0; i < m.gnpoin(); i++)
		disp(i) = bc(i,0);
	for(int i = 0; i < m.gnpoin(); i++)
		disp(m.gnpoin()+i) = bc(i,1);

	m.movemesh(disp);
	m.writeGmsh2(finmesh);*/

	cout << endl;
	return 0;
}
