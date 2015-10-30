#include "alinelast_p2.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	ifstream confile("elast.control");
	string dum;

	double E, nu, ydisp, tol; string inmesh, outmesh, finmesh; int maxiter;

	confile >> dum; confile >> inmesh;
	confile >> dum; confile >> finmesh;
	confile >> dum; confile >> E;
	confile >> dum; confile >> nu;
	confile >> dum; confile >> ydisp;
	confile >> dum; confile >> tol;
	confile >> dum; confile >> maxiter;
	confile.close();

	double lambda = nu*E/((1+nu)*(1-2*nu));
	double mu = E/(2*(1+nu));

	cout << "Linear Elasticity: lambda = " << lambda << ", mu = " << mu << ".\n";

	string meshcheck = "outmesh.msh";

	ifstream meshfile(inmesh);
	UTriMeshCurved m;
	m.readDomn(meshfile);
	meshfile.close();
	//m.writeQuadraticGmsh2(meshcheck);

	LinElastP2 prob(&m, mu, lambda);

	prob.assembleStiffnessMatrix();
	prob.assembleLoadVector();

	//set Dirichlet BCs
	Matrix<double> bdata(2*m.gnface(),1); Matrix<int> extra(m.gnpoin(),1);
	bdata.zeros(); extra.zeros();
	extra(12) = extra(22) = 1;
	bdata(0) = bdata(1) = 0.0; bdata(m.gnface()+0) = bdata(m.gnface()+1) = ydisp;
	prob.dirichletBC_onAllBface(bdata, extra);

	Matrix<double> K = prob.stiffnessMatrix();
	Matrix<double> f = prob.loadVector();

	ofstream outf("check.dat");
	K.fprint(outf);
	outf << '\n';
	f.fprint(outf);
	outf.close();

	Matrix<double> xold(2*m.gnpoin(),1); xold.zeros();
	/*xold(m.gnpoin()+1) = ydisp; xold(m.gnpoin()+3) = ydisp; xold(m.gnpoin()+6) = 0.2; xold(m.gnpoin()+8) = 0.2;
	xold(m.gnpoin()+11) = 0.15; xold(m.gnpoin()+13) = 0.15;*/
	Matrix<double> u(2*m.gnpoin(), 1);		// final displacements

	u = gaussseidel(K, f, xold, tol, maxiter, 'n');
	//u = pointjacobi(K, f, xold, tol, maxiter, 'n');
	//u = cholesky(K,f);
	m.movemesh(u);
	m.writeQuadraticGmsh2(finmesh);

	cout << endl;
	return 0;
}
