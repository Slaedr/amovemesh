
#include <arotation2db.hpp>
#include "arbf0.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

const double pi = 3.14159265;

int main()
{
	string confile = "rbf_wing.control";
	
	string inp, outp, jacs, dum, anglestr;
	int maxiter, rbf_choice, num_steps;
	double support_radius, tol;

	ifstream conf(confile);
	conf >> dum; conf >> inp;
	conf >> dum; conf >> outp;
	conf >> dum; conf >> jacs;
	conf >> dum; conf >> anglestr;
	conf >> dum; conf >> tol;
	conf >> dum; conf >> maxiter;
	conf >> dum; conf >> rbf_choice;
	conf >> dum; conf >> support_radius;
	conf >> dum; conf >> num_steps;
	conf.close();
	double angle = stod(anglestr);

	cout << "RBF choice is " << rbf_choice << endl;

	// read mesh
	ifstream origmesh(inp);
	UTriMesh m(origmesh);
	origmesh.close();

	Matrix<int> flags(m.gnpoin(),1);
	flags.zeros();
	for(int i = 0; i < m.gnbpoin(); i++)
	{
		flags(i) = 1;
	}

	//calculate boundary motion
	Matrix<int> n_rot(1,1);
	n_rot(0) = 3;		// flag of body to be rotated
	MRotation2d r(&m, angle*pi/180, m.gcoords(420,0), m.gcoords(420,1), n_rot);
	Matrix<double> bc = r.rhsvect_rotate();

	// carry out RBF mapping procedure
	Matrix<double>* coords = m.getcoords();
	RBFmove rm(coords, &bc, flags, rbf_choice, support_radius, num_steps, tol, maxiter);
	rm.move();
	Matrix<double> newcoords = rm.getcoords();
	m.setcoords(&newcoords);
	m.writeGmsh2(outp);

	m.compute_jacobians();
	cout << "Checking jacobians\n";
	ofstream ojac(jacs);
	m.detect_negative_jacobians(ojac);
	ojac.close();

	cout << "\nDone.\n" << endl;
	return 0;
}
