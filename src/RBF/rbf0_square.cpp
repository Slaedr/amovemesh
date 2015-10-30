
#include <amesh2.hpp>
#include "arbf0.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

const double pi = 3.14159265;

int main()
{
	string confile = "rbf_square.control";
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
	for(int i = 0; i < m.gnface(); i++)
	{
		for(int j = 0; j < m.gndim(); j++)
			flags(m.gbface(i,j)) = 1;
	}

	Matrix<double> bc(m.gnpoin(),2);
	bc.zeros();
	bc(64,1) = 0.25; bc(79,1) = 0.25;
	bc(4,1) = bc(18,1) = 0.5;
	bc(65,1) = bc(78,1) = 0.75;
	bc(5,1) = bc(17,1) = 0.9;
	bc(66,1) = bc(77,1) = 0.7;
	bc(6,1) = bc(16,1) = 0.5;

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
