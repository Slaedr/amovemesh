#include <amesh2b.hpp>
#include "adg_multistep_mm.hpp"
#include <arotation2db.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

const double pi = 3.14159265;

int main()
{
	string confile = "dg_ms_rotate.control";
	string inp, outp, outdg, jacs, dum, anglestr;
	ifstream conf(confile);
	conf >> dum; conf >> inp;
	conf >> dum; conf >> outp;
	conf >> dum; conf >> outdg;
	conf >> dum; conf >> jacs;
	conf >> dum; conf >> anglestr;
	conf.close();
	double angle = stod(anglestr);

	// read mesh
	ifstream origmesh(inp);
	UTriMesh m(origmesh);
	origmesh.close();

	Matrix<int> flags(m.gnpoin(),1);
	flags.zeros();
	for(int i = 0; i < m.gnbpoin(); i++)
		flags(i) = 1;

	Matrix<double>* coord = m.getcoords();
	Matrix<int> n_rot(1,1); n_rot(0) = 3;		// flag of body to be rotated
	MRotation2d r(&m, angle*pi/180, m.gcoords(420,0), m.gcoords(420,1), n_rot);
	Matrix<double> bc = r.rhsvect_rotate();

	// carry out DG mapping procedure
	DGmove d(2, coord, flags, &bc);
	d.movemesh();
	Matrix<double> newcoords = d.getcoords();
	m.setcoords(&newcoords);
	m.writeGmsh2(outp);
	d.dg.writeGmsh2(outdg);

	// check validity
	/*d.dg.compute_jacobians();
	cout << "Checking DG jacobians\n";
	bool check = d.dg.detect_negative_jacobians();*/

	m.compute_jacobians();
	cout << "Checking jacobians\n";
	ofstream ojac(jacs);
	m.detect_negative_jacobians(ojac);
	ojac.close();

	cout << "Done.\n" << endl;
	return 0;
}
