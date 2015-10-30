#include <aconstants.h>
#include "adgrbf.hpp"
#include <arotation2d.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "dg_rotate.control";
	string inp, outp, outdg, jacs, dum, anglestr;
	double suprad;
	int nmarks;		// number of boundary markers to read
	vector<double> centre(2);
	Matrix<int> n_rot;
	ifstream conf(confile);
	conf >> dum; conf >> inp;
	conf >> dum; conf >> outp;
	conf >> dum; conf >> outdg;
	conf >> dum; conf >> jacs;
	conf >> dum; conf >> anglestr;
	conf >> dum; conf >> centre[0] >> centre[1];
	conf >> dum; conf >> nmarks;

	conf >> dum;
	n_rot.setup(nmarks,1);
	for(int i = 0; i < nmarks; i++)
		conf >> n_rot(i);
	
	conf >> dum; conf >> suprad;
	conf.close();
	double angle = stod(anglestr)*PI/180.0;			// convert to radians
	cout << "support radius is " << suprad << endl;
	cout << "Centre is " << centre[0] << " " << centre[1] << endl;

	// read mesh
	UMesh2d m;
	m.readGmsh2(inp,2);

	Matrix<int> flags(m.gnpoin(),1);
	flags.zeros();
	for(int i = 0; i < m.gnface(); i++)
		for(int j = 0; j < m.gnnofa(); j++)
			flags(m.gbface(i,j)) = 1;

	Matrix<double>* coord = m.getcoords();

	//calculate boundary displacement
	MRotation2d rot(&m, angle, centre[0], centre[1], n_rot);
	Matrix<double> bc = rot.rhsvect_angles();

	// carry out DG mapping procedure
	DGRBFrotate d(2, coord, flags, &bc, centre, 2, suprad);
	d.generateDG();
	d.movemesh();
	Matrix<double> newcoords = d.getcoords();
	m.setcoords(&newcoords);
	m.writeGmsh2(outp);

	d.movedg();
	d.dg.writeGmsh2(outdg);

	// check validity
	d.dg.compute_jacobians();
	cout << "Checking DG jacobians\n";
	bool check = d.dg.detect_negative_jacobians();

	m.compute_jacobians();
	cout << "Checking final mesh jacobians\n";
	ofstream ojac(jacs);
	m.detect_negative_jacobians(ojac);
	ojac.close();

	cout << "Done.\n" << endl;
	return 0;
}
