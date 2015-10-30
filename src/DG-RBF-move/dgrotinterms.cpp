#include "adgrbf_rot_inter_ms.hpp"

using namespace amat;
using namespace acfd;
using namespace std;

const double pi = 3.14159265;

int main()
{
	string conf = "dgrotinterms.control";
	string dummy, inmesh, outmesh, outdg, jacf;
	double angledeg, suprad;
	int rbf_choice;
	ifstream confile(conf);
	confile >> dummy; confile >> inmesh;
	confile >> dummy; confile >> outmesh;
	confile >> dummy; confile >> outdg;
	confile >> dummy; confile >> jacf;
	confile >> dummy; confile >> angledeg;
	confile >> dummy; confile >> rbf_choice;
	confile >> dummy; confile >> suprad;
	confile.close();

	ifstream mmesh(inmesh);
	UTriMesh m(mmesh);
	mmesh.close();

	Matrix<int> nrot(1,1); nrot(0) = 3;
	int rotpoint = 420;

	DGRBFmove dm(&m, rbf_choice, suprad, nrot, angledeg*pi/180, m.gcoords(rotpoint,0), m.gcoords(rotpoint,1));
	dm.movemesh();
	Matrix<double> newcoords = dm.getcoords();
	m.setcoords(&newcoords);
	m.writeGmsh2(outmesh);
	dm.dg.writeGmsh2(outdg);

	m.compute_jacobians();
	cout << "Checking jacobians\n";
	ofstream ojac(jacf);
	m.detect_negative_jacobians(ojac);
	ojac.close();

	cout << endl;
	return 0;
}
