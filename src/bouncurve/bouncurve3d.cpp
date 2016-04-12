#include "abouncurve3d.hpp"

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Please give a control file name!\n";
		return -1;
	}
	string confile(argv[1]);
	ifstream conf(confile);
	string inmesh, intermesh, outmesh, rbf_solver, dum;
	double sup_rad, rbf_tol; int num_steps, numbflags, rbf_maxiter;

	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> outmesh;
	/*conf >> dum; conf >> sup_rad;
	conf >> dum; conf >> num_steps;
	conf >> dum; conf >> rbf_solver;
	conf >> dum; conf >> rbf_tol;
	conf >> dum; conf >> rbf_maxiter;*/
	conf >> dum; conf >> numbflags;
	conf >> dum;
	Matrix<int> bounflags(numbflags,1);
	for(int i = 0; i < numbflags; i++)
		conf >> bounflags(i);
	conf.close();
	//sup_rad = stod(str_sup_rad);
	//cout << sup_rad << endl;
	
	cout << "Number of fixed boundary flags = " << numbflags << endl;
	cout << "They are ";
	for(int i = 0; i < numbflags; i++)
		cout << bounflags(i);
	cout << endl;

	UMesh m;
	m.readGmsh2(inmesh, 3);
	BounCurve cmg(&m, bounflags);
	cmg.compute_boundary_displacement();
	cmg.generate();
	m.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
