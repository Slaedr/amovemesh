#include "acurvedmeshgen-rbf-3d.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "gencurved-rbf-3d.control";
	ifstream conf(confile);
	string inmesh, intermesh, outmesh, dum;
	double sup_rad; int num_steps, numbflags;

	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> outmesh;
	conf >> dum; conf >> sup_rad;
	conf >> dum; conf >> num_steps;
	conf >> dum; conf >> numbflags;
	conf >> dum;
	Matrix<int> bounflags(numbflags,1);
	for(int i = 0; i < numbflags; i++)
		conf >> bounflags(i);
	conf.close();
	//sup_rad = stod(str_sup_rad);
	//cout << sup_rad << endl;

	UMesh m;
	m.readGmsh2(inmesh, 3);
	CurvedMeshGeneration cmg(&m, bounflags, num_steps, 2, sup_rad);
	cmg.compute_boundary_displacement();
	cmg.generate();
	m.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
