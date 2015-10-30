#include "acurvedmeshgen-rbf.hpp"
#include <aoutput.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "gencurved-rbf.control";
	ifstream conf(confile);
	string inmesh, intermesh, outmesh, dum;
	double sup_rad; int num_steps, numbflags;

	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> intermesh;
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

	UMesh2d m;
	m.readGmsh2(inmesh, 2);
	m.compute_topological();
	UMesh2d mq;
	mq = m.convertLinearToQuadratic();
	mq.printmeshstats();

	mq.writeGmsh2(intermesh);
	//writeQuadraticMeshToVtu(intermesh.replace(intermesh.end()-3, intermesh.end(), "vtu"), mq);
	
	cout << "--Now curving the mesh--\n";

	CurvedMeshGeneration cmg(&mq, bounflags, num_steps, 2, sup_rad);
	cmg.compute_boundary_displacement();
	cmg.generate();
	mq.writeGmsh2(outmesh);
	
	//writeQuadraticMeshToVtu(outmesh.replace(outmesh.end()-3, outmesh.end(), "vtu"), mq);

	cout << endl;
	return 0;
}
