#include <ageometry.hpp>
#include "acurvedmeshgen-rbf-hs.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "gencurved.control";
	ifstream conf(confile);
	string inmesh, intermesh, outmesh, splinefile, dum;
	double sup_rad; int num_steps, numbflags;

	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> intermesh;
	conf >> dum; conf >> splinefile;
	conf >> dum; conf >> outmesh;
	conf >> dum; conf >> sup_rad;
	conf >> dum; conf >> num_steps;
	conf >> dum; conf >> numbflags;
	conf >> dum;
	Matrix<int> bounflags(numbflags,1);
	for(int i = 0; i < numbflags; i++)
		conf >> bounflags(i);
	conf.close();
	//cout << sup_rad << endl;
	//cout << bounflags(0) << endl;

	UMesh2d m;
	m.readGmsh2(inmesh, 2);
	m.compute_topological();
	UMesh2d mq;
	mq = m.convertLinearToQuadratic();
	mq.printmeshstats();

	mq.writeGmsh2(intermesh);

	cout << "--Now curving the mesh--\n";

	CurvedMeshGeneration cmg(&mq, bounflags, splinefile, num_steps, 2, sup_rad);
	//cmg.generate();
	//mq.writeGmsh2(outmesh);
	
	cout << endl;
	return 0;
}
