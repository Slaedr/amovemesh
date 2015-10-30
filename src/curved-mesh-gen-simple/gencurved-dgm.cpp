#include "acurvedmeshgen-dg.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "gencurved-dgm.control";
	ifstream conf(confile);
	string inmesh;
	string intermesh;
	string outmesh;
	string dum;

	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> intermesh;
	conf >> dum; conf >> outmesh;

	UMesh2d m;
	m.readGmsh2(inmesh, 2);
	m.compute_topological();
	UMesh2d mq;
	mq = m.convertLinearToQuadratic();
	mq.printmeshstats();
	/*ofstream outf("bface.dat");
	for(int i = 0; i < mq.gnface(); i++)
	{
		for(int j = 0; j < mq.gnnofa()+mq.gnbtag(); j++)
			outf << " " << mq.gbface(i,j);
		outf << '\n';
	}
	outf.close();*/
	mq.writeGmsh2(intermesh);

	Matrix<int> bounflags(2,1); bounflags(0) = 9; bounflags(1) = 1;
	//Matrix<int> bounflags(1,1); bounflags(0) = 3;
	cout << endl;
	//CurvedMeshGeneration cmg(&mq, bounflags);
	CurvedMeshGeneration cmg(&mq, bounflags);
	cmg.compute_boundary_displacement();
	cmg.generate();
	mq.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
