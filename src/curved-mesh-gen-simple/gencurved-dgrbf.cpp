#include "acurvedmeshgen-dgrbf.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "gencurved-dgrbf.control";
	ifstream conf(confile);
	string inmesh;
	string intermesh;
	string outmesh;
	string dum;
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
	//conf >> dum; conf >> num_steps;
	conf.close();

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
	cout << "--Now curving the mesh--\n";
	//CurvedMeshGeneration cmg(&mq, bounflags);
	CurvedMeshGeneration cmg(&mq, bounflags, 2, sup_rad);
	cmg.compute_boundary_displacement();
	cmg.generate();
	mq.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
