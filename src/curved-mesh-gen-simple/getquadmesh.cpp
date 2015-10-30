#include <amesh3d.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "getquadmesh.control";
	ifstream conf(confile);
	string inmesh, intermesh, outmesh, dum;

	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> outmesh;
	conf.close();
	cout << inmesh << endl;

	UMesh m;
	m.readGmsh2(inmesh, 3);
	m.compute_topological();
	UMesh mq;
	mq = m.convertLinearToQuadratic();
	mq.printmeshstats();

	mq.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
