#include <amesh2d.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "quadratize2d.control";
	ifstream conf(confile);
	string inmesh, outmesh, dum;
	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> outmesh;
	conf.close();

	cout << "Reading file " << inmesh << endl;
	UMesh2d m, mq;
	m.readGmsh2(inmesh,2);
	m.compute_topological();
	mq = m.convertLinearToQuadratic();
	mq.writeGmsh2(outmesh);

	cout << "Quadratic mesh written to file " << outmesh << endl;

	cout << endl;
	return 0;
}
