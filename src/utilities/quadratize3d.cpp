#include <amesh3d.hpp>

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "No control file name provided!\n";
		return -1;
	}
	string confile(argv[1]);
	ifstream conf(confile);
	string inmesh, outmesh, dum;
	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> outmesh;
	conf.close();

	cout << "Reading file " << inmesh << endl;
	UMesh m, mq;
	m.readGmsh2(inmesh,3);
	m.compute_topological();
	mq = m.convertLinearToQuadratic();
	mq.writeGmsh2(outmesh);

	cout << "Quadratic mesh written to file " << outmesh << endl;

	cout << endl;
	return 0;
}
