#include <amesh2dh.hpp>

using namespace amat;
using namespace acfd;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 3) {
		cout << "Please give two string arguments - input file and output file." << endl;
		return -1;
	}
	string inmesh(argv[1]);
	string outmesh(argv[2]);
	UMesh2dh m;
	m.readGmsh2(inmesh,2);

	UMesh2dh tm = m.convertQuadToTri();
	tm.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
