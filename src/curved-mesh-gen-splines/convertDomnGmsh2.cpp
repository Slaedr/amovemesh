#include <amesh2d.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string dum, mfile, omfile;
	
	mfile = "../../input/smallmesh.domn";
	omfile = "../../input/smallmesh.msh";

	UMesh2d m;
	m.readDomn(mfile);
	m.writeGmsh2(omfile);

	cout << endl;
	return 0;
}
