#include <amesh2b.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string mname = "feflo.domn.wing.coarse";
	ifstream mfile(mname);
	UTriMesh m(mfile);
	mfile.close();

	m.writeGmsh2("wing.msh");

	cout << endl;
	return 0;
}
