#include <amesh2.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	ifstream infile("rdgflo.domn.32x32");
	UTriMesh m(infile);
	infile.close();
	m.writeGmsh2("rdgflo32x32.msh");

	return 0;
}
