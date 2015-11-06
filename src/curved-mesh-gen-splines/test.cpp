#include <amesh2dh.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	UMesh2dh m;
	m.readGmsh2("../../input/tinquadh.msh", 2);
	m.compute_topological();
	int fac = 27;
	for(int i = 0; i < 4; i++)
		cout << m.gintfac(fac,i) << " ";
	cout << endl;

	int elem = 56-28-1;
	cout << "nfael " << m.gnfael(elem) << endl;
	for(int i = 0; i < m.gnfael(elem); i++)
		cout << m.gesuel(elem,i) << " ";
	cout << endl;

	cout << endl;
	return 0;
}
