#include "amesh3d.hpp"

using namespace amat;
using namespace amc;
using namespace std;

int main()
{
	UMesh m;
	m.readGmsh2("../input/smalltet.msh",3);
	m.compute_topological();
	int offs = 57;
	
	int gp = 7;
	int poin = m.gbpointsinv(gp);
	
	cout << "Faces surrounding b point " << poin+1 << " (point " << gp+1 << ") are" << endl;
	for(int isurr = m.gbfsubp_p(poin); isurr < m.gbfsubp_p(poin+1); isurr++)
		cout << " " << m.gbfsubp(isurr) + 57;
	cout << endl;
	
	cout << endl;
	return 0;
}
