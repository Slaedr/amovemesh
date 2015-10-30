#include "ageometry.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	UMesh2d m;
	m.readGmsh2("tinquad.msh", 2);
	//m.readDomn("feflo.domn.cylinder.coarse");
	m.compute_topological();
	m.compute_boundary_points();
	
	vector<int> fl(8);
	fl[0] = 23; fl[1] = 16; fl[2] = 17; fl[3] = 18; fl[4] = 19; fl[5] = 20; fl[6] = 21; fl[7] = 22;

	bool closed = true;

	CSpline csp;
	csp.setup(&m, fl, closed, true, 1e-6, 1000);
	csp.sequence();
	csp.compute();
	
	ofstream ofile("points.dat");
	int inface;
	for(int iface = 0; iface < fl.size(); iface++)
	{
		for(double u = 0; u <= 1; u+= 0.025)
			ofile << csp.getspline(fl[iface],0,u) << " " << csp.getspline(fl[iface],1,u) << '\n';
	}
	ofile.close();

	cout << endl;
	return 0;
}
