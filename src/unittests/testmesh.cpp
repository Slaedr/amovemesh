#include "amesh2dh.hpp"

using namespace amat;
using namespace acfd;
using namespace std;

int main()
{
	string mfilename = "testhybrid.msh";
	UMesh2dh m;
	m.readGmsh2(mfilename,2);
	cout << endl;
	m.writeGmsh2("testhybrid2.msh");
	m.compute_topological();

	int poin = 8;
	cout << "Els surrounding point " << poin << endl;
	for(int i = m.gesup_p(poin); i < m.gesup_p(poin+1); i++)
		cout << " " << m.gesup(i);
	cout << endl;

	cout << "Points surrounding point " << poin << endl;
	for(int i = m.gpsup_p(poin); i < m.gpsup_p(poin+1); i++)
		cout << " " << m.gpsup(i);
	cout << endl;

	/*int iedge = 92;
	cout << "Data for edge " << iedge << ":" << endl << m.gintedge(iedge,0) << " " << m.gintedge(iedge,1) << endl;
	for(int i = 0; i < m.gelsedsize(iedge); i++)
		cout << " " << m.gelsed(iedge,i);
	cout << endl;*/

	int el = 2;
	cout << "elements surrounding element " << el << endl;
	for(int i = 0; i < m.gnfael(el); i++)
		cout << " " << m.gesuel(el,i);
	cout << endl;

	/*int fac = 26;
	cout << "intfac for face " << fac << endl;
	for(int fac = 0; fac < m.gnaface(); fac++)
	{
		for(int i = 0; i < m.gnnofa()+2; i++)
			cout << " " << m.gintfac(fac,i);
		cout << endl;
	}*/
	m.printmeshstats();

	UMesh2dh mq = m.convertLinearToQuadratic();
	mq.writeGmsh2("testhybrid_quadratic.msh");

	cout << endl;
	return 0;
}
