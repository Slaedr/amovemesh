#include "amesh3d.hpp"

using namespace amat;
using namespace acfd;
using namespace std;

int main()
{
	string mfilename = "tryhex.msh";
	UMesh m;
	m.readGmsh2(mfilename,3);
	cout << endl;
	m.compute_topological();

	/*int poin = 21;
	cout << "Els surrounding point " << poin << endl;
	for(int i = m.gesup_p(poin); i < m.gesup_p(poin+1); i++)
		cout << " " << m.gesup(i);
	cout << endl;

	cout << "Points surrounding point " << poin << endl;
	for(int i = 0; i < m.gpsupsize(poin); i++)
		cout << " " << m.gpsup(poin,i);
	cout << endl;

	int iedge = 92;
	cout << "Data for edge " << iedge << ":" << endl << m.gintedge(iedge,0) << " " << m.gintedge(iedge,1) << endl;
	for(int i = 0; i < m.gelsedsize(iedge); i++)
		cout << " " << m.gelsed(iedge,i);
	cout << endl;

	int el = 2;
	cout << "elements surrounding element " << el << endl;
	for(int i = 0; i < m.gnfael(); i++)
		cout << " " << m.gesuel(el,i);
	cout << endl;

	int fac = 74;
	cout << "intfac for face " << fac << endl;
	for(int i = 0; i < m.gnnofa()+2; i++)
		cout << " " << m.gintfac(fac,i);
	cout << endl;*/

	UMesh mq = m.convertLinearToQuadratic();
	mq.printmeshstats();
	mq.writeGmsh2("trycurvedhex.msh");

	cout << endl;
	return 0;
}
