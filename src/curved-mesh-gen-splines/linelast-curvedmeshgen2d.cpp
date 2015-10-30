#include "alinelast-curvedmeshgen2d.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "linelast-curvedmeshgen2d.control";
	ifstream conf(confile);
	string dum, str_linearmesh, str_qmesh, str_curvedmesh;
	int maxiter, nsplineparts, nflags;
	double ym, pr, tol, cornerangle;
	vector<vector<int>> splineflags;
	vector<int> vec;

	conf >> dum; conf >> str_linearmesh;
	conf >> dum; conf >> str_qmesh;
	conf >> dum; conf >> str_curvedmesh;
	conf >> dum; conf >> cornerangle;
	conf >> dum; conf >> ym;
	conf >> dum; conf >> pr;
	conf >> dum; conf >> tol;
	conf >> dum; conf >> maxiter;
	conf >> dum; conf >> nsplineparts;
	conf >> dum;
	splineflags.reserve(nsplineparts);
	for(int i = 0; i < nsplineparts; i++)
	{
		conf >> nflags;
		vec.resize(nflags);
		for(int j = 0; j < nflags; j++)
			conf >> vec[j];
		splineflags.push_back(vec);
	}

	conf.close();

	cout << "Number of parts to reconstruct = " << nsplineparts << endl;
	for(int i = 0; i < nsplineparts; i++)
		cout << " Number of markers in part " << i << " = " << splineflags[i].size() << endl;
	
	UMesh2d m, mq;
	m.readGmsh2(str_linearmesh,2);
	mq.readGmsh2(str_qmesh,2);

	m.compute_boundary_points();
	//m.compute_topological();

	LinElastP2 rmove;
	Curvedmeshgen2d cu;
	cu.setup(&m, &mq, &rmove, nsplineparts, splineflags, PI/180.0*cornerangle, tol, maxiter, ym, pr);
	cu.compute_boundary_displacements();
	cu.generate_curved_mesh();

	mq.writeGmsh2(str_curvedmesh);

	cout << endl;
	return 0;
}
