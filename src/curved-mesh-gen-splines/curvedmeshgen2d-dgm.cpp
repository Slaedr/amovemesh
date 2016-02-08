#include "acurvedmeshgen2d-dgm.hpp"
#include <ctime>

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Error: Provide control file name.\n";
		return 0;
	}

	string confile = argv[1];
	ifstream conf(confile);
	string dum, str_linearmesh, str_qmesh, str_curvedmesh;
	int nsplineparts, nflags, maxiter;
	double tol, cornerangle;
	vector<vector<int>> splineflags;
	vector<int> vec;

	conf >> dum; conf >> str_linearmesh;
	conf >> dum; conf >> str_qmesh;
	conf >> dum; conf >> str_curvedmesh;
	conf >> dum; conf >> cornerangle;
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

	DGmove dmove;
	Curvedmeshgen2d cu;
	cu.setup(&m, &mq, &dmove, nsplineparts, splineflags, PI/180.0*cornerangle, tol, maxiter);
	cu.compute_boundary_displacements();
	
	clock_t begin = clock();
	cu.generate_curved_mesh();
	clock_t end = clock() - begin;
	cout << "Time taken by DGM is " << (double(end))/CLOCKS_PER_SEC << endl;

	mq.writeGmsh2(str_curvedmesh);

	cout << endl;
	return 0;
}
