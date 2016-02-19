#include "acmg2dh_rbfdg.hpp"

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Give a control file name!\n";
		return 0;
	}
	
	string confile(argv[1]);
	
	ifstream conf(confile);
	string dum, str_linearmesh, str_qmesh, str_curvedmesh, rbf_solver;
	int rbfchoice, splmaxiter, rbfmaxiter, nrbfsteps, nsplineparts, nflags, nlayers;
	double supportradius, spltol, rbftol, cornerangle;
	vector<vector<int>> splineflags;
	vector<int> vec;

	conf >> dum; conf >> str_linearmesh;
	conf >> dum; conf >> str_qmesh;
	conf >> dum; conf >> str_curvedmesh;
	conf >> dum; conf >> cornerangle;
	conf >> dum; conf >> spltol;
	conf >> dum; conf >> splmaxiter;
	conf >> dum; conf >> nlayers;
	conf >> dum; conf >> rbfchoice;
	conf >> dum; conf >> supportradius;
	conf >> dum; conf >> nrbfsteps;
	conf >> dum; conf >> rbftol;
	conf >> dum; conf >> rbfmaxiter;
	conf >> dum; conf >> rbf_solver;
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
	
	UMesh2dh m,mq;
	mq.readGmsh2(str_qmesh, 2);
	
	m.readGmsh2(str_linearmesh, 2);
	m.compute_topological();
	m.compute_boundary_points();

	Curvedmeshgen2d cmg;
	cmg.setup(&m, &mq, nsplineparts, splineflags, cornerangle, spltol, splmaxiter, rbftol, rbfmaxiter, rbf_solver, supportradius, nlayers);
	cmg.compute_boundary_displacements();
	cmg.generate_curved_mesh();
	
	mq.writeGmsh2(str_curvedmesh);

	cout << endl;
	return 0;
}
