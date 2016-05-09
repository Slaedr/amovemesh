#include "acurvedmeshgen2dh_sr.hpp"
#include <ctime>
#include <limits>

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	//string confile = "curvedmeshgen2dh.control";
	if(argc < 2) {
		cout << "Insufficient arguments given!\n";
			return 0;
	}

	string confile(argv[1]);
	ifstream conf(confile);
	string dum, str_linearmesh, str_qmesh, str_curvedmesh, rbf_solver;
	int rbfchoice, rbfmaxiter, splmaxiter, nrbfsteps, nsplineparts, nflags;
	double supportradius, rbftol, spltol, cornerangle;
	vector<vector<int>> splineflags;
	vector<int> vec;

	conf >> dum; conf >> str_linearmesh;
	conf >> dum; conf >> str_qmesh;
	conf >> dum; conf >> str_curvedmesh;
	conf >> dum; conf >> cornerangle;
	conf >> dum; conf >> spltol;
	conf >> dum; conf >> splmaxiter;
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
	
	UMesh2dh m, mq;
	m.readGmsh2(str_linearmesh,2);
	mq.readGmsh2(str_qmesh,2);

	m.compute_boundary_points();
	//m.compute_topological();

	RBFmove rmove;
	Curvedmeshgen2d cu;
	cu.setup(&m, &mq, &rmove, nsplineparts, splineflags, PI/180.0*cornerangle, spltol, splmaxiter, rbftol, rbfmaxiter, rbfchoice, supportradius, nrbfsteps, rbf_solver);
	cu.compute_boundary_displacements();
	
	clock_t begin = clock();
	cu.generate_curved_mesh();
	clock_t end = clock() - begin;
	cout << "Time taken by RBF is " << (double(end))/CLOCKS_PER_SEC << endl;

	mq.writeGmsh2(str_curvedmesh);

	//cout << "Machine epsilon " << numeric_limits<double>::epsilon() << endl;
	cout << endl;
	return 0;
}
