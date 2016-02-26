#include "acmg2dh_linedg.hpp"

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
	
	string inmeshlin, inmeshquad, outmesh, dum, solver;
	double spl_tol, le_tol, suprad, angle_threshold, young, pratio;
	int nlayers, spl_maxiter, le_maxiter, num_parts, nflags;
	vector<int> vec;
	vector<vector<int>> splineflags;

	ifstream fin(argv[1]);
	fin >> dum; fin >> inmeshlin;
	fin >> dum; fin >> inmeshquad;
	fin >> dum; fin >> outmesh;
	fin >> dum; fin >> angle_threshold;
	fin >> dum; fin >> spl_tol;
	fin >> dum; fin >> spl_maxiter;
	fin >> dum; fin >> nlayers;
	fin >> dum; fin >> young;
	fin >> dum; fin >> pratio;
	fin >> dum; fin >> le_tol;
	fin >> dum; fin >> le_maxiter;
	fin >> dum; fin >> solver;
	fin >> dum; fin >> num_parts;
	fin >> dum;
	splineflags.reserve(num_parts);
	for(int i = 0; i < num_parts; i++)
	{
		fin >> nflags;
		vec.resize(nflags);
		for(int j = 0; j < nflags; j++)
		{
			fin >> vec[j];
		}
		splineflags.push_back(vec);
	}
	fin.close();
	
	UMesh2dh m,mq;
	mq.readGmsh2(inmeshquad, 2);
	
	m.readGmsh2(inmeshlin, 2);
	m.compute_topological();
	m.compute_boundary_points();
	
	Curvedmeshgen2d cur;
	cur.setup(&m, &mq, num_parts, splineflags, angle_threshold, spl_tol, spl_maxiter, le_tol, le_maxiter, solver, young, pratio, nlayers);
	cur.compute_boundary_displacements();
	cur.generate_curved_mesh();

	mq.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
