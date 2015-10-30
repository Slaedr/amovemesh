#include <ageometry.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	UMesh2d m;

	string confile = "calc-bmotion.control", dum, mfile, outfile;
	int n_parts, numf, maxiter;
	double angle, tol;
	vector<vector<int>> f2r;
	vector<int> vec;
	ifstream cfile(confile);
	cfile >> dum; cfile >> mfile;
	cfile >> dum; cfile >> outfile;
	cfile >> dum; cfile >> angle;
	cfile >> dum; cfile >> n_parts;
	cfile >> dum;
	f2r.reserve(n_parts);
	for(int i = 0; i < n_parts; i++) {
		cfile >> numf;
		vec.resize(numf);
		for(int j = 0; j < numf; j++)
			cfile >> vec[j];
		f2r.push_back(vec);
	}
	cfile >> dum; cfile >> tol;
	cfile >> dum; cfile >> maxiter;
	cfile.close();

	cout << "Number of parts to reconstruct = " << n_parts << endl;
	for(int i = 0; i < n_parts; i++)
		cout << " Number of markers in part " << i << " = " << f2r[i].size() << endl;

	m.readGmsh2(mfile,2);
	m.compute_boundary_points();

	BoundaryReconstruction2d br;
	br.setup(&m, n_parts, f2r, PI/180.0*angle);
	br.preprocess();
	br.detect_corners();
	br.split_parts();
	br.compute_splines(tol, maxiter);

	ofstream fout(outfile);
	for(int iface = 0; iface < m.gnface(); iface++)
	{
		for(int i = 0; i < f2r.size(); i++)
			for(int j = 0; j < f2r[i].size(); j++)
				if(m.gbface(iface,m.gnnofa()) == f2r[i][j])
				{
					for(double u = 0; u <= 1; u += 0.1)
						fout << br.getcoords(iface,0,u) << " " << br.getcoords(iface,1,u) << '\n';
				}
	}
	fout.close();

	cout << "Spline points written to file " << outfile << endl;

	cout << endl;
	return 0;
}
