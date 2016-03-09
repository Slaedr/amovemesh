#include "acurvedmeshgen2d.hpp"
#include <ctime>

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Please provide a control file!\n";
		return -1;
	}

	ifstream fin(argv[1]);
	if(!fin) {
		cout << "Could not find control file.\n";
		return -2;
	}

	string dum, infile, outfile, solver;
	int maxiter, rbf_choice;
	double tol, sup_rad;

	fin >> dum; fin >> infile;
	fin >> dum; fin >> outfile;
	fin >> dum; fin >> rbf_choice;
	fin >> dum; fin >> sup_rad;
	fin >> dum; fin >> solver;
	fin >> dum; fin >> tol;
	fin >> dum; fin >> maxiter;

	fin.close();

	cout << "RBF choice = " << rbf_choice << ", support radius = " << sup_rad << ", solver = " << solver << ", tolerance = " << tol << ", Max iterations = " << maxiter << ".\n";

	UMesh2dh m;
	m.readGmsh2(infile, 2);

	cout << "CurvedMeshGen starting" << endl;
	CurvedMeshGen cmg(&m, rbf_choice, sup_rad, tol, maxiter, solver);

	clock_t begin = clock();
	cmg.generateCurvedMesh();
	clock_t end = clock() - begin;
	cout << "Time taken by mesh movement is " << (double(end))/CLOCKS_PER_SEC << endl;
	
	m.writeGmsh2(outfile);

	cout << endl;
	return 0;
}
