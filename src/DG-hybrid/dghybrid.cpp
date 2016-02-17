#include "adghybrid.hpp"

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
	double tolerance; int nlayers, maxiter;
	ifstream fin(argv[1]);
	fin >> dum; fin >> inmeshlin;
	fin >> dum; fin >> inmeshquad;
	fin >> dum; fin >> outmesh;
	fin >> dum; fin >> nlayers;
	fin >> dum; fin >> tolerance;
	fin >> dum; fin >> maxiter;
	fin >> dum; fin >> solver;
	fin.close();
	
	UMesh2dh m,mq;
	mq.readGmsh2(inmeshquad, 2);
	
	m.readGmsh2(inmeshlin, 2);
	m.compute_topological();
	
	Matrix<double> bounmotion(mq.gnpoin(), mq.gndim());
	bounmotion.zeros();
	bounmotion(20,0) = 0.25; bounmotion(20,1) = 0.25;
	bounmotion(22,0) = 0; bounmotion(22,1) = 0.3;
	bounmotion(23,0) = -0.25; bounmotion(23,1) = 0.25;

	DGhybrid dgh(&m,&mq,&bounmotion, nlayers, 1e9, 0.4, tolerance, maxiter, solver);
	dgh.compute_backmesh_points();
	dgh.generate_backmesh_and_compute_displacements();
	//dgh.movemesh();

	cout << endl;
	return 0;
}
