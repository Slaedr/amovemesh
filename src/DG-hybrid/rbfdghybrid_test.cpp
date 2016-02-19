#include "arbfdghybrid.hpp"

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
	double tolerance, suprad; int nlayers, maxiter;
	ifstream fin(argv[1]);
	fin >> dum; fin >> inmeshlin;
	fin >> dum; fin >> inmeshquad;
	fin >> dum; fin >> outmesh;
	fin >> dum; fin >> nlayers;
	fin >> dum; fin >> suprad;
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
	/*bounmotion(20,0) = 0.25; bounmotion(20,1) = 0.25;
	bounmotion(22,0) = 0; bounmotion(22,1) = 0.3;
	bounmotion(23,0) = -0.25; bounmotion(23,1) = 0.25;*/
	bounmotion(11,0) = 0.01; bounmotion(11,1) = 0.1;
	bounmotion(12,0) = 0; bounmotion(12,1) = 0.15;
	bounmotion(13,0) = -0.01; bounmotion(13,1) = 0.1;

	DGhybrid dgh(&m,&mq,&bounmotion, nlayers, suprad, tolerance, maxiter, solver);
	dgh.compute_backmesh_points();
	dgh.generate_backmesh_and_compute_displacements();
	dgh.movemesh();

	mq.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
