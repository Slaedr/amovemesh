#include <amesh2d.hpp>
#include "acurvedmeshgen-rbf.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string inmeshg = "rans_bump_liu_2nd_J0.8.msh";
	UMesh2d mg;
	mg.readGmsh2(inmeshg,2);

	string inmesh2 = "rans_bump_liu_base.msh";
	UMesh2d m2;
	m2.readGmsh2(inmesh2,2);

	Matrix<int> fixed_b_flags(2,1);
	fixed_b_flags(0) = 9; fixed_b_flags(1) = 1;

	CurvedMeshGeneration cmg(&m2, fixed_b_flags, 1, 2, 0.002);
	cmg.compute_boundary_displacement_comp(mg);
	cmg.generate();
	m2.writeGmsh2("rans_bump_aditya-0.002.msh");

	double norm = 0;

	for(int i = 0; i < mg.gnpoin(); i++)
	{
		norm += pow(sqrt(mg.gcoords(i,0)*mg.gcoords(i,0) + mg.gcoords(i,1)*mg.gcoords(i,1)) - sqrt(m2.gcoords(i,0)*m2.gcoords(i,0) + m2.gcoords(i,1)*m2.gcoords(i,1)),2);
	}
	norm = sqrt(norm);
	cout << "Norm is " << norm  << endl;

	cout << endl;
	return 0;
}
