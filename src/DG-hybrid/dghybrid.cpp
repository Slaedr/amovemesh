#include "adghybrid.hpp"

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	UMesh2dh m;
	m.readGmsh2("../../input/hairfoil.msh", 2);
	m.compute_topological();

	UMesh2d bm;

	DGhybrid dgh(m,bm,30);
	dgh.compute_backmesh_points();

	cout << endl;
	return 0;
}
