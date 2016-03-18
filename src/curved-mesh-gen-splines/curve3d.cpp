#include <ageometry3d.hpp>

using namespace std;
using namespace amc;

int main()
{
	string meshfile = "../../input/ball-coarse.msh";

	UMesh m;
	m.readGmsh2(meshfile,3);

	VertexCenteredBoundaryReconstruction br(&m,2);
}
