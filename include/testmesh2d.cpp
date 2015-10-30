#include "amesh2d.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	UMesh2d m;
	m.readGmsh2("tinquad.msh",2);
	m.compute_boundary_points();
}
