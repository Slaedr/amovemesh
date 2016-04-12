#include <amesh2d.hpp>
#include "arbf.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string inmesh = "rans_bump_liu_2nd_J0.8.msh";
	UMesh2d mg;
	mg.readGmsh2(inmesh,2);

	UMesh2d m2 = mg;
	Matrix<double> newcoords(m2.gnpoin(), m2.gndim());
	newcoords = *(m2.getcoords());

	for(int ielem = 0; ielem < m2.gnelem(); ielem++)
	{
		for(int inode = 4; inode < 7; inode++)
		{
			for(int idim = 0; idim < 2; idim++)
				newcoords(m2.ginpoel(ielem,inode),idim) = (m2.gcoords(m2.ginpoel(ielem,inode-4),idim) + m2.gcoords(m2.ginpoel(ielem,inode-3),idim))/2.0;
		}

		for(int idim = 0; idim < 2; idim++)
			newcoords(m2.ginpoel(ielem,7),idim) = (m2.gcoords(m2.ginpoel(ielem,3),idim) + m2.gcoords(m2.ginpoel(ielem,0),idim))/2.0;

		for(int idim = 0; idim < 2; idim++)
			newcoords(m2.ginpoel(ielem,8),idim) = (m2.gcoords(m2.ginpoel(ielem,0),idim) + m2.gcoords(m2.ginpoel(ielem,1),idim) + m2.gcoords(m2.ginpoel(ielem,2),idim) + m2.gcoords(m2.ginpoel(ielem,3),idim))/4.0;
	}

	m2.setcoords(&newcoords);

	string outmesh = "rans_bump_liu_base.msh";
	m2.writeGmsh2(outmesh);

	cout << endl;
}
