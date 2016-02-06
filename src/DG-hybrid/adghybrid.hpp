/** \brief Hybrid DG-elasticity method for mesh movement
 */

#ifndef __ADGM_H
#include <adgm.hpp>
#endif

using namespace amat;
using namespace acfd;
using namespace std;

class DGhybrid
{
	/// input mesh
	UMesh2dh* m;
	/// coarse mesh generated from input mesh
	UMesh2d* bm;
	
	DGmove dgm;
	LinElast_P2 linm;
};
