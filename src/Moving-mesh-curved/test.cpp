#include "amesh_c.hpp"
#include "amesh_curved.hpp"
#include <aoutput.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
    string meshn = "testmesh2.domn";
    string origmeshfile = "origmesh-aspect.vtu";
    string refmeshfile = "refmesh-aspect.vtu";
    string movemeshfile = "movemesh-aspect.vtu";
    string cmeshfile = "curvedmesh2.msh";
    string solver = "gaussseidel";
    ifstream meshfile(meshn);
    UTriMesh m(meshfile);
    //writeMeshToVtu(origmeshfile,m);
    m.refine_mesh();
    //writeMeshToVtu(refmeshfile,m);
    //m.esuel.mprint();
    Matrix<double> xb(m.gnpoin(),1), yb(m.gnpoin(),1);
    xb.zeros(); yb.zeros();
    //construct xb and yb - the fixed nodal displacements
    yb(12) = yb(15) = 0.7;
    Matrix<int> extraflags(m.gnpoin(),1); extraflags.zeros();
    extraflags(4) = extraflags(7) = 1;

    //m.movemesh_basicspring(xb,yb,extraflags,solver);
    m.movemesh_farhat(&xb,&yb,extraflags,solver);

    //writeMeshToVtu(movemeshfile,m);

    UTriMeshCurved mc(&m);
    mc.writeQuadraticGmsh2(cmeshfile);

    cout << endl;
    return 0;
}
