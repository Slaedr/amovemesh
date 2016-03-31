#include "acurvedmeshgen3d.hpp"

using namespace std;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Please specify a control file.\n";
		return -1;
	}
#ifdef DEBUG
	cout << "DEBUG!\n";
#endif
	string confile = argv[1], linmesh, cmesh, solver, brtype, stenciltype, dum;
	amc_real tol, angle_limit, suprad;
	int maxiter, rbf_choice, rbf_steps;
	ifstream conf(confile);

	conf >> dum; conf >> linmesh;
	conf >> dum; conf >> cmesh;
	conf >> dum; conf >> brtype;
	conf >> dum; conf >> stenciltype;
	conf >> dum; conf >> angle_limit;
	conf >> dum; conf >> rbf_choice;
	conf >> dum; conf >> suprad;
	conf >> dum; conf >> rbf_steps;
	conf >> dum; conf >> tol;
	conf >> dum; conf >> maxiter;
	conf >> dum; conf >> solver;
	
	conf.close();

	cout << "amc: Generating curved mesh with " << brtype << "-WALF, " << stenciltype << " stencil, " << angle_limit << ", " 
		<< rbf_choice << ", " << suprad << ", " << tol << ", " << maxiter << ", " << solver << endl;

	UMesh m;
	m.readGmsh2(linmesh,3);
	m.compute_topological();
	m.compute_boundary_topological();

	UMesh mq = m.convertLinearToQuadratic();
	
	CurvedMeshGen cmg;
	cmg.setup(&m, &mq, brtype, stenciltype, angle_limit, tol, maxiter, rbf_choice, suprad, rbf_steps, solver);
	cmg.compute_boundary_displacements();
	cmg.generate_curved_mesh();
	mq.writeGmsh2(cmesh);

	std::cout << std::endl;

	return 0;
}
