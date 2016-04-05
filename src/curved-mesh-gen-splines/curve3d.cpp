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

	// compute norm of error for unit ball case
	vector<amc_real> errors(mq.gnbpoin(),0);
	int i,j,k,idim;
	amc_real errnorm = 0, errmax = 0;
	k = 0;
	for(i = 0; i < mq.gnpoin(); i++)
		if(mq.gflag_bpoin(i) == 1)
		{
			for(idim = 0; idim < 3; idim++)
				errors[k] += mq.gcoords(i,idim)*mq.gcoords(i,idim);
			k++;
		}

	cout << setprecision(20);
	for(i = 0; i < mq.gnbpoin(); i++)
	{
		//cout << " " << errors[i] << endl;
		errors[i] = fabs(sqrt(errors[i]) - 1.0);
		
		if(errors[i] > errmax) errmax = errors[i];
		
		errnorm += errors[i]*errors[i];
	}
	errnorm = sqrt(errnorm);
	cout << "The error l2 norm is " << errnorm << ", and the max norm is " << errmax << endl;
	
	std::cout << std::endl;

	std::cout << "Computing the error for straight mesh... \n";
	string stmesh = "../../input/ball-coarse_p2.msh";
	mq.readGmsh2(stmesh, 3);
	errnorm = 0; errmax = 0; k = 0;
	errors.assign(mq.gnbpoin(),0);
	for(i = 0; i < mq.gnpoin(); i++)
		if(mq.gflag_bpoin(i) == 1)
		{
			for(idim = 0; idim < 3; idim++)
				errors[k] += mq.gcoords(i,idim)*mq.gcoords(i,idim);
			k++;
		}

	cout << setprecision(20);
	for(i = 0; i < mq.gnbpoin(); i++)
	{
		//cout << " " << errors[i] << endl;
		errors[i] = fabs(sqrt(errors[i]) - 1.0);
		
		if(errors[i] > errmax) errmax = errors[i];
		
		errnorm += errors[i]*errors[i];
	}
	errnorm = sqrt(errnorm);
	cout << "The error l2 norm for straight mesh is " << errnorm << ", and the max norm is " << errmax << endl;

	return 0;
}
