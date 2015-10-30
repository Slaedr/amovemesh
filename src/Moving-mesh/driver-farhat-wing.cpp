#include <cmath>
//#include <sstream>

#include <arotation2db.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

const double pi = 3.14159265;

int main()
{
    ifstream control("farhat-wing.control");

	string meshfile, initjacfile, finjacfile, movedmeshfile, solver, anglestr;

	getline(control, meshfile);
	cout << meshfile << endl;
	getline(control, initjacfile);
	getline(control, finjacfile);
	getline(control, movedmeshfile);
	getline(control, anglestr);

	double angle = stod(anglestr);
	cout << "Angle = " << angle << endl;
	cout << "Solver: " << solver << endl;
	ifstream mfile(meshfile);
	//ofstream du1(initjacfile);
	ofstream du2(finjacfile);
	cout << "Finished reading control file. Opened mesh file and jacobian file(s)\n";

    UTriMesh m(mfile);
	m.compute_jacobians();
	m.compute_face_data();
	m.allocate_edge_angles();
	m.compute_edge_angles();
    //m.setbflags2();
	//m.detect_negative_jacobians(du1);

	cout << "Going to MRotation2d..\n";
    Matrix<int> n_rot(1,1); n_rot(0) = 3;		// flag of body to be rotated
	MRotation2d r(&m, angle*pi/180, m.gcoords(420,0), m.gcoords(420,1), n_rot);
	Matrix<double> bc = r.rhsvect_rotate();
	//Matrix<double> bc(m.gnpoin(), 2, ROWMAJOR);
	//bc.zeros();
    /*ofstream bcfile("bdata.dat");
    bc.fprint(bcfile);
    bcfile.close();*/

	Matrix<double> bcx = bc.col(0);
	Matrix<double> bcy = bc.col(1);
	//m.movemesh_basicspring(bc.col(0),bc.col(1),solver,1e-5,3000);
	//writeMeshToVtu(movedmeshfile+"-"+solver+anglestr+".vtu",m);
    m.set_boundary_motion(&bcx,&bcy);
	m.movemesh_farhat(1e-6,4000);
    m.writeGmsh2(movedmeshfile+"-farhat-"+anglestr+".msh");

	m.detect_negative_jacobians(du2);

	mfile.close();
	//du1.close();
	du2.close();
	cout << endl;
    return 0;
}
