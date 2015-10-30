#include <cmath>
#include <amesh2.hpp>
#include <aoutput.hpp>
#include <sstream>

#include "amovemesh2d.hpp"

using namespace std;
using namespace acfd;

const double pi = 3.14159265;

int main()
{
    ifstream control("lineal.control");
	
	string meshfile, initjacfile, finjacfile, movedmeshfile, solver, anglestr;
	
	getline(control, meshfile);
	cout << meshfile << endl;
	getline(control, initjacfile);
	getline(control, finjacfile);
	getline(control, movedmeshfile);
	getline(control, solver);
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
	//m.detect_negative_jacobians(du1);
	
	writeMeshToVtu("mymesh.vtu",m);
	cout << "Going to MRotation2d..\n";
	MRotation2d r(&m, angle*pi/180, 0, 0);
	Matrix<double> bc = r.rhsvect_rotate();
	//Matrix<double> bc(m.gnpoin(), 2, ROWMAJOR);
	//bc.zeros();
	
	Matrix<double> bcx = bc.col(0);
	Matrix<double> bcy = bc.col(1);
	//m.movemesh_basicspring(bc.col(0),bc.col(1),solver,1e-5,3000);
	//writeMeshToVtu(movedmeshfile+"-"+solver+anglestr+".vtu",m);
	m.movemesh_lineal(&bcx,&bcy,solver,1e-5,5500);
	writeMeshToVtu(movedmeshfile+"-lineal-"+solver+anglestr+".vtu",m);
	
	m.detect_negative_jacobians(du2);
	
	mfile.close();
	//du1.close(); 
	du2.close();
	cout << endl;
    return 0;
}
