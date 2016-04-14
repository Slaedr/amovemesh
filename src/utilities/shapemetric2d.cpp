/** @brief Computes shape metric of a mesh.
 * The mesh file name is passed as the first command line argument, the reference area is passed as the second, and the output file name is the third argument.
 * @author Aditya Kashi
 */

#ifndef _GLIBCXX_CSTDLIB
#include <cstdlib>
#endif

#include <amesh2d.hpp>
#include <aoutput.hpp>

using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	if(argc < 4) {
		cout << "Insufficient arguments given!\n";
		cout << "Required <mesh file name> <reference area> <output file>\n";
		return 0;
	}

	string mname(argv[1]);
	double area = atof(argv[2]);
	string outname(argv[3]);
	
	UMesh2d m;
	m.readGmsh2(mname,2);
	m.compute_jacobians();

	// get reference area
	double w = m.garea()/m.gnelem();

	m.compute_metric_quantities();
	Matrix<double> shape(m.gnelem(),1);
	m.linearmetric_shape(&shape);
	Matrix<double> shapesize(m.gnelem(),1);
	m.linearmetric_shapesize(&shapesize);

	Matrix<double> metrics(m.gnelem(),2);
	for(int iel = 0; iel < m.gnelem(); iel++)
	{
		metrics(iel,0) = shape(iel,0);
		metrics(iel,1) = shapesize(iel,0);
	}

	// output
	Matrix<double> vect;
	string scaname[] = {"shape_metric","sizeshape_metric"};
	writeScalarsVectorToVtu_CellData(outname,m,metrics,scaname,vect,"none");
	return 0;
}
