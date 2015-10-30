/** Computes skew-size metric of a mesh.
	The mesh file name is passed as the first command line argument, the area of the domain is passed as the second, and the output file name is the third argument.
*/

#ifndef _GLIBCXX_CSTDLIB
#include <cstdlib>
#endif

#include <amesh2d.hpp>
#include <aoutput.hpp>

using namespace amat;
using namespace acfd;

int main(int argc, char* argv[])
{
	if(argc < 4) {
		cout << "Insufficient arguments given!";
		return 0;
	}

	string mname(argv[1]);
	double area = atof(argv[2]);
	string outname(argv[3]);
	
	UMesh2d m;
	m.readGmsh2(mname,2);

	// get reference area
	double w = area/m.gnelem();

	m.compute_metric_quantities();
	Matrix<double> skew(m.gnelem(),1);
	m.linearmetric_skew(&skew);
	Matrix<double> skewsize(m.gnelem(),1);
	m.linearmetric_skewsize(w,&skewsize);

	Matrix<double> metrics(m.gnelem(),2);
	for(int iel = 0; iel < m.gnelem(); iel++)
	{
		metrics(iel,0) = skew(iel,0);
		metrics(iel,1) = skewsize(iel,0);
	}

	// output
	Matrix<double> vect;
	string scaname[] = {"skew_metric","sizeskew_metric"};
	writeScalarsVectorToVtu_CellData(outname,m,metrics,scaname,vect,"none");
	return 0;
}
