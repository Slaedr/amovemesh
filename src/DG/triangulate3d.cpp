/** @file triangulate.cpp
 * @brief Construct Delaunay tetrahedralization of a set of 3D points.
 * @date January 8, 2015
 */

#include <abowyerwatson3d.hpp>

using namespace std;
using namespace amat;

void read_points(string finname, Delaunay3d& d3d)
{
	string dum; int indum; char chdum;
	
	ifstream fin(finname);
	getline(fin,dum);
	fin >> d3d.npoints;
	fin >> chdum;
	getline(fin,dum);

	d3d.points.setup(d3d.npoints,3);
	d3d.nodes.reserve(d3d.npoints+3);

	for(int i = 0; i < d3d.npoints; i++)
	{
		fin >> indum;
		for(int j = 0; j < 3; j++)
			fin >> d3d.points(i,j);
	}

	fin.close();
}

int main(int argc, char* argv[])
{
	Delaunay3d d3;
	read_points("../../input/delau.domn.points", d3);
	d3.bowyer_watson();
	d3.writeGmsh2("../../output/dg.msh");
	d3.compute_jacobians();
	bool negjac = d3.detect_negative_jacobians();
	cout << endl;
	return 0;
}
