/** @file testd.cpp
 * @author Aditya Kashi
 */

#include <abowyerwatson3d.hpp>

using namespace amat;
using namespace std;

int main()
{
	int np = 6;
	Matrix<double> pts(np,3);
	pts(0,0) = 0; pts(0,1) = 0; pts(0,2) = 0;
	pts(1,0) = 1; pts(1,1) = 0; pts(1,2) = 0;
	pts(2,0) = 0; pts(2,1) = 1; pts(2,2) = 0;
	pts(3,0) = 0; pts(3,1) = 0; pts(3,2) = 1;
	pts(4,0) = 1; pts(4,1) = 1; pts(4,2) = 1;
	pts(5,0) = 1; pts(5,1) = 3; pts(5,2) = 3;

	Delaunay3d d(&pts, np);
	d.bowyer_watson();
	d.writeGmsh2("tryd.msh");
	cout << endl;

	return 0;
}
