/** @file testd.cpp
 * @author Aditya Kashi
 */

#include <abowyerwatson3d.hpp>

using namespace amat;
using namespace std;

int main()
{
	int np = 12;
	Matrix<double> pts(np,3);
	/*pts(0,0) = 0; pts(0,1) = 0; pts(0,2) = 0;
	pts(1,0) = 1; pts(1,1) = 0; pts(1,2) = 0;
	pts(2,0) = 0; pts(2,1) = 1; pts(2,2) = 0;
	pts(3,0) = 0; pts(3,1) = 0; pts(3,2) = 1;
	pts(4,0) = 1; pts(4,1) = 1; pts(4,2) = 1;
	pts(5,0) = 1; pts(5,1) = 3; pts(5,2) = 3;*/
	
	pts(0,0) = 0; pts(0,1) = 0; pts(0,2) = 0;
	pts(1,0) = 0; pts(1,1) = 1; pts(1,2) = 0;
	pts(2,0) = 1; pts(2,1) = 0; pts(2,2) = 0;
	pts(3,0) = 1; pts(3,1) = 1; pts(3,2) = 0;
	pts(4,0) = 0; pts(4,1) = 0; pts(4,2) = 1;
	pts(5,0) = 0; pts(5,1) = 1; pts(5,2) = 1;
	pts(6,0) = 1; pts(6,1) = 0; pts(6,2) = 1;
	pts(7,0) = 1; pts(7,1) = 1; pts(7,2) = 1;
	pts(8,0) = 1; pts(8,1) = 1.1; pts(8,2) = 1.5;
	pts(9,0) = 0; pts(9,1) = 1.01; pts(9,2) = 1.5;
	pts(10,0) = 1; pts(10,1) = 0; pts(10,2) = 1.5;
	pts(11,0) = 0; pts(11,1) = 0; pts(11,2) = 1.5;

	Delaunay3d d(&pts, np);
	d.bowyer_watson();
	d.writeGmsh2("tryd.msh");
	d.compute_jacobians();
	bool decneg = d.detect_negative_jacobians();
	cout << endl;

	return 0;
}
