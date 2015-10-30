#include "amatrix2.hpp"
#include "alinalg.hpp"

using namespace amat;
using namespace std;

int main()
{
	Matrix<double> A(4,4);
	double Ad[] = {5/2.0,2.5,4,4/5.0, 1.0/4,42,3.1,0, 4.2,11,0,3, 10,3.2,0,2};
	A.setdata(Ad, 16);
	A.mprint();
	
	Matrix<double> b(4,1);
	double bd[] = {0.5,4.0/3,3.0/2,1.0};
	b.setdata(bd, 4);
	b.mprint();
	
	double& var = A(1,1);
	var = 6;
	cout << "** " << A(1,1) << endl;
	
	Matrix<double> x = gausselim(A,b);
	
	x.mprint();
	cout << endl;
	
	return 0;
}
