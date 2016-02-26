#include "alinalg.hpp"

using namespace amat;
using namespace std;

int main()
{
	SpMatrix A(5,5);
	A.set(0,0, 1);
	A.set(1,1, 1);
	A.set(2,2, 1);
	A.set(3,3, 1);
	A.set(4,4, 1);
	
	Matrix<double> rhs(5,2);
	rhs.zeros();
	rhs(1,0) = 2.5;
	rhs(2,0) = 1.5;
	rhs(1,1) = 4.5;
	rhs(2,1) = 3.5;

	Matrix<double> x = gausselim(A,rhs);

	x.mprint();

	return 0;
}
