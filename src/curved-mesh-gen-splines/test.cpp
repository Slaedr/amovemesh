#include <alinalg.hpp>

using namespace std;
using namespace amat;

int main()
{
	SpMatrix A(5,5);
	A.set(0,0,1);
	A.set(1,1,1.5);
	A.set(2,2,2.5);
	A.set(3,3,1.1);
	A.set(4,4,0.6);
	
	Matrix<double> b(5,1);
	b.zeros();
	b(0) = 2.5;
	
	Matrix<double> x(5,1);
	
	superLU_solve(&A, &b, &x);

	cout << endl;
	return 0;
}
