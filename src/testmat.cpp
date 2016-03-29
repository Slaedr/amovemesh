#include "amatrix2.hpp"

using namespace amat;
using namespace std;

int main()
{
	int m = 4, n = 3;
	Matrix<double> A(4,3);
	double Ad[] = {1,-1,4, 1,4,-2, 1,4,2, 1,-1,0};
	A.setdata(Ad, 12);
	A.mprint();
	cout << "Matrix printed above\n";
	cout << "1-norm " << A.matrixNorm_1() << endl;
	
	/*Matrix<double> b(4,1);
	double bd[] = {0.5,4.0/3,3.0/2,1.0};
	b.setdata(bd, 4);
	b.mprint();
	
	double& var = A(1,1);
	var = 6;
	cout << "** " << A(1,1) << endl;
	
	Matrix<double> x = gausselim(A,b);
	
	x.mprint();*/
	
	/*vector<double>* v = new vector<double>[n];
	for(int i = 0; i < n; i++)
		v[i].resize(m-i);
	
	qr(A, v);
	
	cout << "QR done.\n";
	A.mprint();

	delete [] v;*/
	cout << endl;
	return 0;
}
