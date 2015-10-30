#include <amesh3d.hpp>

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "diffmesh.control";
	ifstream conf(confile);
	string mesh1, mesh2, dum;

	conf >> dum; conf >> mesh1;
	conf >> dum; conf >> mesh2;
	conf.close();
	
	UMesh m1, m2;
	m1.readGmsh2(mesh1,3);
	m2.readGmsh2(mesh2,3);
	
	if(m1.gnpoin() != m2.gnpoin()) 
	{
		cout << "No. of points is not same!\n";
		return 0;
	}
	
	Matrix<double> mag1(m1.gnpoin(), 1);
	Matrix<double> mag2(m1.gnpoin(), 1);
	
	for(int i = 0; i < m1.gnpoin(); i++)
	{
		mag1(i) = 0;
		mag2(i) = 0;
		for(int j = 0; j < m1.gndim(); j++)
		{
			mag1(i) += m1.gcoords(i,j)*m1.gcoords(i,j);
			mag2(i) += m2.gcoords(i,j)*m2.gcoords(i,j);
		}
		mag1(i) = sqrt(mag1(i));
		mag2(i) = sqrt(mag2(i));
	}
	
	double tot = 0;
	for(int i = 0; i < m1.gnpoin(); i++)
	{
		tot += (mag1(i) - mag2(i))*(mag1(i) - mag2(i));
	}
	tot = sqrt(tot);
	
	cout << "\nl2 norm of difference is " << tot << endl;
	
	cout << endl;
	return 0;
}
	