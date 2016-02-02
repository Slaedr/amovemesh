#include <amesh2d.hpp>

using namespace amat;
using namespace acfd;
using namespace std;

int main()
{
	string inmesh = "../input/naca0012-coarse.domn";
	string outfile = "../input/naca0012-coarse.points";

	UMesh2d m;
	m.readDomn(inmesh);

	Matrix<int> low(m.gnpoin(),1);
	low.zeros();
	for(int i = 0; i < m.gnelem(); i++)
		for(int j = 0; j < 8; j++)
			low(m.ginpoel(i,j)) = 1;

	int np = 0, k = 1;
	for(int i = 0; i < m.gnpoin(); i++)
		if(low.get(i) == 1)
			np++;
	
	ofstream fout(outfile);
	fout << "npoin\n" << np << '\n';
	for(int i = 0; i < m.gnpoin(); i++)
		if(low.get(i) == 1)
		{
			fout << k << " ";
			for(int idim = 0; idim < m.gndim(); idim++)
				fout << " " << m.gcoords(i,idim);
			fout << '\n';
			k++;
		}

	fout.close();

	cout << endl;
	return 0;
}
