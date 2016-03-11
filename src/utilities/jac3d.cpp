#include <amesh3d.hpp>

using namespace amc;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Please give a file name\n";
		return -1;
	}
	
	ifstream fin(argv[1]);
	string dum, infile, outfile; double thresh;
	fin >> dum; fin >> infile;
	fin >> dum; fin >> outfile;
	fin >> dum; fin >> thresh;
	fin.close();
	
	UMesh m;
	m.readGmsh2(infile,3);
	m.compute_jacobians();
	
	ofstream fout(outfile);
	amc_int njac = 0;
	for(int i = 0; i < m.gnelem(); i++)
	{
		fout << m.gjacobians(i) << '\n';
		if(m.gjacobians(i) < thresh)
			njac++;
	}
	
	fout.close();
	cout << "Number of elements with jacobian less than " << thresh << " is " << njac << endl;
	return 0;
}
	
	