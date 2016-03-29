/** @brief Converts a mesh file from Dr Luo's Domn format (containing no bface data) to Gmsh 2 format
 */
#include <amesh3d.hpp>

using namespace amat;
using namespace amc;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		cout << "Please give an input file and an output file!";
		return -1;
	}
	string inmeshname(argv[1]);
	string outmeshname(argv[2]);

	UMesh m;
	m.readDomn(inmeshname);
	m.writeGmsh2(outmeshname);

	cout << endl;
	return 0;
}
