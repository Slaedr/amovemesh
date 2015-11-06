#include <amesh2dh.hpp>

using namespace amat;
using namespace acfd;
using namespace std;

int main()
{
	string confilename = "convertformat.control";
	ifstream conf(confilename);
	string dum, inmesh, informat, outmesh, outformat;
	conf >> dum; conf >> inmesh;
	conf >> dum; conf >> informat;
	conf >> dum; conf >> outmesh;
	conf >> dum; conf >> outformat;
	conf.close();

	cout << "Input file is of type " << informat << ". Writing as " << outformat << ".\n";

	UMesh2dh m;
	if(informat == "domn")
		m.readDomn(inmesh);
	else if(informat == "msh")
		m.readGmsh2(inmesh,2);
	else {
		cout << "Invalid format. Exiting." << endl;
		return -1;
	}

	m.writeGmsh2(outmesh);

	cout << endl;
	return 0;
}
