#include <abowyerwatson.hpp>
#include <adistbinsort.hpp>

using namespace std;
using namespace amat;
using namespace amc;

int main(int argc, char* argv[])
{
	string confile(argv[1]);
	string inp;
	string outp;
	string dum;
	ifstream conf(confile);
	conf >> dum; conf >> inp;
	conf >> dum; conf >> outp;
	conf.close();

	ifstream infile(inp);
	Matrix<double>* points;
	UMesh2dh m(infile);
	infile.close();
	points = m.getcoords();

	/*ofstream outfile("coarse_points.dat");
	points.fprint(outfile);
	outfile.close();*/

	/*ifstream infile(inp);
	Matrix<double> points;
	points.fread(infile);
	infile.close();*/

	Delaunay2D prob(points);
	prob.bowyer_watson();
	prob.writeGmsh2(outp);

	cout << endl;
	return 0;
}
