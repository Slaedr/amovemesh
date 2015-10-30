#include "abowyerwatson.hpp"

using namespace std;
using namespace amat;
using namespace acfd;

int main()
{
	string confile = "triangulation.control";
	string inp;
	string outp;
	string dum;
	ifstream conf(confile);
	conf >> dum; conf >> inp;
	conf >> dum; conf >> outp;
	conf.close();

	ifstream infile(inp);
	Matrix<double> points;
	UTriMesh m(infile);
	infile.close();
	points = m.return_boundary_points();

	/*ofstream outfile("coarse_points.dat");
	points.fprint(outfile);
	outfile.close();*/

	/*ifstream infile(inp);
	Matrix<double> points;
	points.fread(infile);
	infile.close();*/

	Delaunay2D prob(&points, points.rows());
	prob.bowyer_watson();
	prob.writeGmsh2(outp);

	cout << endl;
	return 0;
}
