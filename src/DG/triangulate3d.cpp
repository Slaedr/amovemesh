/** @file triangulate3d.cpp
 * @brief Construct Delaunay tetrahedralization of a set of 3D points.
 * @date May 25, 2016
 */

#include <adistbinsort.hpp>
#include <abowyerwatson3d.hpp>

using namespace std;
using namespace amat;
using namespace amc;

void read_points(string finname, Matrix<amc_real>& points, amc_int& npoints)
{
	string dum; int indum; char chdum;
	
	ifstream fin(finname);
	getline(fin,dum);
	fin >> npoints;
	fin >> chdum;
	getline(fin,dum);

	points.setup(npoints,3);

	for(int i = 0; i < npoints; i++)
	{
		fin >> indum;
		for(int j = 0; j < 3; j++)
			fin >> points(i,j);
	}

	fin.close();
}

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Insufficient arguments. Please provide a control file name.\n";
		return -1;
	}
	
	string inpfile, outfile, dum;
	int numbins[3];
	ifstream fin(argv[1]);
	fin >> dum; fin >> inpfile;
	fin >> dum; fin >> outfile;
	fin >> dum; fin >> numbins[0];
	fin >> dum; fin >> numbins[1];
	fin >> dum; fin >> numbins[2];
	fin.close();
	
	amc_int npoints;
	Matrix<amc_real> points, *sortedpoints;

	read_points(inpfile, points, npoints);
	
	DistBinSort<3>* sorter;
	if(numbins[0] == 0 || numbins[1] == 0 || numbins[2] == 0)
	{
		sortedpoints = &points;
		cout << "No bin sort" << endl;
	}
	else
	{
		sortedpoints = new Matrix<amc_real>(npoints,3);
		cout << "Bin sort with x-, y- and z-bins " << numbins[0] << " " << numbins[1] << " " << numbins[2] << endl;
		sorter = new DistBinSort<3>(&points, numbins, sortedpoints);
	}
	
	Delaunay3d d3(sortedpoints);
	d3.bowyer_watson();

	//d3.compute_jacobians();
	//bool negjac = d3.detect_negative_jacobians();
	d3.check();

	d3.writeGmsh2(outfile);
	cout << "Output written to " << outfile << endl;
	
	if(numbins[0] != 0 && numbins[1] != 0 && numbins[2] != 0)
	{
		delete sortedpoints;
		delete sorter;
	}
	
	cout << endl;
	return 0;
}
