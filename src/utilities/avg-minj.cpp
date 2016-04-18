/** @file avg-minj.cpp
 * @brief Computes average value of cell-centered scalar data from a Gmsh pos (msh) file
 * @author Aditya Kashi
 * @date April 18, 2016
 * 
 * The input format is
 * 
 * 		number_of_elements
 * 		elem_number        value
 * 		elem_number        value
 * 		.
 * 		.
 */

#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Please provide a file from which to read values." << endl;
		return -1;
	}
	
	string mfile(argv[1]);
	int nelem, dum;
	double jac, avg = 0;
	
	ifstream fin(mfile);
	fin >> nelem;
	
	for(int i = 0; i < nelem; i++)
	{
		fin >> dum;
		fin >> jac;
		
		if(jac <= 0.0) cout << "Negative jacobian for element " << i << "!" << endl;
		
		avg += jac;
	}
	fin.close();
	
	avg /= static_cast<double>(nelem);
	
	cout << "The average of minimum Jacobians over all cells is " << avg << endl;
	
	cout << endl;
	return 0;
}