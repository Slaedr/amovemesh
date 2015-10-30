#ifndef _GLIBCXX_STRING
#include <string>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef __AMATRIX2_H
#include <amatrix2.h>
#endif

#define __AINPUT_H

using namespace amat;
using namespace std;

namespace acfd {

class LoadInitData
{
	string file;
	Matrix<double> data;
public:
	LoadInitData(string fname)
	{
		file = fname;
		ifstream ifile(fname);
		
		//TODO: load initial data into data
		
		ifile.close();
	}
	
	Matrix<double> loadinitdata()
	{ return data; }
};

class LoadBoundaryData
{
	string file;
};

} // end namespace acfd