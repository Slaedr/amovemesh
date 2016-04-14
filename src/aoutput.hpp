/** @file aoutput.hpp
 * @brief A collection of subroutines to write mesh data to various kinds of output formats
 */

#ifndef __AOUTPUT_H

#ifndef __AMATRIX2_H
#include <amatrix2.hpp>
#endif

#ifndef __AMESH2DGENERAL_H
#include <amesh2d.hpp>
#endif

#ifndef __AMESH2DHYBRID_H
#include <amesh2dh.hpp>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLIBCXX_STRING
#include <string>
#endif

#define __AOUTPUT_H 1

using namespace std;
using namespace amat;
using namespace amc;

void writeScalarsVectorToVtu_CellData(string fname, const UMesh2d& m, const Matrix<double>& x, string scaname[], const Matrix<double>& y, string vecname);

// -------------- Implemetations -----------------------

/** Writes multiple scalar data sets and one vector data set, all cell-centered data, to a file in VTU format.
 * If either x or y is a 0x0 matrix, it is ignored.
 * \param fname is the output vtu file name
 */
void writeScalarsVectorToVtu_CellData(string fname, const UMesh2d& m, const Matrix<double>& x, string scaname[], const Matrix<double>& y, string vecname)
{
	int elemcode = 5;
	if(m.gnnode() == 4)
		elemcode = 9;
	else if(m.gnnode() == 6)
		elemcode = 22;
	else if(m.gnnode() == 8)
		elemcode = 23;
	else if(m.gnnode() == 9)
		elemcode = 28;
	
	cout << "aoutput: Writing vtu output...\n";
	ofstream out(fname);

	int nscalars = x.cols();

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//out << "\t\t<CellData Scalars=\""<<scaname[0]<< "\" Vectors=\"" << vecname << "\">\n";
	
	if(x.msize()>0 || y.msize()>0) {
		out << "\t\t<CellData ";
		if(x.msize() > 0)
			out << "Scalars=\"" << scaname[0] << "\" ";
		if(y.msize() > 0)
			out << "Vectors=\"" << vecname << "\"";
		out << ">\n";
	}
	
	//enter cell scalar data if available
	if(x.msize() > 0) {
		//cout << "aoutput: Writing scalars..\n";
		for(int in = 0; in < nscalars; in++)
		{
			out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << scaname[in] << "\" Format=\"ascii\">\n";
			for(int i = 0; i < m.gnelem(); i++)
				out << "\t\t\t\t" << x.get(i,in) << '\n';
			out << "\t\t\t</DataArray>\n";
		}
		//cout << "aoutput: Scalars written.\n";
	}

	//enter vector cell data if available
	if(y.msize() > 0) {
		out << "\t\t\t<DataArray type=\"Float64\" Name=\"" << vecname << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		for(int i = 0; i < m.gnelem(); i++)
		{
			out << "\t\t\t\t";
			for(int idim = 0; idim < y.cols(); idim++)
				out << y.get(i,idim) << " ";
			if(y.cols() == 2)
				out << "0.0 ";
			cout << '\n';
		}
		out << "\t\t\t</DataArray>\n";
	}
	if(x.msize() > 0 || y.msize() > 0)
		out << "\t\t</CellData>\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
	{
		out << "\t\t\t";
		for(int idim = 0; idim < m.gndim(); idim++)
			out << m.gcoords(i,idim) << " ";
		if(m.gndim() == 2)
			out << "0.0 ";
		out << '\n';
	}
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++) {
		out << "\t\t\t\t"; 
		for(int inode = 0; inode < m.gnnode(); inode++)	
			out << m.ginpoel(i,inode) << " ";
		out << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.gnnode()*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << elemcode << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}


/// Writes a quadratic mesh in VTU format.
/** VTK does not have a 9-node quadrilateral, so we ignore the cell-centered note for output.
 */
void writeQuadraticMeshToVtu(string fname, UMesh2d& m)
{
	int elemcode = 5;
	if(m.gnnode() == 4)
		elemcode = 9;
	else if(m.gnnode() == 6)
		elemcode = 22;
	else if(m.gnnode() == 8)
		elemcode = 23;
	else if(m.gnnode() == 9)
		elemcode = 23;

	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
	{
		out << "\t\t\t\t";
		for(int j = 0; j < m.gnnode(); j++)
			out << m.ginpoel(i,j) << " ";
		out << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.gnnode()*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << elemcode << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

/// Writes a hybrid mesh in VTU format.
/** VTK does not have a 9-node quadrilateral, so we ignore the cell-centered note for output.
 */
void writeMeshToVtu(string fname, UMesh2dh& m)
{
	

	cout << "Writing vtu output...\n";
	ofstream out(fname);

	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "<UnstructuredGrid>\n";
	out << "\t<Piece NumberOfPoints=\"" << m.gnpoin() << "\" NumberOfCells=\"" << m.gnelem() << "\">\n";

	//enter points
	out << "\t\t<Points>\n";
	out << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnpoin(); i++)
		out << "\t\t\t" << m.gcoords(i,0) << " " << m.gcoords(i,1) << " " << 0.0 << '\n';
	out << "\t\t</DataArray>\n";
	out << "\t\t</Points>\n";

	//enter cells
	out << "\t\t<Cells>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
	{
		out << "\t\t\t\t";
		for(int j = 0; j < m.gnnode(i); j++)
			out << m.ginpoel(i,j) << " ";
		out << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
		out << "\t\t\t\t" << m.gnnode(i)*(i+1) << '\n';
	out << "\t\t\t</DataArray>\n";
	out << "\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
	for(int i = 0; i < m.gnelem(); i++)
	{
		int elemcode = 5;
		if(m.gnnode(i) == 4)
			elemcode = 9;
		else if(m.gnnode(i) == 6)
			elemcode = 22;
		else if(m.gnnode(i) == 8)
			elemcode = 23;
		else if(m.gnnode(i) == 9)
			elemcode = 23;
			out << "\t\t\t\t" << elemcode << '\n';
	}
	out << "\t\t\t</DataArray>\n";
	out << "\t\t</Cells>\n";

	//finish upper
	out << "\t</Piece>\n";
	out << "</UnstructuredGrid>\n";
	out << "</VTKFile>";
	out.close();
	cout << "Vtu file written.\n";
}

#endif