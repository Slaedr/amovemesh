/*
This file defines a class Matrix to store a matrix in column-major or row-major storage.
It provides functions to replace rows or columns with vectors.
It provides a function to isolate sub-matrices given a starting i,j and offsets to both.

Aditya Kashi
April 17, 2015

Notes:
If A is a column-major matrix, A[i][j] == A[j * nrows + i] where i is the row-index and j is the column index.

--> TODO: Make 'storage' a template parameter rather than class member

TODO: Make a Vector class, maybe as a sub-class of Matrix
*/

// Comment out the following line for removing array index-range checks and thus improving performance
#define DEBUG 1

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_IOMANIP
#include <iomanip>
#endif

#ifndef _GLIBCXX_ARRAY
#include <array>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#define __AMATRIX2_H

//using namespace std;

const int WIDTH = 10;		// width of field for printing matrices

namespace amat {

double dabs(double x)
{
	if(x < 0) return (-1.0)*x;
	else return x;
}

enum MStype {ROWMAJOR, COLMAJOR};

template <class T, MSType storage>
class Matrix
{
private:
	int nrows;
	int ncols;
	int size;
	T* elems;

public:
	//No-arg constructor. Note: no memory allocation! Make sure Matrix::setup(int,int,MStype) is used.
	Matrix()
	{
		nrows = 1; ncols = 1;
		elems = new T;
	}

	// Constructor which sets storage to row-major
	Matrix(int nr, int nc)
	{
		if(nc==0)
		{
			std::cout << "\nError: Number of columns is zero. Setting it to 1.";
			nc=1;
		}
		if(nr==0)
		{
			std::cout << "\nError: Number of rows is zero. Setting it to 1.";
			nr=1;
		}
		nrows = nr; ncols = nc; size = nrows*ncols; storage = ROWMAJOR;
		elems = new T[nrows*ncols];
	}

	// Full-arg constructor
	Matrix(int nr, int nc)
	{
		if(nc==0)
		{
			std::cout << "\nError: Number of columns is zero. Setting it to 1.";
			nc=1;
		}
		if(nr==0)
		{
			std::cout << "\nError: Number of rows is zero. Setting it to 1.";
			nr=1;
		}
		nrows = nr; ncols = nc;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
	}

	Matrix(const Matrix<T,storage>& other)
	{
		nrows = other.nrows;
		ncols = other.ncols;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
		for(int i = 0; i < nrows*ncols; i++)
		{
			elems[i] = other.elems[i];
		}
	}

	~Matrix()
	{
		delete [] elems;
	}

	Matrix<T>& operator=(Matrix<T> rhs)
	{
#ifdef DEBUG
		if(this==&rhs) return *this;		// check for self-assignment
#endif
		nrows = rhs.nrows;
		ncols = rhs.ncols;
		storage = rhs.storage;
		size = nrows*ncols;
		delete [] elems;
		elems = new T[nrows*ncols];
		for(int i = 0; i < nrows*ncols; i++)
		{
			elems[i] = rhs.elems[i];
		}
		return *this;
	}

	//Separate setup function in case no-arg constructor has to be used
	void setup(int nr, int nc, MStype st)
	{
		if(nc==0)
		{
			std::cout << "\nError: Number of columns is zero. Setting it to 1.";
			nc=1;
		}
		if(nr==0)
		{
			std::cout << "\nError: Number of rows is zero. Setting it to 1.";
			nr=1;
		}
		nrows = nr; ncols = nc; storage = st;
		size = nrows*ncols;
		delete [] elems;
		elems = new T[nrows*ncols];
	}

	// no deleting earlier allocation: use in case of Matrix<t>* (pointer to Matrix<t>)
	void setupraw(int nr, int nc, MStype st)
	{
		//std::cout << "\nEntered setupraw";
		if(nc==0)
		{
			std::cout << "\nError: Number of columns is zero. Setting it to 1.";
			nc=1;
		}
		if(nr==0)
		{
			std::cout << "\nError: Number of rows is zero. Setting it to 1.";
			nr=1;
		}
		nrows = nr; ncols = nc; storage = st;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
	}

	void zeros()
	{
		for(int i = 0; i < size; i++)
			elems[i] = (T)(0.0);
	}

	void ones()
	{
		for(int i = 0; i < size; i++)
			elems[i] = 1;
	}

	void identity()
	{
		T one = (T)(1);
		T zero = (T)(0);
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				if(i==j) operator()(i,j) = one;
				else operator()(i,j) = zero;
	}

	// function to set matrix elements from a ROW-MAJOR std::vector
	void setdata(const T* A, int sz)
	{
#ifdef DEBUG
		if(sz != size)
		{
			std::cout << "\nError in setdata: argument size does not match matrix size";
			return;
		}
#endif
		if(storage == COLMAJOR)
			for(int i = 0; i < nrows; i++)
				for(int j = 0; j < ncols; j++)
					elems[j*nrows+i] = A[i*ncols+j];
		else
			for(int i = 0; i < nrows; i++)
				for(int j = 0; j < ncols; j++)
					elems[i*ncols+j] = A[i*ncols+j];
	}

	T get(int i, int j=0)
	{
#ifdef DEBUG
		if(i>=nrows || j>=ncols) { std::cout << "Matrix: get(): Index beyond array size(s)\n"; return 0; }
#endif
		return elems[i*ncols + j];
	}

	void set(T data, int i, int j)
	{
#ifdef DEBUG
		if(i>=nrows || j>=ncols) { std::cout << "Matrix: set(): Index beyond array size(s)\n"; return; }
#endif
		//if(i < 0 || j < 0) {std::cout << "Matrix: get(): Negative index!\n"; return 0; }
		elems[i*ncols + j] = data;
	}

	int rows() { return nrows; }
	int cols() { return ncols; }

	void mprint()
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << elems[i*ncols+j];
			std::cout << std::endl;
		}
	}

	void fprint(std::ofstream& outfile)
	{
		//outfile << '\n';
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << elems[i*ncols+j];
			outfile << '\n';
		}
	}

	void fread(std::ifstream& infile)
	{
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				infile >> elems[i*ncols + j];
	}

	// For expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141
	T& operator()(int x, int y=0)
	{
#ifdef DEBUG
		if(x>=nrows || y>=ncols) { std::cout << "Matrix (): Index beyond array size(s)\n"; return elems[0]; }
#endif
		return elems[x*ncols + y];
	}

	T maxincol(int j)
	{
		T max = get(0,j);
		for(int i = 0; i < nrows; i++)
			if(max < get(i,j)) max = get(i,j);
		return max;
	}

	T maxinrow(int i)
	{
		T max = get(i,0);
		for(int j = 0; j < nrows; j++)
			if(max < get(i,j)) max = get(i,j);
		return max;
	}

	T max()
	{
		T max = elems[0];
		for(int i = 0; i < size; i++)
			if(elems[i] > max) max = elems[i];
		return max;
	}

	T absmax()
	{
		T max = abs(elems[0]);
		for(int i = 0; i < size; i++)
			if(abs(elems[i]) > max) max = abs(elems[i]);
		return max;
	}

	double dabsmax()
	{
		double max = dabs((double)elems[0]);
		for(int i = 0; i < size; i++)
			if(dabs(elems[i]) > max) max = dabs(elems[i]);
		return max;
	}

	T minincol(int j)
	{
		T min = get(0,j);
		for(int i = 0; i < nrows; i++)
			if(min > get(i,j)) min = get(i,j);
		return min;
	}

	T mininrow(int i)
	{
		T max = get(i,0);
		for(int j = 0; j < nrows; j++)
			if(max > get(i,j)) max = get(i,j);
		return max;
	}

	T min()
	{
		T max = elems[0];
		for(int i = 0; i < size; i++)
			if(elems[i] < max) max = elems[i];
		return max;
	}

	// function to return a sub-matrix of this matrix
	Matrix<T,storage> sub(int startr, int startc, int offr, int offc)
	{
		Matrix<T,storage> B(offr, offc);
		for(int i = 0; i < offr; i++)
			for(int j = 0; j < offc; j++)
				B(i,j) = elems[(startr+i)*ncols + startc + j];
		return B;
	}

	//Function that returns a given column of the matrix as a row-major matrix
	Matrix<T,storage> col(int j)
	{
		Matrix<T,storage> b(nrows, 1);
		for(int i = 0; i < nrows; i++)
			b(i,0) = elems[i*ncols + j];
		return b;
	}

	Matrix<T,storage> row(int i)
	{
		Matrix<T,storage> b(1, ncols);
		for(int j = 0; j < ncols; j++)
			b(0,j) = elems[i*ncols + j];
		return b;
	}

	/*//Function to return a reference to a given column of the matrix
	Matrix<T>& colr(int j)
	{
		//Matrix<T>* b(nrows, 1);
		Matrix<T>* b; b->elems.reserve(nrows);
		if(storage==COLMAJOR)
			for(int i = 0; i < nrows; i++)
				b->elems[i] = elems[j*nrows + i];
		else
			for(int i = 0; i < nrows; i++)
				b.elems[i] = &elems[i*ncols + j];
		return *b;
	} */

	// Function for replacing a column of the matrix with a vector. NOTE: No check for whether b is really a vector - which it must be.
	void replacecol(int j, Matrix<T> &b)
	{
#ifdef DEBUG
		if(b.cols() != 1 || b.rows() != nrows) { std::cout << "\nSize error in replacecol"; return; }
#endif
		for(int i = 0; i < nrows; i++)
			elems[i*ncols + j] = b.elems[i];
	}

	//Function for replacing a row
	void replacerow(int i, Matrix<T> b)
	{
#ifdef DEBUG
		if(b.cols() != ncols || b.rows() != 1) { std::cout << "\nSize error in replacerow"; return; }
#endif
		for(int j = 0; j < ncols; j++)
			elems[i*ncols + j] = b.elems[j];
	}

	//transpose
	Matrix<T,storage> trans()
	{
		Matrix<T,storage> t(ncols, nrows, storage);
		for(int i = 0; i < ncols; i++)
			for(int j = 0; j < nrows; j++)
				t(i,j) = get(j,i);
		return t;
	}

	// Multiply a matrix by a scalar. Note: only expressions of type A*3 work, not 3*A
	Matrix<T,storage> operator*(T num)
	{
		Matrix<T,storage> A(nrows,ncols,storage);
		for(int i = 0; i < size; i++)
			A.elems[i] = elems[i] * num;
		return A;
	}

	Matrix<T> operator+(Matrix<T> B)
	{
#ifdef DEBUG
		if(nrows != B.rows() || ncols != B.cols())
		{
			std::cout << "! Matrix: Addition cannot be performed due to incompatible sizes\n";
			Matrix<T> C(1,1);
			return C;
		}
#endif
		Matrix<T> C(nrows, ncols);

		for(int i = 0; i < C.size; i++)
			C.elems[i] = elems[i] + B.elems[i];
		return C;
	}

	Matrix<T> operator-(Matrix<T> B)
	{
#ifdef DEBUG
		if(nrows != B.rows() || ncols != B.cols())
		{
			std::cout << "! Matrix: Subtraction cannot be performed due to incompatible sizes\n";
			Matrix<T> C(1,1);
			return C;
		}
#endif
		Matrix<T> C(nrows, ncols);

		for(int i = 0; i < C.size; i++)
			C.elems[i] = elems[i] - B.elems[i];
		return C;
	}

	Matrix<T> operator*(Matrix<T> B)
	{
		Matrix<double> C(nrows, B.cols(), storage);
		C.zeros();
#ifdef DEBUG
		if(ncols != B.rows())
		{
			std::cout << "! Matrix: Multiplication cannot be performed - incompatible sizes!\n";
			return C;
		}
#endif
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < B.cols(); j++)
				for(int k = 0; k < ncols; k++)
					C(i,j) += get(i,k) * B(k,j);

		return C;
	}
};


//Specialization for COLAMJOR
template <class T>
T Matrix<T, COLMAJOR>::get(int i, int j=0)
{
#ifdef DEBUG
	if(i>=nrows || j>=ncols) { std::cout << "Matrix: get(): Index beyond array size(s)\n"; return 0; }
#endif
	return elems[j*nrows + i];
}

void Matrix<T, COLAMJOR>::set(T data, int i, int j)
{
#ifdef DEBUG
	if(i>=nrows || j>=ncols) { std::cout << "Matrix: set(): Index beyond array size(s)\n"; return; }
#endif
	//if(i < 0 || j < 0) {std::cout << "Matrix: get(): Negative index!\n"; return 0; }
	elems[j*nrows + i] = data;
}

template <class T>
void Matrix<T, COLAMJOR>::mprint()
{
	std::cout << "\n";
	for(int i = 0; i < nrows; i++)
	{
		for(int j = 0; j < ncols; j++)
			std::cout << std::setw(WIDTH) << elems[j*nrows+i];
		std::cout << std::endl;
	}
}

template <class T>
void Matrix<T, COLAMJOR>::fprint(std::ofstream& outfile)
{
	//outfile << '\n';
	for(int i = 0; i < nrows; i++)
	{
		for(int j = 0; j < ncols; j++)
			outfile << " " << elems[j*nrows+i];
		outfile << '\n';
	}
}

template <class T>
void Matrix<T, COLAMJOR>::fread(std::ifstream& infile)
{
	for(int i = 0; i < nrows; i++)
		for(int j = 0; j < ncols; j++)
			infile >> elems[j*nrows + i];
}

// For expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141
template <class T>
T& Matrix<T, COLAMJOR>::operator()(int x, int y=0)
{
#ifdef DEBUG
	if(x>=nrows || y>=ncols) { std::cout << "Matrix (): Index beyond array size(s)\n"; return elems[0]; }
#endif
	return elems[y*nrows + x];
}

template <class T>
Matrix<T,COLMAJOR> Matrix<T, COLAMJOR>::sub(int startr, int startc, int offr, int offc)
{
	Matrix<T,COLMAJOR> B(offr, offc);
	for(int i = 0; i < offr; i++)
		for(int j = 0; j < offc; j++)
			B(i,j) = elems[(startc+j)*nrows + startr + i];
	return B;
}

//Function that returns a given column of the matrix as a row-major matrix
template <class T>
Matrix<T,COLMAJOR> Matrix<T, COLAMJOR>::col(int j)
{
	Matrix<T,COLMAJOR> b(nrows, 1);
	for(int i = 0; i < nrows; i++)
		b(i,0) = elems[j*nrows + i];
	return b;
}

template <class T>
Matrix<T,COLMAJOR> Matrix<T, COLMAJOR>::row(int i)
{
	Matrix<T,COLMAJOR> b(1, ncols);
	for(int j = 0; j < ncols; j++)
		b(0,j) = elems[j*nrows + i];
	return b;
}

// Function for replacing a column of the matrix with a vector. NOTE: No check for whether b is really a vector - which it must be.
template <class T>
void Matrix<T, COLAMJOR>::replacecol(int j, Matrix<T> &b)
{
#ifdef DEBUG
	if(b.cols() != 1 || b.rows() != nrows) { std::cout << "\nSize error in replacecol"; return; }
#endif
	for(int i = 0; i < nrows; i++)
		elems[j*nrows + i] = b.elems[i];
}

//Function for replacing a row
template <class T>
void Matrix<T, COLAMJOR>::replacerow(int i, Matrix<T> b)
{
#ifdef DEBUG
	if(b.cols() != ncols || b.rows() != 1) { std::cout << "\nSize error in replacerow"; return; }
#endif
	for(int j = 0; j < ncols; j++)
		elems[j*nrows + i] = b.elems[j];
}


} //end namespace amat
