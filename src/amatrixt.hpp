/**
 * @file amatrixt.hpp
 * @brief Defines a classes to manipulate matrices.
 *
 * The type of scalars to be stored and the storage order are template parameters.
 * Since the storage order is also a template parameter here, as opposed to amatrix2.hpp, there is no need to check storage order
 * for get/set operations during runtime. Hence, this implementation is faster.
 * @author Aditya Kashi
 * @date April 3, 2016
 */

/**
 * \namespace amat
 * \brief Includes all array and matrix storage classes, as well as linear algebra.
 */

/**
 * \class Matrix
 * \brief Stores a 2D array or dense matrix.
 * 
 * Notes:
 * If A is a column-major matrix, A[i][j] == A[j * nrows + i] where i is the row-index and j is the column index.
 */

#ifndef __AMATRIXT_H

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

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef __ACONSTANTS_H
#include <aconstants.h>
#endif

#define __AMATRIXT_H

#ifndef MATRIX_DOUBLE_PRECISION
#define MATRIX_DOUBLE_PRECISION 14
#endif

namespace amat {

/// Real type
using amc::amc_real;

// Integer type
using amc::amc_int;

const int WIDTH = 10;		// width of field for printing matrices

inline double dabs(double x)
{
	if(x < 0) return (-1.0)*x;
	else return x;
}
inline double minmod(double a, double b)
{
	if(a*b>0 && dabs(a) <= dabs(b)) return a;
	else if (a*b>0 && dabs(b) < dabs(a)) return b;
	else return 0.0;
}

/// Matrix storage type. Note that its usage is deprecated.
enum MStype {ROWMAJOR, COLMAJOR};

template <class T, MStype storage>
class Matrix
{
private:
	int nrows;
	int ncols;
	int size;
	T* elems;
	bool isalloc;

public:
	///No-arg constructor. Note: no memory allocation! Make sure Matrix::setup(int,int,MStype) is used.
	Matrix()
	{
		nrows = 0; ncols = 0; size = 0;
		isalloc = false;
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
		isalloc = true;
	}

	~Matrix()
	{
		if(isalloc == true)	
			delete [] elems;
		isalloc = false;
	}
};

/// Specialization of the Matrix class to row-major storage
template <class T> 
class Matrix<T, ROWMAJOR>
{
private:
	int nrows;
	int ncols;
	int size;
	T* elems;
	bool isalloc;

public:
	///No-arg constructor. Note: no memory allocation! Make sure Matrix::setup(int,int,MStype) is used.
	Matrix()
	{
		nrows = 0; ncols = 0; size = 0;
		isalloc = false;
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
		isalloc = true;
	}

	/// Initialize the array from another row-major array
	Matrix(const Matrix<T,ROWMAJOR>& other)
	{
		nrows = other.nrows;
		ncols = other.ncols;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows*ncols; i++)
		{
			elems[i] = other.elems[i];
		}
	}
	
	/// Initialize the array from a column-major array
	Matrix(const Matrix<T,COLMAJOR>& other)
	{
		nrows = other.nrows;
		ncols = other.ncols;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				elems[i*ncols+j] = other.get(i,j);
	}

	~Matrix()
	{
		if(isalloc == true)	
			delete [] elems;
		isalloc = false;
	}

	Matrix<T,ROWMAJOR>& operator=(const Matrix<T,ROWMAJOR>& rhs)
	{
#ifdef DEBUG
		if(this==&rhs) return *this;		// check for self-assignment
#endif
		nrows = rhs.nrows;
		ncols = rhs.ncols;
		size = nrows*ncols;
		if(isalloc == true)
			delete [] elems;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows*ncols; i++)
		{
			elems[i] = rhs.elems[i];
		}
		return *this;
	}

	Matrix<T,ROWMAJOR>& operator=(const Matrix<T,COLMAJOR>& rhs)
	{
#ifdef DEBUG
		if(this==&rhs) return *this;		// check for self-assignment
#endif
		nrows = rhs.nrows;
		ncols = rhs.ncols;
		size = nrows*ncols;
		if(isalloc == true)
			delete [] elems;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				elems[i*ncols+j] = rhs.get(i,j);
		return *this;
	}
	/// Separate setup function in case no-arg constructor has to be used
	void setup(int nr, int nc)
	{
		if(nc==0)
		{
			std::cout << "Matrix: setup(): Error: Number of columns is zero. Setting it to 1.\n";
			nc=1;
		}
		if(nr==0)
		{
			std::cout << "Matrix(): setup(): Error: Number of rows is zero. Setting it to 1.\n";
			nr=1;
		}
		nrows = nr; ncols = nc;
		size = nrows*ncols;
		if(isalloc == true)
			delete [] elems;
		elems = new T[nrows*ncols];
		isalloc = true;
	}

	/// Setup without deleting earlier allocation: use in case of Matrix<t>* (pointer to Matrix<t>)
	void setupraw(int nr, int nc)
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
		isalloc = true;
	}
	
	/// Fill the matrix with zeros.
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

	/// function to set matrix elements from a ROW-MAJOR array
	void setdata_rowmajor(const T* A, int sz)
	{
#ifdef DEBUG
		if(sz != size)
		{
			std::cout << "\nError in setdata: argument size does not match matrix size";
			return;
		}
#endif
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				elems[i*ncols+j] = A[i*ncols+j];
	}

	T get(const int i, const int j=0) const
	{
#ifdef DEBUG
		if(i>=nrows || j>=ncols) { std::cout << "Matrix: get(): Index beyond array size(s)\n"; return 0; }
		if(i < 0 || j < 0) { std::cout << "Matrix: get(): Index less than 0!\n"; return 0; }
#endif
		return elems[i*ncols + j];
	}

	void set(int i, int j, T data)
	{
#ifdef DEBUG
		if(i>=nrows || j>=ncols) { std::cout << "Matrix: set(): Index beyond array size(s)\n"; return; }
		if(i < 0 || j < 0) {std::cout << "Matrix: get(): Negative index!\n"; return 0; }
#endif
		elems[i*ncols + j] = data;
	}

	int rows() const { return nrows; }
	int cols() const { return ncols; }
	int msize() const { return size; }

	/// Prints the matrix to standard output.
	void mprint() const
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << std::setprecision(WIDTH/2+1) << elems[i*ncols+j];
			std::cout << std::endl;
		}
	}

	/// Prints the matrix to file
	void fprint(std::ofstream& outfile) const
	{
		outfile << std::setprecision(MATRIX_DOUBLE_PRECISION);
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << elems[i*ncols+j];
			outfile << '\n';
		}
	}

	/// Reads matrix from file
	void fread(std::ifstream& infile)
	{
		infile >> nrows; infile >> ncols;
		size = nrows*ncols;
		delete [] elems;
		elems = new T[nrows*ncols];
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				infile >> elems[i*ncols + j];
	}

	/// Getter/setter function for expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141
	T& operator()(const int x, const int y=0)
	{
#ifdef DEBUG
		if(x>=nrows || y>=ncols) { std::cout << "Matrix (): Index beyond array size(s)\n"; return elems[0]; }
		if(x < 0 || y < 0) {std::cout << "Matrix: (): Negative index!\n"; return 0; }
#endif
		return elems[x*ncols + y];
	}
	
	/// Const getter/setter function for expressions like x = A(1,2)
	const T& operator()(const int x, const int y=0) const
	{
#ifdef DEBUG
		if(x>=nrows || y>=ncols) { std::cout << "Matrix (): Index beyond array size(s)\n"; return elems[0]; }
		if(x < 0 || y < 0) {std::cout << "Matrix: (): Negative index!\n"; return 0; }
#endif
		return elems[x*ncols + y];
	}

	T maxincol(int j) const
	{
		T max = get(0,j);
		for(int i = 0; i < nrows; i++)
			if(max < get(i,j)) max = get(i,j);
		return max;
	}

	T maxinrow(int i) const
	{
		T max = get(i,0);
		for(int j = 0; j < nrows; j++)
			if(max < get(i,j)) max = get(i,j);
		return max;
	}

	T max() const
	{
		T max = elems[0];
		for(int i = 0; i < size; i++)
			if(elems[i] > max) max = elems[i];
		return max;
	}

	T absmax() const
	{
		T max = abs(elems[0]);
		for(int i = 0; i < size; i++)
			if(abs(elems[i]) > max) max = abs(elems[i]);
		return max;
	}

	/// Returns the magnitude of the element with largest magnitude
	double dabsmax() const
	{
		double max = dabs((double)elems[0]);
		for(int i = 0; i < size; i++)
			if(dabs(elems[i]) > max) max = dabs(elems[i]);
		return max;
	}

	T minincol(int j) const
	{
		T min = get(0,j);
		for(int i = 0; i < nrows; i++)
			if(min > get(i,j)) min = get(i,j);
		return min;
	}

	T mininrow(int i) const
	{
		T max = get(i,0);
		for(int j = 0; j < nrows; j++)
			if(max > get(i,j)) max = get(i,j);
		return max;
	}

	T min() const
	{
		T max = elems[0];
		for(int i = 0; i < size; i++)
			if(elems[i] < max) max = elems[i];
		return max;
	}

	T average() const
	{
		T avg = 0;
		for(int i = 0; i < size; i++)
			avg += elems[i];
		avg = avg/size;
		return avg;
	}

	// sums the square of all elements in the matrix and returns the square root of this sum
	T l2norm() const		
	{
		T tot = 0;
		for(int i = 0; i < size; i++)
		{
			tot += elems[i]*elems[i];
		}
		tot = std::sqrt(tot);
		return tot;
	}

	/// function to return a sub-matrix of this matrix
	template <MStype storage> Matrix<T,storage> sub(int startr, int startc, int offr, int offc) const
	{
		Matrix<T,storage> B(offr, offc);
		for(int i = 0; i < offr; i++)
			for(int j = 0; j < offc; j++)
				B(i,j) = get(startr+i,startc+j);
		return B;
	}

	/// Function that returns a given column of the matrix as a column (mx1) matrix
	template <MStype storage> Matrix<T,storage> col(int j) const
	{
		Matrix<T,storage> b(nrows, 1);
		for(int i = 0; i < nrows; i++)
			b(i,0) = get(i,j);
		return b;
	}

	template <MStype storage> Matrix<T,storage> row(int i) const
	{
		Matrix<T,storage> b(1, ncols);
		for(int j = 0; j < ncols; j++)
			b(0,j) = get(i,j);
		return b;
	}

	//Function to return a pointer to a given row of the matrix
	T* row_p(int i)
	{
		T* b = &elems[i*ncols];
		return b;
	}

	/// Function for replacing a column of the matrix with a vector. NOTE: No check for whether b is really a vector - which it must be.
	/** Inefficient, since this is a row-major array.
	 */
	template <MStype storage> void replacecol(int j, const Matrix<T,storage>& b)
	{
#ifdef DEBUG
		if(b.cols() != 1 || b.rows() != nrows) { std::cout << "\nSize error in replacecol"; return; }
#endif
		for(int i = 0; i < nrows; i++)
			elems[i*ncols+j] = b.get(i);
	}

	/// Function for replacing a row
	template <MStype storage> void replacerow(int i, Matrix<T,storage> b)
	{
#ifdef DEBUG
		if(b.cols() != ncols || b.rows() != 1) { std::cout << "\nSize error in replacerow"; return; }
#endif
		for(int j = 0; j < ncols; j++)
			elems[i*ncols + j] = b.get(0,j);
	}

	/// Returns the transpose
	template <MStype storage> Matrix<T,storage> trans() const
	{
		Matrix<T,storage> t(ncols, nrows);
		for(int i = 0; i < ncols; i++)
			for(int j = 0; j < nrows; j++)
				t(i,j) = get(j,i);
		return t;
	}

	/// Multiply a matrix by a scalar. Note: only expressions of type A*3 work, not 3*A
	template <MStype storage> Matrix<T,storage> operator*(T num) const
	{
		Matrix<T,storage> A(nrows,ncols);
		int i,j;

		for(i = 0; i < nrows; i++)
			for(j = 0; j < ncols; j++)
				A(i,j) = get(i,j) * num;
		return A;
	}

	template <MStype storage1, MStype storage2> Matrix<T,storage1> operator*(const Matrix<T,storage2>& B) const
	{
		Matrix<T, storage1> C(nrows, B.cols());
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
					C(i,j) += get(i,k) * B.get(k,j);

		return C;
	}

	/// Returns sum of products of respective elements of flattened arrays containing matrix elements of this and A
	template <MStype storage> T dot_product(const Matrix<T,storage>& A) const
	{
#		ifdef DEBUG
		if(A.rows() != nrows || A.cols() != ncols)
			std::cout << "! Matrix rowmajor: Size mismatch!" << std::endl;
#		endif
		int i,j;
		double ans = 0;
		for(i = 0; i < nrows; i++)
			for(j = 0; j < ncols; j++)
			{
				T temp = get(i,j)*A.get(i,j);
				ans += temp;
			}
		return ans;
	}

	/// Computes 1-norm (max column-sum norm) of the matrix
	T matrixNorm_1() const
	{
		T max = 0, sum;
		int i,j;
		for(j = 0; j < ncols; j++)
		{
			sum = 0;
			for(i = 0; i < nrows; i++)
			{
				sum += (T)( fabs(get(i,j)) );
			}
			if(max < sum) max = sum;
		}
		return max;
	}
};

/// Specialization of the Matrix class to column-major storage
template <class T> 
class Matrix<T, COLMAJOR>
{
private:
	int nrows;
	int ncols;
	int size;
	T* elems;
	bool isalloc;

public:
	///No-arg constructor. Note: no memory allocation! Make sure Matrix::setup(int,int,MStype) is used.
	Matrix()
	{
		nrows = 0; ncols = 0; size = 0;
		isalloc = false;
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
		isalloc = true;
	}

	/// copy-initialize from another column-major matrix
	Matrix(const Matrix<T,COLMAJOR>& other)
	{
		nrows = other.nrows;
		ncols = other.ncols;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows*ncols; i++)
		{
			elems[i] = other.elems[i];
		}
	}
	
	/// Copy-initialize from a row-major matrix; not very efficient
	Matrix(const Matrix<T,ROWMAJOR>& other)
	{
		nrows = other.nrows;
		ncols = other.ncols;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				elems[j*nrows+i] = other.elems[i*ncols+j];
	}

	~Matrix()
	{
		if(isalloc == true)	
			delete [] elems;
		isalloc = false;
	}

	Matrix<T,COLMAJOR>& operator=(const Matrix<T,COLMAJOR>& rhs)
	{
#ifdef DEBUG
		if(this==&rhs) return *this;		// check for self-assignment
#endif
		nrows = rhs.nrows;
		ncols = rhs.ncols;
		size = nrows*ncols;
		if(isalloc == true)
			delete [] elems;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows*ncols; i++)
			elems[i] = rhs.elems[i];
		return *this;
	}
	
	Matrix<T,COLMAJOR>& operator=(const Matrix<T,ROWMAJOR>& rhs)
	{
#ifdef DEBUG
		if(this==&rhs) return *this;		// check for self-assignment
#endif
		nrows = rhs.nrows;
		ncols = rhs.ncols;
		size = nrows*ncols;
		if(isalloc == true)
			delete [] elems;
		elems = new T[nrows*ncols];
		isalloc = true;
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				elems[j*nrows+i] = rhs.elems[i*ncols+j];
		return *this;
	}

	/// Separate setup function in case no-arg constructor has to be used
	void setup(int nr, int nc)
	{
		if(nc==0)
		{
			std::cout << "Matrix: setup(): Error: Number of columns is zero. Setting it to 1.\n";
			nc=1;
		}
		if(nr==0)
		{
			std::cout << "Matrix(): setup(): Error: Number of rows is zero. Setting it to 1.\n";
			nr=1;
		}
		nrows = nr; ncols = nc;
		size = nrows*ncols;
		if(isalloc == true)
			delete [] elems;
		elems = new T[nrows*ncols];
		isalloc = true;
	}

	/// Setup without deleting earlier allocation: use in case of Matrix<t>* (pointer to Matrix<t>)
	void setupraw(int nr, int nc)
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
		nrows = nr; ncols = nc;
		size = nrows*ncols;
		elems = new T[nrows*ncols];
		isalloc = true;
	}
	
	/// Fill the matrix with zeros.
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

	/// function to set matrix elements from a ROW-MAJOR array
	void setdata_rowmajor(const T* A, int sz)
	{
#ifdef DEBUG
		if(sz != size)
		{
			std::cout << "\nError in setdata: argument size does not match matrix size";
			return;
		}
#endif
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				elems[j*nrows+i] = A[i*ncols+j];
	}

	T get(const int i, const int j=0) const
	{
#ifdef DEBUG
		if(i>=nrows || j>=ncols) { std::cout << "Matrix: get(): Index beyond array size(s)\n"; return 0; }
		if(i < 0 || j < 0) { std::cout << "Matrix: get(): Index less than 0!\n"; return 0; }
#endif
		return elems[j*nrows + i];
	}

	void set(const int i, const int j, const T data)
	{
#ifdef DEBUG
		if(i>=nrows || j>=ncols) { std::cout << "Matrix: set(): Index beyond array size(s)\n"; return; }
		if(i < 0 || j < 0) {std::cout << "Matrix: get(): Negative index!\n"; return 0; }
#endif
		elems[j*nrows + i] = data;
	}

	int rows() const { return nrows; }
	int cols() const { return ncols; }
	int msize() const { return size; }

	/// Prints the matrix to standard output.
	void mprint() const
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << elems[j*nrows+i];
			std::cout << std::endl;
		}
	}

	/// Prints the matrix to file
	void fprint(std::ofstream& outfile) const
	{
		outfile << std::setprecision(MATRIX_DOUBLE_PRECISION);
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << elems[j*nrows+i];
			outfile << '\n';
		}
	}

	/// Reads matrix from file
	void fread(std::ifstream& infile)
	{
		infile >> nrows; infile >> ncols;
		size = nrows*ncols;
		delete [] elems;
		elems = new T[nrows*ncols];
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				infile >> elems[j*nrows + i];
	}

	/// Getter/setter function for expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141
	T& operator()(const int x, const int y=0)
	{
#ifdef DEBUG
		if(x>=nrows || y>=ncols) { std::cout << "Matrix (): Index beyond array size(s)\n"; return elems[0]; }
#endif
		return elems[y*nrows + x];
	}
	
	/// Const Getter/setter function for expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141
	const T& operator()(const int x, const int y=0) const
	{
#ifdef DEBUG
		if(x>=nrows || y>=ncols) { std::cout << "Matrix (): Index beyond array size(s)\n"; return elems[0]; }
#endif
		return elems[y*nrows + x];
	}

	T maxincol(int j) const
	{
		T max = get(0,j);
		for(int i = 0; i < nrows; i++)
			if(max < get(i,j)) max = get(i,j);
		return max;
	}

	T maxinrow(int i) const
	{
		T max = get(i,0);
		for(int j = 0; j < nrows; j++)
			if(max < get(i,j)) max = get(i,j);
		return max;
	}

	T max() const
	{
		T max = elems[0];
		for(int i = 0; i < size; i++)
			if(elems[i] > max) max = elems[i];
		return max;
	}

	T absmax() const
	{
		T max = abs(elems[0]);
		for(int i = 0; i < size; i++)
			if(abs(elems[i]) > max) max = abs(elems[i]);
		return max;
	}

	/// Returns the magnitude of the element with largest magnitude
	double dabsmax() const
	{
		double max = dabs((double)elems[0]);
		for(int i = 0; i < size; i++)
			if(dabs(elems[i]) > max) max = dabs(elems[i]);
		return max;
	}

	T minincol(int j) const
	{
		T min = get(0,j);
		for(int i = 0; i < nrows; i++)
			if(min > get(i,j)) min = get(i,j);
		return min;
	}

	T mininrow(int i) const
	{
		T max = get(i,0);
		for(int j = 0; j < nrows; j++)
			if(max > get(i,j)) max = get(i,j);
		return max;
	}

	T min() const
	{
		T max = elems[0];
		for(int i = 0; i < size; i++)
			if(elems[i] < max) max = elems[i];
		return max;
	}

	T average() const
	{
		T avg = 0;
		for(int i = 0; i < size; i++)
			avg += elems[i];
		avg = avg/size;
		return avg;
	}

	// sums the square of all elements in the matrix and returns the square root of this sum
	T l2norm() const		
	{
		T tot = 0;
		for(int i = 0; i < size; i++)
		{
			tot += elems[i]*elems[i];
		}
		tot = std::sqrt(tot);
		return tot;
	}

	/// function to return a sub-matrix of this matrix
	template <MStype storage> Matrix<T,storage> sub(int startr, int startc, int offr, int offc) const
	{
		Matrix<T,storage> B(offr, offc);
		for(int i = 0; i < offr; i++)
			for(int j = 0; j < offc; j++)
				B(i,j) = get(startr+i,startc+j);
		return B;
	}

	/// Function that returns a given column of the matrix as a row-major matrix
	template <MStype storage> Matrix<T,storage> col(int j) const
	{
		Matrix<T,storage> b(nrows, 1);
		for(int i = 0; i < nrows; i++)
			b(i,0) = elems[j*nrows + i];
		return b;
	}

	template <MStype storage> Matrix<T,storage> row(int i) const
	{
		Matrix<T,storage> b(1, ncols);
		for(int j = 0; j < ncols; j++)
			b(0,j) = elems[j*nrows + i];
		return b;
	}

	///Function to return a pointer to a given column of the matrix
	T* col_p(int j)
	{
		T* b = &elems[j*nrows];
		return *b;
	}

	/// Function for replacing a column of the matrix with a vector. NOTE: No check for whether b is really a vector - which it must be.
	template <MStype storage> void replacecol(const int j, const Matrix<T,storage>& b)
	{
#ifdef DEBUG
		if(b.cols() != 1 || b.rows() != nrows) { std::cout << "\nSize error in replacecol"; return; }
#endif
		for(int i = 0; i < nrows; i++)
			elems[j*nrows + i] = b.get(i,0);
	}

	/// Function for replacing a row; inefficient since this is a column-major array. b should be a single row
	template <MStype storage> void replacerow(int i, const Matrix<T,storage>& b)
	{
#ifdef DEBUG
		if(b.cols() != ncols || b.rows() != 1) { std::cout << "\nSize error in replacerow"; return; }
#endif
		for(int j = 0; j < ncols; j++)
			elems[j*nrows + i] = b.get(0,j);
	}

	/// Returns the transpose
	template <MStype storage> Matrix<T,storage> trans() const
	{
		Matrix<T,storage> t(ncols, nrows);
		for(int i = 0; i < ncols; i++)
			for(int j = 0; j < nrows; j++)
				t(i,j) = get(j,i);
		return t;
	}

	/// Multiply a matrix by a scalar. Note: only expressions of type A*3 work, not 3*A
	template <MStype storage> Matrix<T,storage> operator*(T num)
	{
		Matrix<T,storage> A(nrows,ncols);
		int i,j;
		for(i = 0; i < nrows; i++)
			for(j = 0; j < ncols; j++)
				A(i,j) = elems[j*nrows+i] * num;
		return A;
	}

	template <MStype storage1, MStype storage2> Matrix<T,storage1> operator*(const Matrix<T,storage2>& B) const
	{
		Matrix<T, storage1> C(nrows, B.cols());
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
					C(i,j) += get(i,k) * B.get(k,j);

		return C;
	}

	/// Returns sum of products of respective elements of flattened arrays containing matrix elements of this and A
	template <MStype storage> T dot_product(const Matrix<T,storage>& A) const
	{
#		ifdef DEBUG
		if(A.rows() != nrows || A.cols() != ncols)
			std::cout << "! Matrix rowmajor: Size mismatch!" << std::endl;
#		endif
		int i,j;
		double ans = 0;
		for(i = 0; i < nrows; i++)
			for(j = 0; j < ncols; j++)
			{
				T temp = get(i,j)*A.get(i,j);
				ans += temp;
			}
		return ans;
	}

	/// Computes 1-norm (max column-sum norm) of the matrix
	T matrixNorm_1() const
	{
		T max = 0, sum;
		int i,j;
		for(j = 0; j < ncols; j++)
		{
			sum = 0;
			for(i = 0; i < nrows; i++)
			{
				sum += (T)( fabs(get(i,j)) );
			}
			if(max < sum) max = sum;
		}
		return max;
	}
};

} //end namespace amat

#endif
