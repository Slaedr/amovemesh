/**
@file aarrray2d.hpp
@brief This file defines a class Array2d to store a 2D array in column-major or row-major storage.
@author Aditya Kashi
@date November 16, 2015
Notes:
This class is almost a copy of the Matrix class.
Matrix multiplication has been removed; it was inefficient anyway.
If A is a column-major matrix, A[i][j] == A[j * nrows + i] where i is the row-index and j is the column index.
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

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifdef _OPENMP
	#ifndef OMP_H
		#include <omp.h>
		#define nthreads_m 4
	#endif
#endif

#ifndef __ACONSTANTS_H
#include <aconstants.h>
#endif

#define __AARRAY2D_H

#ifndef MATRIX_DOUBLE_PRECISION
#define MATRIX_DOUBLE_PRECISION 14
#endif

const int WIDTH = 10;		// width of field for printing matrices

namespace amat {

double dabs(double x)
{
	if(x < 0) return (-1.0)*x;
	else return x;
}
double minmod(double a, double b)
{
	if(a*b>0 && dabs(a) <= dabs(b)) return a;
	else if (a*b>0 && dabs(b) < dabs(a)) return b;
	else return 0.0;
}

enum MStype {ROWMAJOR, COLMAJOR};

template <class T>
class Array2d
{
private:
	int nrows;
	int ncols;
	int size;
	MStype storage;
	T* elems;
	bool isalloc;

public:
	/// No-arg constructor. Note: no memory allocation! Make sure Array2d<T>::setup(int,int,MStype) is used.
	Array2d();

	/// Full-arg constructor.
	Array2d(int nr, int nc, MStype st = ROWMAJOR);

	Array2d(const Array2d<T>& other);

	~Array2d();

	Array2d<T>& operator=(Array2d<T> rhs);

	/// Separate setup function in case no-arg constructor has to be used.
	void setup(int nr, int nc, MStype st=ROWMAJOR);

	/// Same as [setup](@ref setup), but without deleting earlier allocation: use in case of Array2d<t>* (pointer to Array2d<t>)
	void setupraw(int nr, int nc, MStype st);

	void zeros();

	void ones();

	void identity();

	/// Function to set matrix elements from a ROW-MAJOR std::vector
	void setdata(const T* A, int sz);

	inline T get(int i, int j=0) const;

	void set(int i, int j, T data);

	int rows() const;
	int cols() const;
	int msize() const;
	MStype storetype() const;

	void mprint() const;

	void fprint(std::ofstream& outfile) const;

	void fread(std::ifstream& infile);

	/// For expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141.
	T& operator()(int x, int y=0);

	T maxincol(int j) const;

	T maxinrow(int i) const;

	T max() const;

	T absmax() const;

	double dabsmax() const;

	T minincol(int j) const;

	T mininrow(int i) const;

	T min() const;

	T average() const;

	T l2norm() const;		// sums the square of all elements in the matrix and returns the square root of this sum

	/// function to return a sub-matrix of this matrix
	Array2d<T> sub(int startr, int startc, int offr, int offc) const;

	/// Function that returns a given column of the matrix as a row-major matrix
	Array2d<T> col(int j) const;

	Array2d<T> row(int i) const;

	/*//Function to return a reference to a given column of the matrix
	Array2d<T>& colr(int j); */

	/// Function for replacing a column of the matrix with a vector. 
	/** NOTE: No check for whether b is really a vector - which it must be.
	*/
	void replacecol(int j, Array2d<T> b);

	/// Function for replacing a row
	void replacerow(int i, Array2d<T> b);

	/// Computes transpose.
	Array2d<T> trans() const;

	/// Multiply a matrix by a scalar. Note: only expressions of type A*3 work, not 3*A
	Array2d<T> operator*(T num);

	/**	The matrix addition and subtraction operators are inefficient! Do not use in long loops. */
	Array2d<T> operator+(Array2d<T> B) const;

	Array2d<T> operator-(Array2d<T> B) const;

	/// Returns sum of products of respective elements of flattened arrays containing matrix elements of this and A.
	T dot_product(const Array2d<T>& A);
};

//No-arg constructor. Note: no memory allocation! Make sure Array2d<T>::setup(int,int,MStype) is used.
template <typename T>
Array2d<T>::Array2d()
{
	nrows = 0; ncols = 0; size = 0;
	isalloc = false;
}

// Full-arg constructor
template <typename T>
Array2d<T>::Array2d(int nr, int nc, MStype st)
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
	elems = new T[nrows*ncols];
	isalloc = true;
}

template <typename T>
Array2d<T>::Array2d(const Array2d<T>& other)
{
	nrows = other.nrows;
	ncols = other.ncols;
	storage = other.storage;
	size = nrows*ncols;
	elems = new T[nrows*ncols];
	isalloc = true;
	for(int i = 0; i < nrows*ncols; i++)
	{
		elems[i] = other.elems[i];
	}
}

template <typename T>
Array2d<T>::~Array2d()
{
	if(isalloc == true)	
		delete [] elems;
	isalloc = false;
}

template <typename T>
Array2d<T>& Array2d<T>::operator=(Array2d<T> rhs)
{
#ifdef DEBUG
	if(this==&rhs) return *this;		// check for self-assignment
#endif
	nrows = rhs.nrows;
	ncols = rhs.ncols;
	storage = rhs.storage;
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

//Separate setup function in case no-arg constructor has to be used
template <typename T>
void Array2d<T>::setup(int nr, int nc, MStype st)
{
	if(nc==0)
	{
		std::cout << "Array2d: setup(): Error: Number of columns is zero. Setting it to 1.\n";
		nc=1;
	}
	if(nr==0)
	{
		std::cout << "Array2d(): setup(): Error: Number of rows is zero. Setting it to 1.\n";
		nr=1;
	}
	nrows = nr; ncols = nc; storage = st;
	size = nrows*ncols;
	if(isalloc == true)
		delete [] elems;
	elems = new T[nrows*ncols];
	isalloc = true;
}

// no deleting earlier allocation: use in case of Array2d<t>* (pointer to Array2d<t>)
template <typename T>
void Array2d<T>::setupraw(int nr, int nc, MStype st)
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
	isalloc = true;
}

template <typename T>
inline void Array2d<T>::zeros()
{
	for(int i = 0; i < size; i++)
		elems[i] = (T)(0.0);
}

template <typename T>
void Array2d<T>::ones()
{
	for(int i = 0; i < size; i++)
		elems[i] = 1;
}

template <typename T>
void Array2d<T>::identity()
{
	T one = (T)(1);
	T zero = (T)(0);
	for(int i = 0; i < nrows; i++)
		for(int j = 0; j < ncols; j++)
			if(i==j) operator()(i,j) = one;
			else operator()(i,j) = zero;
}

// function to set matrix elements from a ROW-MAJOR std::vector
template <typename T>
void Array2d<T>::setdata(const T* A, int sz)
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

template <typename T>
inline T Array2d<T>::get(int i, int j) const
{
#ifdef DEBUG
	if(i>=nrows || j>=ncols) { std::cout << "Array2d: get(): Index beyond array size(s)\n"; return 0; }
	if(i < 0 || j < 0) { std::cout << "Array2d: get(): Index less than 0!\n"; return 0; }
#endif
	if(storage == COLMAJOR)
		return elems[j*nrows + i];
	else
		return elems[i*ncols + j];
}

template <typename T>
void Array2d<T>::set(int i, int j, T data)
{
#ifdef DEBUG
	if(i>=nrows || j>=ncols) { std::cout << "Array2d: set(): Index beyond array size(s)\n"; return; }
#endif
	//if(i < 0 || j < 0) {std::cout << "Array2d: get(): Negative index!\n"; return 0; }
	if(storage == COLMAJOR)
		elems[j*nrows + i] = data;
	else
		elems[i*ncols + j] = data;
}

template <typename T>
inline int Array2d<T>::rows() const { return nrows; }
template <typename T>
inline int Array2d<T>::cols() const { return ncols; }
template <typename T>
inline int Array2d<T>::msize() const { return size; }
template <typename T>
MStype Array2d<T>::storetype() const { return storage;}

template <typename T>
void Array2d<T>::mprint() const
{
	std::cout << "\n";
	if(storage==COLMAJOR)
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << elems[j*nrows+i];
			std::cout << std::endl;
		}
	else
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << std::setprecision(WIDTH/2+1) << elems[i*ncols+j];
			std::cout << std::endl;
		}
}

template <typename T>
void Array2d<T>::fprint(std::ofstream& outfile) const
{
	//outfile << '\n';
	outfile << std::setprecision(MATRIX_DOUBLE_PRECISION);
	if(storage==COLMAJOR)
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << elems[j*nrows+i];
			outfile << '\n';
		}
	else
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << elems[i*ncols+j];
			outfile << '\n';
		}
}

template <typename T>
void Array2d<T>::fread(std::ifstream& infile)
{
	infile >> nrows; infile >> ncols;
	size = nrows*ncols;
	storage = ROWMAJOR;
	delete [] elems;
	elems = new T[nrows*ncols];
	if(storage==ROWMAJOR)
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				infile >> elems[i*ncols + j];
	else
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < ncols; j++)
				infile >> elems[j*nrows + i];
}

// For expressions like A(1,2) = 141 to set the element at 1st row and 2nd column to 141
template <typename T>
inline T& Array2d<T>::operator()(int x, int y)
{
#ifdef DEBUG
	if(x>=nrows || y>=ncols) { std::cout << "Array2d (): Index beyond array size(s)\n"; return elems[0]; }
#endif
	if(storage == COLMAJOR) return elems[y*nrows + x];
	else return elems[x*ncols + y];
}

template <typename T>
T Array2d<T>::maxincol(int j) const
{
	T max = get(0,j);
	for(int i = 0; i < nrows; i++)
		if(max < get(i,j)) max = get(i,j);
	return max;
}

template <typename T>
T Array2d<T>::maxinrow(int i) const
{
	T max = get(i,0);
	for(int j = 0; j < nrows; j++)
		if(max < get(i,j)) max = get(i,j);
	return max;
}

template <typename T>
T Array2d<T>::max() const
{
	T max = elems[0];
	for(int i = 0; i < size; i++)
		if(elems[i] > max) max = elems[i];
	return max;
}

template <typename T>
T Array2d<T>::absmax() const
{
	T max = abs(elems[0]);
	for(int i = 0; i < size; i++)
		if(abs(elems[i]) > max) max = abs(elems[i]);
	return max;
}

template <typename T>
inline double Array2d<T>::dabsmax() const
{
	double max = dabs((double)elems[0]);
	for(int i = 0; i < size; i++)
		if(dabs(elems[i]) > max) max = dabs(elems[i]);
	return max;
}

template <typename T>
T Array2d<T>::minincol(int j) const
{
	T min = get(0,j);
	for(int i = 0; i < nrows; i++)
		if(min > get(i,j)) min = get(i,j);
	return min;
}

template <typename T>
T Array2d<T>::mininrow(int i) const
{
	T max = get(i,0);
	for(int j = 0; j < nrows; j++)
		if(max > get(i,j)) max = get(i,j);
	return max;
}

template <typename T>
T Array2d<T>::min() const
{
	T max = elems[0];
	for(int i = 0; i < size; i++)
		if(elems[i] < max) max = elems[i];
	return max;
}

template <typename T>
T Array2d<T>::average() const
{
	T avg = 0;
	for(int i = 0; i < size; i++)
		avg += elems[i];
	avg = avg/size;
	return avg;
}

template <typename T>
T Array2d<T>::l2norm() const		// sums the square of all elements in the matrix and returns the square root of this sum
{
	T tot = 0;
	for(int i = 0; i < size; i++)
	{
		tot += elems[i]*elems[i];
	}
	tot = std::sqrt(tot);
	return tot;
}

// function to return a sub-matrix of this matrix
template <typename T>
Array2d<T> Array2d<T>::sub(int startr, int startc, int offr, int offc) const
{
	Array2d<T> B(offr, offc);
	if(storage == COLMAJOR)
		for(int i = 0; i < offr; i++)
			for(int j = 0; j < offc; j++)
				B(i,j) = elems[(startc+j)*nrows + startr + i];
	else
		for(int i = 0; i < offr; i++)
			for(int j = 0; j < offc; j++)
				B(i,j) = elems[(startr+i)*ncols + startc + j];
	return B;
}

//Function that returns a given column of the matrix as a row-major matrix
template <typename T>
Array2d<T> Array2d<T>::col(int j) const
{
	Array2d<T> b(nrows, 1);
	if(storage==COLMAJOR)
		for(int i = 0; i < nrows; i++)
			b(i,0) = elems[j*nrows + i];
	else
		for(int i = 0; i < nrows; i++)
			b(i,0) = elems[i*ncols + j];
	return b;
}

template <typename T>
Array2d<T> Array2d<T>::row(int i) const
{
	Array2d<T> b(1, ncols);
	if(storage==COLMAJOR)
		for(int j = 0; j < ncols; j++)
			b(0,j) = elems[j*nrows + i];
	else
		for(int j = 0; j < ncols; j++)
			b(0,j) = elems[i*ncols + j];
	return b;
}

/*//Function to return a reference to a given column of the matrix
template <typename T>
Array2d<T>& Array2d<T>::colr(int j)
{
	//Array2d<T>* b(nrows, 1);
	Array2d<T>* b; b->elems.reserve(nrows);
	if(storage==COLMAJOR)
		for(int i = 0; i < nrows; i++)
			b->elems[i] = elems[j*nrows + i];
	else
		for(int i = 0; i < nrows; i++)
			b.elems[i] = &elems[i*ncols + j];
	return *b;
} */

// Function for replacing a column of the matrix with a vector. NOTE: No check for whether b is really a vector - which it must be.
template <typename T>
void Array2d<T>::replacecol(int j, Array2d<T> b)
{
#ifdef DEBUG
	if(b.cols() != 1 || b.rows() != nrows) { std::cout << "\nSize error in replacecol"; return; }
#endif
	if(storage == COLMAJOR)
		for(int i = 0; i < nrows; i++)
			elems[j*nrows + i] = b.elems[i];
	if(storage == ROWMAJOR)
		for(int i = 0; i < nrows; i++)
			elems[i*ncols + j] = b.elems[i];
}

//Function for replacing a row
template <typename T>
void Array2d<T>::replacerow(int i, Array2d<T> b)
{
#ifdef DEBUG
	if(b.cols() != ncols || b.rows() != 1) { std::cout << "\nSize error in replacerow"; return; }
#endif
	if(storage == COLMAJOR)
		for(int j = 0; j < ncols; j++)
			elems[j*nrows + i] = b.elems[j];
	else
		for(int j = 0; j < ncols; j++)
			elems[i*ncols + j] = b.elems[j];
}

//transpose
template <typename T>
Array2d<T> Array2d<T>::trans() const
{
	Array2d<T> t(ncols, nrows, storage);
	for(int i = 0; i < ncols; i++)
		for(int j = 0; j < nrows; j++)
			t(i,j) = get(j,i);
	return t;
}

// Multiply a matrix by a scalar. Note: only expressions of type A*3 work, not 3*A
template <typename T>
Array2d<T> Array2d<T>::operator*(T num)
{
	Array2d<T> A(nrows,ncols,storage);
	int i;

	for(i = 0; i < A.size; i++)
		A.elems[i] = elems[i] * num;
	return A;
}

/**	The matrix addition and subtraction operators are inefficient! Do not use in long loops. */
template <typename T>
Array2d<T> Array2d<T>::operator+(Array2d<T> B) const
{
#ifdef DEBUG
	if(nrows != B.rows() || ncols != B.cols())
	{
		std::cout << "! Array2d: Addition cannot be performed due to incompatible sizes\n";
		Array2d<T> C(1,1);
		return C;
	}
#endif
	Array2d<T> C(nrows, ncols);
	int i;

	for(i = 0; i < C.size; i++)
		C.elems[i] = elems[i] + B.elems[i];
	return C;
}

template <typename T>
Array2d<T> Array2d<T>::operator-(Array2d<T> B) const
{
#ifdef DEBUG
	if(nrows != B.rows() || ncols != B.cols())
	{
		std::cout << "! Array2d: Subtraction cannot be performed due to incompatible sizes\n";
		Array2d<T> C(1,1);
		return C;
	}
#endif
	Array2d<T> C(nrows, ncols);

	for(int i = 0; i < C.size; i++)
		C.elems[i] = elems[i] - B.elems[i];
	return C;
}

template <typename T>
T Array2d<T>::dot_product(const Array2d<T>& A)
/* Returns sum of products of respective elements of flattened arrays containing matrix elements of this and A */
{
	T* elemsA = A.elems;
	#ifdef _OPENMP
	T* elems = this->elems;
	int size = this->size;
	#endif
	int i;
	double ans = 0;
	//#pragma omp parallel for if(size >= 4) default(none) private(i) shared(elems,elemsA,size) reduction(+: ans) num_threads(nthreads_m)
	for(i = 0; i < size; i++)
	{
		T temp = elems[i]*elemsA[i];
		ans += temp;
	}
	return ans;
}
} //end namespace amat
