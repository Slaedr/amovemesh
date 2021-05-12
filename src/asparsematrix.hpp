/** \file asparsematrix.hpp
 * \brief Sparse storage schemes for matrices (header-only library)
 * \author Aditya Kashi
 * \date 21 July 2015
 */

#ifndef AMC_SPARSEMATRIX_H
#define AMC_SPARSEMATRIX_H

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "amatrix.hpp"

namespace amat {

typedef Matrix<amc_real> Mat;

/// The type of dense storage to use for vectors in matrix-vector multiplication, etc.
template <typename T>
using DenseMatrix = Matrix<T>;

/// Structure for matrix in compressed row storage
/** \note Set [allocated](@ref allocated) to "true" right after reserving memory for the arrays.
 * Use matrix type "SLU_NR" while using this with SuperLU.
 */
template <class T>
struct SMatrixCRS
{
	int nnz;			///< Number of non-zero entries in the matrix
	T* val;				///< Non-zero values of the matrix
	int* col_ind;		///< Column indices of non-zero values of the matrix
	int* row_ptr;		///< Element indices of val and col_ind in which each row begins
	bool allocated;		///< Stores whether the 3 arrays have been allocated

	SMatrixCRS()
	{
		allocated = false;
	}
	
	~SMatrixCRS()
	{
		if(allocated)
		{
			delete [] val;
			delete [] col_ind;
			delete [] row_ptr;
			allocated = false;
		}
	}
};


/// Abstract class for sparse storage of matrices
template <class T>
class SparseMatrix
{
protected:
	int nnz;
	int nrows;
	int ncols;
public:
	SparseMatrix() {}
	SparseMatrix(int num_rows, int num_cols) : nrows(num_rows), ncols(num_cols)
	{
		nnz = 0;
	}
	int rows() const { return nrows; }
	int cols() const { return ncols; }

	virtual void set(const int x, const int y, const double value) = 0;
	virtual T get(const int x, const int y) const = 0;
	virtual void mprint() const = 0;

	/// Multiplication of this sparse matrix with a row-major vector
	virtual void multiply(const Mat& x, Mat* const a, const char paralel) const = 0;
	virtual void multiply_parts(const Mat* x, const Mat* y, Mat* const ans, const int p) const = 0;
	virtual double getelem_multiply_parts(const int rownum, const Mat* x, const Mat* y, const int p, const double num) const = 0;
};

/** \brief Implements sparse matrix storage in row-storage format,
 * but with separate arrays for each row of the matrix.
 */
template<class T>
class MatrixCRS : public SparseMatrix<T>
{
	std::vector<T>* val;
	std::vector<int>* col_ind;			///< col_ind[k][j] contains the column number of val[k][j]
	std::vector<int> rsize;				///< Number of non-zero elements in each row

	using SparseMatrix<T>::nnz;
	using SparseMatrix<T>::nrows;
	using SparseMatrix<T>::ncols;

public:
	MatrixCRS();

	MatrixCRS(int num_rows, int num_cols);

	MatrixCRS(const MatrixCRS& other);

	MatrixCRS& operator=(const MatrixCRS& other);

	void setup(int num_rows, int num_cols);

	~MatrixCRS();

	void set(const int x, const int y, const T value)
	{
		if(dabs(double(value)) > ZERO_TOL)
		{
			// search row x for the value
			int pos = -1;
			for(int j = 0; j < rsize[x]; j++)
				if(col_ind[x][j] == y)
					pos = j;
			
			if(pos == -1)
			{
				//#pragma omp critical (omp_setval)
				{
					val[x].push_back(value);
					col_ind[x].push_back(y);
					rsize[x]++;
				}
			}
			else
			{
				val[x][pos] = value;
			}
			nnz++;
		}
	}
	

	T get(const int x, const int y) const
	{
		T retval = T(0);
		for(int j = 0; j < rsize[x]; j++)
			if(col_ind[x][j] == y)
				retval = val[x][j];
		return retval;
	}

	/// Combined getter/setter method.
	/// Usage of this method might lead to zeros being stored.
	T& operator()(const int x, const int y)
	{
		T* retval;
		// search row x for the value
		int pos = -1;
		for(int j = 0; j < rsize[x]; j++)
			if(col_ind[x][j] == y)
				pos = j;
		
		if(pos == -1)
		{
			val[x].push_back(T(0));
			col_ind[x].push_back(y);
			rsize[x]++;
			retval = &val[x][rsize[x]-1];
			//nnz++;
		}
		else
		{
			retval = &val[x][pos];
		}
		return *retval;
	}

	/** Shrinks allocated memory of val and col_ind to fit the current number of non-zero elements.*/
	void trim();

	void mprint() const;

	void fprint(std::ofstream& outfile) const;

	/**	Sorts each row by column index using insertion sort.
	 * We expect each row to have no more than tens of elements, so insertion sort is a good candidate.
	 */
	void sort_rows();

	/// Returns product of sparse matrix with vector x and stores it in a.
	void multiply(const Mat& x, Mat* const a, const char paralel = 'n') const;

	/// Like the multiply() method, except the argument matrix is considered x for the first p-1 rows, 0 in the pth row and y for the remaining rows.
	void multiply_parts(const Mat* x, const Mat* y, Mat* const ans, const int p) const;

	/** Returns dot product of rownum'th row of this matrix with a certain vector.
	 *
	 * This vector is composed of x for for the first p-1 elements and y for p+1'th element onwards,
	 * with num as pth element.
	 * \note NOTE: Make sure dimensions of x and y are same.
	 */
	T getelem_multiply_parts(const int rownum, const Mat* x, const Mat* y, const int p, const T num) const;

	/// D is returned as a column-vector containing diagonal elements of this sparse matrix.
	void get_diagonal(Mat* D) const;
	
	void get_lower_triangle(MatrixCRS& L);
	
	void get_upper_triangle(MatrixCRS& L);

	MatrixCRS<T> transpose();

	/// Combines A, B, C, and D (4 n-by-n matrices) into one 2n-by-2n matrix. 
	/** CAUTION: All 4 input matrices must have same size! */
	void combine_sparse_matrices(const MatrixCRS& A11, const MatrixCRS& A12, 
			const MatrixCRS& A21, const MatrixCRS& A22);

	///	Converts to dense matrix.
	Mat toDense() const;

	/// Returns the sparse matrix in compressed row format
	void get_CRS_matrix(SMatrixCRS<T>& rmat) const;

#	ifdef COMPILE_STUPID_STUFF
	/// Creates an Eigen3 sparse matrix in row major format
	void get_Eigen3_rowmajor_matrix( Eigen::MappedSparseMatrix<T, Eigen::RowMajor>& A ) const;
#	endif

};

class MatrixCRS_traditional : public SparseMatrix<double>
{
	double* val;
	int* col_ind;
	int* row_ptr;
	using SparseMatrix::nnz;
	using SparseMatrix::nrows;
	using SparseMatrix::ncols;
	int valsize;		// total size allocated to array val, and therefore also to array col_ind

public:
	MatrixCRS_traditional(int num_rows, int num_cols);

	~MatrixCRS_traditional();

	void set(int x, int y, double value)
	{
		if(col_ind[valsize-1] != -1)		// re-allocate if no space is available
		{
			double* oldval = val;
			int* oldcol = col_ind;
			val = new double[2*valsize];
			col_ind = new int[2*valsize];
			for(int i = 0; i < valsize; i++)		// copy data over
			{
				val[i] = oldval[i];
				col_ind[i] = oldcol[i];
			}
			for(int i = valsize; i < 2*valsize; i++)		// initialize extra space with dummy data
			{
				val[i] = 0;
				col_ind[i] = -1;
			}
			delete [] oldval;
			delete [] oldcol;
		}

		// if value is not zero

		int ind = row_ptr[x];
		int nind = row_ptr[x+1];
		if(ind == nind)				// only zeros in this row
		{
			for(int i = nind; i < nnz; i++)
			{
				val[i+1] = val[i];				// create space for this element
				col_ind[i+1] = col_ind[i];
			}
		}
	}

	double get(const int x, const int y)
	{
		int ind = row_ptr[x];
		//if(col_ind[ind] > y) return 0.0;
		int nind = row_ptr[x+1];
		for(int i = ind; i <= nind-1; i++)
		{
			if(col_ind[i]==y) return val[i];
		}
		return 0.0;
	}
};

#ifndef SPARSE_MATRIX_TO_USE

/** We use MatrixCRS as the default sparse matrix implementation.
 *
 * \note NOTE that the default sparse matrix type stores floating-point entries of type amc_real
 */
typedef MatrixCRS<amc::amc_real> SpMatrix;

#define SPARSE_MATRIX_TO_USE 1
#endif

}	// end namespace
#endif
