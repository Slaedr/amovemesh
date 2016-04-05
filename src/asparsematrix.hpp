/** \file asparsematrix.hpp
 * \brief Sparse storage schemes for matrices (header-only library)
 * \author Aditya Kashi
 * \date 21 July 2015
 */

#ifndef __ASPARSEMATRIX_H
#define __ASPARSEMATRIX_H

// for memcpy() etc
#ifndef _GLIBCXX_CSTRING
#include <cstring>
#endif

// for malloc()
#ifndef _GLIBCXX_CSTDLIB
#include <cstdlib>
#endif

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_FSTREAM
#include <fstream>
#endif

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#ifndef __AMATRIX2_H
#include "amatrix2.hpp"
#endif

#ifdef _OPENMP
#ifndef OMP_H
#include <omp.h>
#endif
#endif

#ifndef ZERO_TOL
#define ZERO_TOL 1e-15
#endif

namespace amat {

typedef Matrix<double> Mat;

/// Structure for matrix in compressed row storage
/** \note Set [allocated](@ref allocated) to "true" right after reserving memory for the arrays.
 * Use matrix type "SLU_NR" while using this with SuperLU.
 */
template <class T> struct SMatrixCRS
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
template <class T> class SparseMatrix
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

	virtual void multiply(const Mat& x, Mat* const a, const char paralel) const = 0;
	virtual void multiply_parts(const Mat* x, const Mat* y, Mat* const ans, const int p) const = 0;

	virtual double getelem_multiply_parts(const int rownum, const Mat* x, const Mat* y, const int p, const double num) const = 0;
};

#define INITIAL_ROW_SIZE 20

/** \brief Implements sparse matrix storage in row-storage format, but with separate arrays for each row of the matrix.
 */
template<class T> class MatrixCRS : public SparseMatrix<T>
{
	std::vector<T>* val;
	std::vector<int>* col_ind;			///< col_ind[k][j] contains the column number of val[k][j]
	std::vector<int> rsize;				///< Number of non-zero elements in each row

	using SparseMatrix<T>::nnz;
	using SparseMatrix<T>::nrows;
	using SparseMatrix<T>::ncols;

public:
	MatrixCRS() 
	{ 
		val = new std::vector<T>[2];  
		col_ind = new std::vector<int>[2];
	}

	MatrixCRS(int num_rows, int num_cols) : SparseMatrix<T>(num_rows,num_cols)
	{
		rsize.resize(nrows,0);
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];
		for(int i = 0; i < nrows; i++)
		{
			val[i].reserve(INITIAL_ROW_SIZE);
			col_ind[i].reserve(INITIAL_ROW_SIZE);
		}
	}

	MatrixCRS(const MatrixCRS& other)
	{
		rsize = other.rsize;
		nnz = other.nnz;
		nrows = other.nrows; ncols = other.ncols;

		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];

		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < rsize[i]; j++)
			{
				val[i].push_back(other.val[i][j]);
				col_ind[i].push_back(other.col_ind[i][j]);
			}
	}

	MatrixCRS& operator=(const MatrixCRS& other)
	{
		rsize = other.rsize;
		nnz = other.nnz;
		nrows = other.nrows; ncols = other.ncols;

		delete [] val;
		delete [] col_ind;
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];

		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < rsize[i]; j++)
			{
				val[i].push_back(other.val[i][j]);
				col_ind[i].push_back(other.col_ind[i][j]);
			}
		return *this;
	}

	void setup(int num_rows, int num_cols)
	{
		delete [] val;
		delete [] col_ind;

		nrows = num_rows; ncols = num_cols;
		rsize.resize(nrows,0);
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];
		for(int i = 0; i < nrows; i++)
		{
			val[i].reserve(INITIAL_ROW_SIZE);
			col_ind[i].reserve(INITIAL_ROW_SIZE);
		}
		nnz = 0;
	}

	~MatrixCRS()
	{
		delete [] val;
		delete [] col_ind;
	}

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
	void trim()
	{
		for(int i = 0; i < nrows; i++) {
			val[i].shrink_to_fit();
			col_ind[i].shrink_to_fit();
		}
	}


	/*void print_data()
	{
		for(int i = 0; i < nnz; i++)
			std::cout << val[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < nnz; i++)
			std::cout << row_ind[i] << " ";
		std::cout << std::endl;
		for(int i = 0; i < nnz; i++)
			std::cout << col_ind[i] << " ";
		std::cout << std::endl;
	}*/

	void mprint() const
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << get(i,j);
			std::cout << std::endl;
		}
	}

	void fprint(std::ofstream& outfile) const
	{
		//outfile << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				outfile << " " << get(i,j);
			outfile << std::endl;
		}
	}

	/**	Sorts each row by column index using insertion sort.
	 * We expect each row to have no more than tens of elements, so insertion sort is a good candidate.
	 */
	void sort_rows()
	{
		// insertion sort
		int tempc = 0; T tempv = 0; int j = 0;

		std::vector<int>* rsize = &(this->rsize);
		#ifdef _OPENMP
		std::vector<T>* val = this->val;
		std::vector<int>* col_ind = this->col_ind;
		int nrows = this->nrows;
		int ncols = this->ncols;
		#endif

		for(int l = 0; l < nrows; l++)			// for each row
		{
			for(int i = 1; i < (*rsize)[l]; i++)
			{
				j = i;
				while(j>0 && col_ind[l][j] < col_ind[l][j-1])
				{
					// swap col_ind[l][j] and col_ind[l][j-1], same for vp
					tempc = col_ind[l][j];
					col_ind[l][j] = col_ind[l][j-1];
					col_ind[l][j-1] = tempc;
					tempv = val[l][j];
					val[l][j] = val[l][j-1];
					val[l][j-1] = tempv;
					j--;
				}
			}
		}
	}

	/// Returns product of sparse matrix with std::vector x and stores it in a.
	void multiply(const Mat& x, Mat* const a, const char paralel = 'n') const
	{
		a->zeros();
		
		//std::cout << "MatrixCOO: multiply(): Rows and columns of a: " << a->rows() << " " << a->cols() << std::endl;
		double temp;
		#ifdef _OPENMP
		const std::vector<T>* val = this->val; const std::vector<int>* col_ind = this->col_ind;
		#endif

 		const std::vector<int>* rsize = &(this->rsize);

		for(int k = 0; k < x.cols(); k++)		// for each column of x
		{
			int i;
			//#pragma omp parallel for default(none) private(i,temp) shared(val,col_ind,rsize,k,a,x) if(paralel=='y')
			for(i = 0; i < nrows; i++)
			{
				for(int j = 0; j < (*rsize)[i]; j++)
				{
					temp = val[i][j]*x.get(col_ind[i][j],k);
					(*a)(i,k) += temp;
				}
			}
		}
	}

	/// Like the multiply() method, except the argument matrix is considered x for the first p-1 rows, 0 in the pth row and y for the remaining rows.
	void multiply_parts(const Mat* x, const Mat* y, Mat* const ans, const int p) const
	{
		ans->zeros();
		#ifdef _OPENMP
		const std::vector<T>* val = MatrixCRS::val; const std::vector<int>* col_ind = MatrixCRS::col_ind;
		#endif
		const std::vector<int>* rsize = &this->rsize;

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			int i;
			//#pragma omp parallel for default(none) private(i) shared(val,row_ind,col_ind,nnz,k,ans,x,y,p)
			for(i = 0; i < nrows; i++)
			{
				for(int j = 0; j < (*rsize)[i]; j++)
				{
					if(col_ind[i][j] < p)
						(*ans)(i,k) += val[i][j] * x->get(col_ind[i][j],k);
					else if(col_ind[i][j] > p)
						(*ans)(i,k) += val[i][j] * y->get(col_ind[i][j],k);
				}
			}
		}
	}

	/** Returns dot product of rownum'th row of this matrix with a certain std::vector. 
	 * This vector is composed of x for for the first p-1 elements and y for p+1'th element onwards, with num as pth element.
	 * \note NOTE: Make sure dimensions of x and y are same.
	 */
	T getelem_multiply_parts(const int rownum, const Mat* x, const Mat* y, const int p, const T num) const
	{
		T ans = 0;
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n";

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			#ifdef _OPENMP
			const std::vector<T>* val = MatrixCRS::val; const std::vector<int>* col_ind = MatrixCRS::col_ind;
			#endif
			const std::vector<int>* rsize = &this->rsize;
			int j;
			//#pragma omp parallel for default(none) private(j) shared(val,col_ind,rsize,k,rownum,p,x,y,num) reduction(+:ans)
			for(j = 0; j < (*rsize)[rownum]; j++)
			{
					if(col_ind[rownum][j] < p)
						ans += val[rownum][j] * x->get(col_ind[rownum][j],k);
					else if(col_ind[rownum][j] > p)
						ans += val[rownum][j] * y->get(col_ind[rownum][j],k);
					else
						ans += val[rownum][j] * num;
			}
		}
		return ans;
	}

	/// D is returned as a column-vector containing diagonal elements of this sparse matrix.
	void get_diagonal(Mat* D) const
	{
		D->zeros();
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < rsize[i]; j++)	
				if(i == col_ind[i][j])
				{
					//(*D)(row_ind[i]) = val[i];
					D->set(i,0, val[i][j]);
				}
		}
	}
	
	void get_lower_triangle(MatrixCRS& L)
	{
		delete [] L.val;
		delete [] L.col_ind;
		L.val = new std::vector<T>[nrows];
		L.col_ind = new std::vector<int>[nrows];
		
		L.nrows = nrows;
		L.ncols = ncols;
		L.rsize.resize(nrows,0);
		
		
		for(int i = 0; i < nrows; i++)
		{
		
			for(int j = 0; j < rsize[i]; j++)
				if(i > col_ind[i][j])
				{
					L.val[i].push_back(val[i][j]);
					L.col_ind[i].push_back(col_ind[i][j]);
					L.rsize[i]++;
				}
		}
		
		// shrink to fit
		/*if(L.nnz < nnz)
		{
			for(int i = nnz-1; i >= L.nnz; i--)
			{
				delete L.val+i;
			}
		}*/
	}
	
	void get_upper_triangle(MatrixCRS& L)
	{
		delete [] L.val;
		delete [] L.col_ind;
		L.val = new std::vector<T>[nrows];
		L.col_ind = new std::vector<int>[nrows];
		
		L.nrows = nrows;
		L.ncols = ncols;
		L.rsize.resize(nrows,0);
		
		
		for(int i = 0; i < nrows; i++)
		{
		
			for(int j = 0; j < rsize[i]; j++)
				if(i < col_ind[i][j])
				{
					L.val[i].push_back(val[i][j]);
					L.col_ind[i].push_back(col_ind[i][j]);
					L.rsize[i]++;
				}
		}
	}

	MatrixCRS<T> transpose()
	{
		std::cout << "MatrixCRS: transpose(): Transposing matrix" << std::endl;
		MatrixCRS<T> mat;
		delete [] mat.val;
		delete [] mat.col_ind;
		mat.val = new std::vector<T>[nrows];
		mat.col_ind = new std::vector<int>[nrows];
		mat.rsize.resize(nrows);
		mat.nrows = nrows; mat.ncols = ncols;

		int i;

		for(i = 0; i < nrows; i++)
		{
			for(int j = 0; j < rsize[i]; j++)
			{
				mat.val[col_ind[i][j]].push_back(val[i][j]);
				mat.col_ind[col_ind[i][j]].push_back(i);
				mat.rsize[col_ind[i][j]]++;
			}
		}

		return mat;
	}

	/// Combines A, B, C, and D (4 n-by-n matrices) into one 2n-by-2n matrix. 
	/** CAUTION: All 4 input matrices must have same size! */
	void combine_sparse_matrices(const MatrixCRS& A11, const MatrixCRS& A12, const MatrixCRS& A21, const MatrixCRS& A22)
	{
		std::cout << "MatrixCOO: combine_sparse_matrices(): Combining matrices" << std::endl;
		//std::cout << "MatrixCOO: combine_sparse_matrices(): A11 dimensions: " << A11.nrows << " " << A11.ncols << std::endl;

		nrows = A11.nrows*2; ncols = A11.ncols*2;
		delete [] val; delete [] col_ind;
		val = new std::vector<T>[nrows];
		col_ind = new std::vector<int>[nrows];

		// reserve some amount of space
		int maxrowsize = 6;
		for(int ir = 0; ir < nrows; ir++)
			if(A11.rsize[ir] > maxrowsize) maxrowsize = A11.rsize[ir];

		for(int ir = 0; ir < nrows; ir++) {
			val[ir].reserve(INITIAL_ROW_SIZE*2);
			col_ind[ir].reserve(INITIAL_ROW_SIZE*2);
		}

		rsize.resize(nrows,0);

		// start filling up the matrix

		int i;
		for(i = 0; i < A11.nrows; i++)
		{
			for(int j = 0; j < A11.rsize[i]; j++)
			{
				val[i].push_back(A11.val[i][j]);
				col_ind[i].push_back(A11.col_ind[i][j]);
				rsize[i]++;
			}
		}

		for(i = 0; i < A12.nrows; i++)
		{
			for(int j = 0; j < A12.rsize[i]; j++)
			{
				val[i].push_back(A12.val[i][j]);
				col_ind[i].push_back(A12.col_ind[i][j]+nrows/2);
				rsize[i]++;
			}
		}

		for(i = 0; i < A21.nrows; i++)
		{
			for(int j = 0; j < A21.rsize[i]; j++)
			{
				val[i+nrows/2].push_back(A21.val[i][j]);
				col_ind[i+nrows/2].push_back(A21.col_ind[i][j]);
				rsize[i+nrows/2]++;
			}
		}

		for(i = 0; i < A22.nrows; i++)
		{
			for(int j = 0; j < A22.rsize[i]; j++)
			{
				val[i+nrows/2].push_back(A22.val[i][j]);
				col_ind[i+nrows/2].push_back(A22.col_ind[i][j]+nrows/2);
				rsize[i+nrows/2]++;
			}
		}
	}

	///	Converts to dense matrix.
	Mat toDense() const
	{
		Matrix<T> dense(nrows, ncols);
		for(int i = 0; i < nrows; i++)
			for(int j = 0; j < rsize[i]; j++)
				dense(i,col_ind[i][j]) = val[i][j];
		return dense;
	}

	/// Returns the sparse matrix in compressed row format
	void get_CRS_matrix(SMatrixCRS<T>& rmat) const
	{
		rmat.nnz = 0;
		for(int i = 0; i < nrows; i++)
			rmat.nnz += rsize[i];

		if(rmat.allocated)
		{
			delete [] rmat.val;
			delete [] rmat.col_ind;
			delete [] rmat.row_ptr;
		}

		rmat.val = new T[rmat.nnz];
		rmat.col_ind = new int[rmat.nnz];
		rmat.row_ptr = new int[nrows+1];
		rmat.allocated = true;
		rmat.row_ptr[0] = 0;

		int i, j, k = 0;
		for(i = 0; i < nrows; i++)
		{
			for(j = 0; j < rsize[i]; j++)
			{
				rmat.val[k] = val[i][j];
				rmat.col_ind[k] = col_ind[i][j];
				k++;
			}
			rmat.row_ptr[i+1] = k;
		}
	}

#	ifdef COMPILE_STUPID_STUFF
	/// Creates an Eigen3 sparse matrix in row major format
	void get_Eigen3_rowmajor_matrix( Eigen::MappedSparseMatrix<T, Eigen::RowMajor>& A ) const
	{
		SMatrixCRS<T> rmat;
		rmat.nnz = 0;
		for(int i = 0; i < nrows; i++)
			rmat.nnz += rsize[i];

		rmat.val = new T[rmat.nnz];
		rmat.col_ind = new int[rmat.nnz];
		rmat.row_ptr = new int[nrows+1];
		rmat.allocated = true;
		rmat.row_ptr[0] = 0;

		int i, j, k = 0;
		for(i = 0; i < nrows; i++)
		{
			for(j = 0; j < rsize[i]; j++)
			{
				rmat.val[k] = val[i][j];
				rmat.col_ind[k] = col_ind[i][j];
				k++;
			}
			rmat.row_ptr[i+1] = k;
		}

		/*int i, j, k = 0;
		for(i = 0; i < nrows; i++)
		{
			for(j = 0; j < rsize[i]; j++)
			{
				A.valuePtr()[k] = val[i][j];
				A.innerIndexPtr()[k] = col_ind[i][j];
				k++;
			}
			A.outerIndexPtr()[i+1] = k;
		}*/
	}
#	endif

	/* Computes the LU factorization with partial pivoting.
	 * L is the unit lower triangular matrix, U is the upper triangular matrix.
	 */
	/*void LUfactor(MatrixCRS<double>& L, MatrixCRS<double>& U, MatrixCRS<int>& P)
	{
		// copy this matrix (A) into U
		delete [] U.val;
		delete [] U.col_ind;
		U.val = new vector<double>[nrows];
		U.col_ind = new vector<int>[ncols];
		U.rsize = rsize;
		for(int i = 0; i < nrows; i++)
		{
			U.val[i].resize(rsize[i]);
			U.col_ind[i].resize(rsize[i]);
			for(int j = 0; j < rsize[i]; j++)
			{
				U.val[i][j] = val[i][j];
				U.col_ind[i][j] = col_ind[i][j];
			}
		}

		// Make L and P identity matrices
		delete [] L.val;
		delete [] L.col_ind;
		delete [] P.val;
		delete [] P.col_ind;
		L.val = new vector<double>[nrows];
		L.col_ind = new vector<int>[nrows];
		P.val = new vector<int>[nrows];
		P.col_ind = new vector<int>[nrows];
		for(int i = 0; i < nrows; i++)
		{
			L.rsize.push_back(1);
			P.rsize.push_back(1);
			L.val[i].resize(1);
			L.col_ind[i].resize(1);
			P.val[i].resize(1);
			P.col_ind[i].resize(1);

			L.val[i].reserve(10);
			L.col_ind[i].reserve(10);

			for(int j = 0; j < rsize[i]; j++)
			{
				L.val[i][j] = 1.0;
				L.col_ind[i][j] = i;
				P.val[i][j] = 1;
				P.col_ind[i][j] = i;
			}
		}

		// We can now start Gaussian elimination with partial pivoting
		double max = 0;
		int maxi = 0;
		for(int k = 0; k < nrows-1; k++)
		{
			// first find maximum element in the kth column.
			max = U.val[k][0];
			maxi = 0;
			for(int i = k+1; i < nrows; i++)
			{
				// we need a binary search of the ith row to look for col_ind k.
			}
		}
	}*/
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
	MatrixCRS_traditional(int num_rows, int num_cols) : SparseMatrix<double>(num_rows, num_cols)
	{
		val = new double[num_rows*2];
		col_ind = new int[num_rows*2];
		valsize = num_rows*2;
		for(int i = 0; i < valsize; i++)
		{
			val[i] = 0.0;
			col_ind[i] = -1;
		}
		row_ptr = new int[nrows+1];
		for(int i = 0; i < nrows; i++) row_ptr[i] = -1;
		row_ptr[nrows+1] = 0;
		nnz = 0;
	}
	~MatrixCRS_traditional()
	{
		delete [] val;
		delete [] col_ind;
		delete [] row_ptr;
	}

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
typedef MatrixCRS<double> SpMatrix;
#define SPARSE_MATRIX_TO_USE 1
#endif

}	// end namespace
#endif
