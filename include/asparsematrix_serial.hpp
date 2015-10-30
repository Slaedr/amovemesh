/* Implemetation of sparse storage schemes for matrices.
   Aditya Kashi
   21 July 2015
*/

#define __ASPARSEMATRIX_H

// for memcpy() etc
#ifndef _GLIBCXX_CSTRING
#include <cstring>
#endif

#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLBCXX_CMATH
#include <cmath>
#endif

#ifndef __AMATRIX2_H
#include "amatrix2.hpp"
#endif

namespace amat {

typedef Matrix<double> Mat;

class SparseMatrix
/* Abstract class for sparse storage of matrices */
{
protected:
	double* val;
	int nnz;
	int nrows;
	int ncols;
	int valsize;
public:
	SparseMatrix(int num_rows, int num_cols) : nrows(num_rows), ncols(num_cols)
	{
		nnz = 0;
	}
	int rows() { return nrows; }
	int cols() { return ncols; }
	virtual void set(int x, int y, double value) = 0;
	virtual double get(int x, int y) = 0;
	virtual void mprint() = 0;

	virtual Mat multiply(Mat& x) = 0;
	virtual void multiply_parts(Mat* x, Mat* y, Mat* ans, int p) = 0;
};

class MatrixCOO : public SparseMatrix
/*
 Implements sparse matrix storage in coordinate format. val contains non-zero values, col_ind and row_ind contain corresponding
 column and row indices; all three have length nnz after all values are stored.
*/
{
	int* col_ind;
	int* row_ind;
	bool sorted;		// a boolean flag that is true if row_ind is sorted in increasing order, and for each row index, col_ind is sorted in increasing order
public:
	MatrixCOO(int num_rows, int num_cols) : SparseMatrix(num_rows, num_cols)
	{
		// Size of the 3 arrays is initially 2*num_rows
		valsize = 2*nrows;
		val = new double[valsize];
		col_ind = new int[valsize];
		row_ind = new int[valsize];
		for(int i = 0; i < valsize; i++)
		{
			val[i] = 0;
			col_ind[i] = -1;
			row_ind[i] = -1;
		}
		sorted = false;
	}

	~MatrixCOO()
	{
		delete [] val;
		delete [] col_ind;
		delete [] row_ind;
	}

	void set(int x, int y, double value)
	{
		if(row_ind[valsize-1] != -1)		// if last location of row_ind is not unused, allocate more space
		{
			std::cout << "MatrixCOO: set(): Increasing size of arrays\n";
			/*double* tval = new double[valsize];
			int* trow_ind = new int[valsize];
			int* tcol_ind = new int[valsize];

			for(int i = 0; i < valsize; i++)
			{
				tval[i] = val[i];
				trow_ind[i] = row_ind[i];
				tcol_ind[i] = col_ind[i];
			}

			delete [] val;
			delete [] row_ind;
			delete [] col_ind;

			val = new double[2*valsize];
			row_ind = new int[2*valsize];
			col_ind = new int[2*valsize];

			for(int i = 0; i < valsize; i++)
			{
				val[i] = tval[i];
				row_ind[i] = trow_ind[i];
				col_ind[i] = tcol_ind[i];
			}*/

			double* tval = new double[2*valsize];
			int* trow_ind = new int[2*valsize];
			int* tcol_ind = new int[2*valsize];

			memcpy(tval, val, valsize*sizeof(double));
			memcpy(trow_ind, row_ind, valsize*sizeof(int));
			memcpy(tcol_ind, col_ind, valsize*sizeof(int));

			delete [] val;
			delete [] row_ind;
			delete [] col_ind;

			val = tval;
			row_ind = trow_ind;
			col_ind = tcol_ind;

			// set default values to newly created empty space
			for(int i = valsize; i < 2*valsize; i++)
			{
				val[i] = 0;
				row_ind[i] = -1;
				col_ind[i] = -1;
			}

			valsize *= 2;
		}
	
		// now add the new element if it's non-zero
		if(value != 0.0)
		{
			int index = -1;
			for(int i = 0; i < nnz; i++)
			{
				if(row_ind[i]==x)
				{
					if(col_ind[i]==y)
					{
						index = i;
						break;
					}
				}
			}
			if(index == -1)
			{
				val[nnz] = value;
				row_ind[nnz] = x;
				col_ind[nnz] = y;
				nnz++;
				sorted = false;
			}
			else
			{	// a non-zero element already exists in the position (x,y)
				val[index] = value;
			}
		}
	}

	double get(int x, int y)
	{
		//std::cout << "MatrixCOO: get(): nnz is " << nnz << '\n';
		int index = -1;
		if(sorted == false)
		{
			for(int i = 0; i < nnz; i++)
			{
				if(row_ind[i]==x)
				{
					if(col_ind[i]==y)
					{
						index = i;
						break;
					}
				}
			}
			if(index == -1)
			{
				//std::cout << "! MatrixCOO: get(): Element at " << x << ", " << y << " not found!\n";
				return 0;
			}
			else return val[index];
		}
		else
		{
			// binary search
			int left = 0, right = nnz-1; int mid = (left+right)/2;
			while(left <= right)
			{
				mid = (left+right)/2;
				if(row_ind[mid] == x) break;
				else if(row_ind[mid] < x) left = mid+1;
				else right = mid-1;
			}
			//std::cout << "* get(): " << mid << '\n';

			if(row_ind[mid] != x)
				return 0;

			int i;
			bool inflag1 = false, inflag2 = false;		// flags to prevent infinite jumping between consecutive elements
			for(i = mid; i >= 0 && i < nnz && row_ind[i] == x;)
			{
				if(col_ind[i] == y)
				{
					index = i;
					break;
				}
				else if(col_ind[i] < y)
				{
					if(inflag2 == true) return 0;
					inflag1 = true;
					i++;
				}
				else
				{
					if(inflag1 == true) return 0;
					inflag2 = true;
					i--;
				}
			}
			if(index == -1)
			{
				//std::cout << "! MatrixCOO: get(): Element at " << x << ", " << y << " not found!\n";
				return 0;
			}
			return val[index];
		}
	}

	void row_sort(double* v, int* r, int* c, int len)
	{
		if(len==1) return;

		int len1, len2;
		len1 = len/2; len2 = len - len1;

		// NOTE: Can this be done in-place, perhaps using 1 or 2 temp locations per array?

		// But for now, copy the two halves of the arrays into newly-allocated arrays
		double* v1 = new double[len1];
		double* v2 = new double[len2];
		int* r1 = new int[len1];
		int* r2 = new int[len2];
		int* c1 = new int[len1];
		int* c2 = new int[len2];

		for(int i = 0; i < len1; i++)
		{
			v1[i] = v[i];
			r1[i] = r[i];
			c1[i] = c[i];
		}
		for(int i = len1; i < len; i++)
		{
			v2[i-len1] = v[i];
			r2[i-len1] = r[i];
			c2[i-len1] = c[i];
		}

		row_sort(v1,r1,c1,len1);
		row_sort(v2,r2,c2,len2);

		// now merge the two.

		int i = 0, j = 0, k = 0;
		while(k < len)
		{
			if(i == len1)
			{
				for(int x = j; x < len2; x++)
				{
					r[k] = r2[x];
					c[k] = c2[x];
					v[k] = v2[x];
					k++;
				}
				break;
			}

			if(j == len2)
			{
				for(int x = i; x < len1; x++)
				{
					r[k] = r1[x];
					c[k] = c1[x];
					v[k] = v1[x];
					k++;
				}
				break;
			}

			if(r1[i] <= r2[j])
			{
				r[k] = r1[i];
				c[k] = c1[i];
				v[k] = v1[i];
				k++;
				i++;
			}
			else
			{
				r[k] = r2[j];
				c[k] = c2[j];
				v[k] = v2[j];
				k++;
				j++;
			}
		}

		delete [] v1; delete [] r1; delete [] c1;
		delete [] v2; delete [] r2; delete [] c2;
	}

	void col_sort()
	/* Sorts elements in each row (after row-sort) by column index. Currently uses insertion sort, which is best for sorting small arrays.
	   Call only after sorting by row! */
	{
		double** vp = new double*[nrows];
		int** cp = new int*[nrows];		// array of pointers to start of each row
		int* lens = new int[nrows];		// length of each row
		for(int i = 0; i < nrows; i++)
			lens[i] = 0;

		// now calculate number of non-zero entries in each row
		int temp=-1, k = -1, iprev = 0;
		for(int i = 0; i < nnz; i++)
		{
			if(temp != row_ind[i])			// when index moves from one row to the next
			{
				k++;
				cp[k] = col_ind + i;
				vp[k] = val + i;
				lens[k]++;
			}
			else lens[k]++;
			temp = row_ind[i];
		}
		/*std::cout << "Lengths of rows: ";
		for(int i = 0; i < nrows; i++)
			std::cout << lens[i] << " ";
		std::cout << std::endl;*/

		// insertion sort
		int tempc = 0; double tempv = 0; int j = 0;
		for(int l = 0; l < nrows; l++)			// for each row
		{
			for(int i = 1; i < lens[l]; i++)
			{
				j = i;
				while(j>0 && cp[l][j] < cp[l][j-1])
				{
					// swap cp[l][j] and cp[l][j-1], same for vp
					tempc = cp[l][j];
					cp[l][j] = cp[l][j-1];
					cp[l][j-1] = tempc;
					tempv = vp[l][j];
					vp[l][j] = vp[l][j-1];
					vp[l][j-1] = tempv;
					j--;
				}
			}
		}

		delete [] vp;
		delete [] cp;
		delete [] lens;
	}

	void sort()
	{
		std::cout << "MatrixCOO: sort(): Sorting the matrix.\n";
		row_sort(val, row_ind, col_ind, nnz);
		//std::cout << "MatrixCOO: sort(): Sorting by column\n";
		col_sort();
		sorted = true;
	}

	void print_data()
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
	}

	void mprint()
	{
		std::cout << "\n";
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < ncols; j++)
				std::cout << std::setw(WIDTH) << get(i,j);
			std::cout << std::endl;
		}
	}

	Mat multiply(Mat& x)
	{
		Mat a(x.rows(),x.cols());
		a.zeros();
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n"; We don't need our arrays to be sorted.

		for(int k = 0; k < x.cols(); k++)		// for each column of x
		{
			for(int i = 0; i < nnz; i++)
			{
				a(row_ind[i],k) += val[i]*x(col_ind[i],k);
			}
		}
		return a;
	}

	void multiply_parts(Mat* x, Mat* y, Mat* ans, int p)
	/* Like the multiply() method, except the argument matrix is considered x for the first p-1 rows, 0 in the pth row and y for the remaining rows. */
	{
		ans->zeros();
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n";

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			for(int i = 0; i < nnz; i++)
			{
				if(col_ind[i] < p)
					(*ans)(row_ind[i],k) += val[i] * x->get(col_ind[i],k);
				else if(col_ind[i] > p)
					(*ans)(row_ind[i],k) += val[i] * y->get(col_ind[i],k);
			}
		}
	}

	double getelem_multiply_parts(int rownum, Mat* x, Mat* y, int p)
	/* Returns dot product of rownum'th row of this matrix with a certain vector. This vector is composed of x for for the first p-1 elements and y for p+1'th element onwards, with zero as pth element.
	NOTE: Make sure dimensions of x and y are same. */
	{
		double ans = 0;
		//if(sorted == false) std::cout << "! MatrixCOO: multiply(): Matrix not sorted yet!\n";

		for(int k = 0; k < x->cols(); k++)		// for each column of x
		{
			for(int i = 0; i < nnz; i++)
			{
				if(row_ind[i] == rownum)
				{
					if(col_ind[i] < p)
						ans += val[i] * x->get(col_ind[i],k);
					else if(col_ind[i] > p)
						ans += val[i] * y->get(col_ind[i],k);
				}
			}
		}
	}

	void get_diagonal(Mat* D)
	{
		//if(sorted == false) { std::cout << "! MatrixCOO: get_diagonal(): Not sorted!!\n"; return; } does not depend on sorting
		D->zeros();
		for(int i = 0; i < nnz; i++)
		{
			if(row_ind[i] == col_ind[i])
			{
				(*D)(row_ind[i]) = val[i];
			}
		}
	}
};

class MatrixCRS : public SparseMatrix
{
	int* col_ind;
	int* row_ptr;
	int nnz;
	int nrows;
	int ncols;
	int valsize;		// total size allocated to array val, and therefore also to array col_ind

public:
	MatrixCRS(int num_rows, int num_cols) : SparseMatrix(num_rows, num_cols)
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
	~MatrixCRS()
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

	double get(int x, int y=0)
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

	double& operator()(int x, int y=0)
	{

	}
};

}
