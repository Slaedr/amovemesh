/** \file adistbinsort.hpp
 * \brief Bin sort of set of points.
 * \author Aditya Kashi
 * \date April 20, 2016
 */

#ifndef __ADISTBINSORT_H

#define __ADISTBINSORT_H 1

#ifndef __AMATRIX_H
#include <amatrix.hpp>
#endif

#ifndef _GLIBCXX_CMATH
#include <cmath>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

namespace amc {

/// Sorts a set of points into bins such that each bin contains points close to each other
template<int _ndim> class DistBinSort
{
	amc_int npoin;								///< total number of points
	const amat::Matrix<amc_real>* pointlist;	///< List of unsorted points
	amat::Matrix<amc_real>* sortedlist;			///< List of sorted points
	const int* ndbins;							///< number of bins in each direction
	int nbins;									///< Total number of bins
	std::vector<std::vector<amc_real>> rbin;	///< coordinates of bin-divisions; contains an array of numbers for each coordinate direction
	std::vector<amc_int> bmap;					///< Stores the bin-sorted index of each point in the input array
	std::vector<amc_int> invbmap;				///< For each bin-sorted point index, stores its index in the original array
	const amc_real epsilon;						///< A tolerance, used for increasing the size of the bounding cell slightly

	/// Returns the bin number given the 'position' of a bin in the coordinate directions
	int binfunc(const int* r) const;

public:
	/// Set data and compute bin sort
	/** \param[in] point_list is the initial list of unsorted points
	 * \param[in] num_bins is an array containing the number of bins to use in each coordinate direction
	 * \param[out] sorted_list is a pre-allocated list which will be filled with the sorted points
	 */
	DistBinSort(const amat::Matrix<amc_real>* const point_list, const int* num_bins, amat::Matrix<amc_real>* const sorted_list);

	inline amc_int gbmap(const amc_int index) const
	{ return bmap[index]; }
	
	inline amc_int ginvbmap(const amc_int index) const
	{ return invbmap[index]; }
};

}
#endif
