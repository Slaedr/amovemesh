#include <adistbinsort.hpp>

namespace amc {

template<int _ndim>
inline int DistBinSort<_ndim>::binfunc(const int* r) const
{
	if(_ndim == 1) return r[0];
	else if(_ndim == 2)
	{
		int tot = r[1]*ndbins[0];
		if(r[1]%2==1)
			tot += (ndbins[0]-r[0]-1);
		else
			tot += r[0];
		return tot;
	}
	else if(_ndim == 3)
	{
		int tot = r[2]*ndbins[0]*ndbins[1];
		if(r[2]%2 == 0)
		{
			tot += r[1]*ndbins[0];
			if(r[1]%2==1)
				tot += (ndbins[0]-r[0]-1);
			else
				tot += r[0];
			return tot;
		}
		else
		{
			tot += (ndbins[1]-r[1]-1);
			if(r[1]%2==1)
				tot += r[0];
			else
				tot += (ndbins[0]-r[0]-1);
			return tot;
		}
	}
	else
	{
#ifdef DEBUG
		std::cout << "DistBinSort: Dimensions greater than 3 not supported!" << std::endl;
#endif
		return 0;
	}
}

template int DistBinSort<2>::binfunc(const int* r) const;
template int DistBinSort<3>::binfunc(const int* r) const;

template<int _ndim>
DistBinSort<_ndim>::DistBinSort(const amat::Matrix<amc_real>* const point_list, const int* num_bins, amat::Matrix<amc_real>* const sorted_list) 
	: pointlist(point_list), ndbins(num_bins), sortedlist(sorted_list),  epsilon(1e-10)
{
	int idim, ibin;
	amc_int ipoin;

	npoin = pointlist->rows();
	std::cout << "DistBinSort: Dimension = " << pointlist->cols() << ", number of points = " << npoin << ", number if bins in\n";
	for(idim = 0; idim < _ndim; idim++)
		std::cout << "	Dim " << idim << " : " << ndbins[idim] << std::endl;
	
	if(pointlist->cols() != _ndim)
	{
		std::cout << "DistBinSort: ! Input data should be a list of points of dimension " << _ndim << ", not " << pointlist->cols() << "!" << std::endl;
		return;
	}

	nbins = 1;
	for(idim = 0; idim < _ndim; idim++)
		nbins *= ndbins[idim];
	std::cout << "DistBinSort: Total number of bins = " << nbins << std::endl;

	// allocate
	rbin.resize(_ndim);
	for(idim = 0; idim < _ndim; idim++)
		rbin[idim].resize(ndbins[idim]+1);
	bmap.resize(npoin);
	invbmap.resize(npoin);

	// compute bounds of the point set
	std::vector<amc_real> rmax(_ndim,0), rmin(_ndim,0);
	for(ipoin = 0; ipoin < npoin; ipoin++)
	{
		for(idim = 0; idim < _ndim; idim++)
		{
			if(pointlist->get(ipoin,idim) < rmin[idim]) rmin[idim] = pointlist->get(ipoin,idim);
			if(pointlist->get(ipoin,idim) > rmax[idim]) rmax[idim] = pointlist->get(ipoin,idim);
		}
	}
	for(idim = 0; idim < _ndim; idim++)
	{
		rmin[idim] -= epsilon;
		rmax[idim] += epsilon;
	}
	std::cout << "DistBinSort: Bounds of the point set are\n";
	for(idim = 0; idim < _ndim; idim++)
		std::cout << rmin[idim] << " ... " << rmax[idim] << '\n';
	std::cout << std::endl;

	std::vector<amc_real> binsize(_ndim);
	for(idim = 0; idim < _ndim; idim++)
	{
		binsize[idim] = (rmax[idim]-rmin[idim])/ndbins[idim];
	}

	// set bin boundaries; rbin stores the location of the starting point of a bin
	for(idim = 0; idim < _ndim; idim++)
		for(ibin = 0; ibin < ndbins[idim]+1; ibin++)
			rbin[idim][ibin] = rmin[idim] + ibin*binsize[idim];

	for(idim = 0; idim < _ndim; idim++)	
		if(fabs(rmax[idim] - rbin[idim][ndbins[idim]]) > ZERO_TOL)
			std::cout << "DistBinSort: ! Error in computing bin sizes!" << std::endl;

	std::vector<std::vector<amc_int>> ptrs(nbins);		// will store, for each point in each bin, its index in the original list
	for(ibin = 0; ibin < nbins; ibin++)
		ptrs[ibin].reserve(npoin/nbins+1);

	int binloc[_ndim], curbin, first, last;

	// start loop over points
	for(ipoin = 0; ipoin < npoin; ipoin++)
	{
		for(idim = 0; idim < _ndim; idim++)
		{
			// do a binary search to find the bin for this point
			//std::cout << "Started searching for point " << ipoin << ", in dir " << idim << std::endl;
			/*first = 0; last = ndbins[idim]-1;
			while(true)
			{
				if(first == last)
				{
					binloc[idim] = first;
					break;
				}
				curbin = (first+last)/2;
				if(pointlist->get(ipoin,idim) < rbin[idim][curbin+1])
					last = curbin;
				else
					first = curbin+1;
			}*/
			// or linear search
			ibin = 0;
			while(pointlist->get(ipoin,idim) >= rbin[idim][ibin])
			{
				ibin++;
#ifdef DEBUG
				if(ibin >= ndbins[idim])
				{
					std::cout << "DistBinSort: ! Something's wrong. Linear search couldn't locate point " << ipoin << std::endl;
					break;
				}
#endif
			}
			binloc[idim] = ibin-1;
		}

		// get global bin number from the dimensional bin numbers
		curbin = binfunc(binloc);
		std::cout << "DistBinSort: Point " << ipoin << " is located in bin " << curbin << std::endl;
		ptrs[curbin].push_back(ipoin);
	}

	amc_int i, k = 0;
	for(ibin = 0; ibin < nbins; ibin++)
	{
		for(i = 0; i < ptrs[ibin].size(); i++)
		{
			for(idim = 0; idim < _ndim; idim++)
				(*sortedlist)(k,idim) = pointlist->get(ptrs[ibin][i],idim);
			
			// update bin map and inverse bin map
			bmap[ptrs[ibin][i]] = k;
			invbmap[k] = ptrs[ibin][i];
			
			k++;
		}
	}
	std::cout << "DistBinSort: Number of points stored in the sorted list =  " << k << std::endl;
}

template DistBinSort<2>::DistBinSort(const amat::Matrix<amc_real>* const point_list, const int* num_bins, amat::Matrix<amc_real>* const sorted_list);
template DistBinSort<3>::DistBinSort(const amat::Matrix<amc_real>* const point_list, const int* num_bins, amat::Matrix<amc_real>* const sorted_list);

} // end namespace
