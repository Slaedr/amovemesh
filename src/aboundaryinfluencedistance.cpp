#include <aboundaryinfluencedistance.hpp>

namespace amc{

/// Function for the adjustment factor
amc_real g(amc_real y)
{
	if(y < 0.5)
		return exp(y*y);
	else
	{
		amc_real val = exp(0.25);
		return val*y + val*0.5;
	}
}

void boundaryInfluenceDist2D(const UMesh2dh* const m, const amat::Matrix<amc_real>* const pointdisps, amat::Matrix<amc_real>* const radii)
{
	amc_int iface;
	int inofa, idim;
	amc_real l, y[NDIM2], h;
	radii->zeros();

	for(iface = 0; iface < m->gnface(); iface++)
	{
		l = (m->gcoords(m->gbface(iface,0),0)-m->gcoords(m->gbface(iface,1),0))*(m->gcoords(m->gbface(iface,0),0)-m->gcoords(m->gbface(iface,1),0));
		l += (m->gcoords(m->gbface(iface,0),1)-m->gcoords(m->gbface(iface,1),1))*(m->gcoords(m->gbface(iface,0),1)-m->gcoords(m->gbface(iface,1),1));
		l = sqrt(l);
		for(idim = 0; idim < NDIM2; idim++)
			y[idim] = pointdisps->get(m->gbface(iface,2),idim);
		h = sqrt(y[0]*y[0] + y[1]*y[1]);
		radii(iface) = l*0.5 + g(h/l);
	}
}

}
