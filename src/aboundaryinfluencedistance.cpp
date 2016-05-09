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
	amc_int iface, ipoin;
	int inofa, idim;
	amc_real l, temp, disp;
	radii->zeros();

	for(iface = 0; iface < m->gnface(); iface++)
	{
		l = (m->gcoords(m->gbface(iface,0),0)-m->gcoords(m->gbface(iface,1),0))*(m->gcoords(m->gbface(iface,0),0)-m->gcoords(m->gbface(iface,1),0));
		l += (m->gcoords(m->gbface(iface,0),1)-m->gcoords(m->gbface(iface,1),1))*(m->gcoords(m->gbface(iface,0),1)-m->gcoords(m->gbface(iface,1),1));
		l = sqrt(l);
		disp = pointdisps->get(m->gbface(iface,2),0)*pointdisps->get(m->gbface(iface,2),0) + pointdisps->get(m->gbface(iface,2),1)*pointdisps->get(m->gbface(iface,2),1);
		disp = sqrt(disp);

		temp = l/2.0 + g(disp/l);
		(*radii)(m->gbface(iface,2)) = 2.0*temp;

		(*radii)(m->gbface(iface,0)) += temp;
		(*radii)(m->gbface(iface,1)) += temp;
	}
	for(ipoin = 0; ipoin < m->gnpoin(); ipoin++)
	{
		if(!m->gflag_bpoin(ipoin)) continue;
		(*radii)(ipoin) /= 2.0;
	}
}

}
