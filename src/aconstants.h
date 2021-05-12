#ifndef AMC_CONSTANTS_H
#define AMC_CONSTANTS_H
	
#define PI 3.14159265358979323846
#define SQRT3 1.73205080756887729353

/// tolerance to check if something is zero, ie, macine epsilon for double; can use: numeric_limits<floating_type>::epsilon() from the file limits
#define ZERO_TOL 2.22045e-16

/// A small number likely smaller than most convergence tolerances
#define A_SMALL_NUMBER 1e-12

#ifndef MESHDATA_DOUBLE_PRECISION
#define MESHDATA_DOUBLE_PRECISION 20
#endif

#define NDIM2 2
#define NDIM3 3

namespace amc 
{
	/// the floating point type to be used
	typedef double amc_real;
	
	/// the index (counting) type to be used
	typedef int amc_int;
}

#endif
