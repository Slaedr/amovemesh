/** @brief Class to serve as abstract generic mesh class
 * @author Aditya Kashi
 * @date Feb 24, 2016
 */

#ifndef __AMESHBASE_H
#define __AMESHBASE_H 1

namespace amc {

/// Possible mesh types
enum MeshType
{
	STRUCTURED,
	UNSTRUCTURED,
	TRIANGULAR,
	QUADRANGULAR,
	TETRAHEDRAL,
	HEXAHEDRAL,
	HYBRID
};

typedef int Integer;

/// Abstract generic mesh class
class Mesh
{
protected:
	MeshType meshtype;		/// A number describing the [kind of mesh](@ref MeshType)
	Integer ndim;			/// Dimensionality of mesh - 1D, 2D or 3D
};

} // end namespace

#endif