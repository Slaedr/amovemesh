AMoCurve
=========

A C++ library for mesh movement and curved mesh generation.

Several methods of mesh movement have been implemented.

Building
--------
Create a directory for binaries, eg, bin, and then run CMake for the src directory. This creates makefiles in the bin directory.
So assuming you are currently in the top level AMoCurve directory,

	mkdir bin
	cd bin
	cmake ../src -DCMAKE_BUILD_TYPE=Release
	make

will build an optimized version of the binaries. To compile a debug version with all warnings during compile time, use -DCMAKE_BUILD_TYPE=Debug instead.

Using
-----
If you do not have access to displacements for high-order ndoes, use the executable amc for 3D meshes and curveh for hybrid 2D meshes. If you have access to a quadratic mesh with boundaries already curved (but perhaps containing invalid elements), use cmg or cmg2d for 2D and 3D meshes respectively. Note that the 3D mesh class currently only supports meshes with a single type of element, while a hybrid mesh capability is available for 2D.

See test-cases/ for example control files for the programs. You *must* stick to the format specified in the example control files, until such time as more advanced control-file-parsing is incorporated.

Documentation
-------------
For documentation related to the methodology used for this project, see the docs/thesis directory. To build the code's reference documentation, issue
	doxygen builddoc.cfg
from within the docs/ directory. This will generate html documentation from source code comments.

Things to do
------------
- Create a wiki
- Write 3D hybrid mesh class
- Remove old classes such as UTriMeshCurved etc
- Test the 3D Delaunay triangulation more extensively
- Use the templated array class from amatrixt.hpp in all downstream code - this will improve performance.
- Speed up the 3D Delaunay triangulation implemenation by using more efficient data structures and by bin-sorting the points.
- Make ndim a template parameter to the mesh classes - this will result in performance gains through loop-unrolling, hopefully.
- Use YAML-cpp (or something similar) for processing control files
- In the mesh classes, extend compute_topological() to high-order elements
