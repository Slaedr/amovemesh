AMoCurve
=========

C++ programs for mesh movement and curved mesh generation.

Several methods of mesh movement have been implemented, including spring analogies, linear elasticity, Delaunay graph mapping, radial basis function mapping etc.

Curved meshes can be generated, from linear meshes, for both 2D and 3D cases. Currently, 3D cases require specification of boundary displacements.

Building
--------
Requires the [Eigen](http://eigen.tuxfamily.org/) matrix library version 3.2.8 or newer. An environment variable `EIGEN_DIR` needs to be set to the top-level Eigen directory, if it is not in a standard include location. 

Create a directory for binaries, eg, build, and then run CMake for the AMoCurve/src directory. This creates makefiles in the bin directory. The compiler to be used should be set in the environment variable CXX, or can be passed as an argument to the CMake command via -DCMAKE_CXX_COMPILER=path_to_compiler. If the compiler is not GNU, you may have to edit the top-level CMakeLists.txt file to set debug and optimization flags. C++ 11 support is required. GCC 5 series is known to work.

So assuming you are currently in the top level AMoCurve directory, and a (GNU) compiler is set in CXX,

	mkdir build
	cd build
	cmake ../src -DCMAKE_BUILD_TYPE=Release
	make

will build an optimized version of the binaries. To compile a debug version with all warnings during compile time, and without aggressive optimization for a quick and light compilation, use `-DCMAKE_BUILD_TYPE=Debug' instead. Note that the release version runs significantly faster than the debug version, though. For a parallel build using, for example, 4 threads, use
	
	make -j4

instead of just `make'.

Usage
-----
If you do not have access to displacements for high-order nodes, use the executable amc for 3D meshes and curveh for hybrid 2D meshes. If you have access to a quadratic mesh with boundaries already curved (but perhaps containing invalid elements), use cmg or cmg2d for 2D and 3D meshes respectively. Note that the 3D mesh class currently only supports meshes with a single type of element, while a hybrid mesh capability is available for 2D.

See test-cases/ for example control files for the programs. You *must* stick to the format specified in the example control files, until such time as more advanced control-file-parsing is incorporated.

Documentation
-------------
For documentation related to the methodology used for this project, see the docs/thesis directory. To build the code's reference documentation, issue

	doxygen builddoc.cfg

from within the docs/ directory. This will generate html documentation from source code comments in docs/html/.

Things to do
------------
- Write 3D hybrid mesh class
- Remove old classes such as UTriMeshCurved etc
- Test the 3D Delaunay triangulation more extensively
- Speed up the 3D Delaunay triangulation implemenation by using more efficient data structures and by bin-sorting the points.
- Make ndim a template parameter to the mesh classes - this will result in performance gains through loop-unrolling, hopefully.
- Use YAML-cpp (or something similar) for processing control files
- In the mesh classes, extend compute_topological() to high-order elements
