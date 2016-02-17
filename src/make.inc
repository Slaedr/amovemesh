OPTIONS = -std=c++11

# if Eigen3 path is set in the env var EIGEN_DIR, append a macro definition to OPTIONS
ifdef EIGEN_DIR
OPTIONS = $(OPTIONS) -DEIGEN_LIBRARY
endif 
