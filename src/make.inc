# set generic options, used by all programs, here
OPTIONSS = -std=c++11 -ggdb
PREFIX = $(PRMM)/bin

# if Eigen3 path is set in the env var EIGEN_DIR, append a macro definition and the include path to OPTIONS
ifdef EIGEN_DIR
OPTIONS = $(OPTIONSS) -DEIGEN_LIBRARY -I$(EIGEN_DIR)
else
OPTIONS = $(OPTIONSS)
endif 
