OPTIONS = -std=c++11 -ggdb
SUPERLU_REQ = -I$(SUPERLU_INC_DIR) -l$(SUPERLU_LIB) -L$(SUPERLU_LIB_DIR)

splinec:
	$(GPP) calc-bmotion.cpp -o cb.out -std=c++11 -I$(PRMM)/include -ggdb

crbf:
	g++-5 gencurved.cpp -o crbf.out -std=c++11 -I/media/aditya/Data/include -ggdb

convert:
	$(GPP) convertDomnGmsh2.cpp -o conv.out -std=c++11 -I$(PRMM)/include -ggdb

quadratize:
	$(GPP) quadratize2d.cpp -o quadratize2d.out -std=c++11 -ggdb -I$(PRMM)/include

curve:
	$(GPP) curvedmeshgen2d.cpp -o curve.out -std=c++11 -ggdb -I$(PRMM)/include

curvedgm:
	$(GPP) curvedmeshgen2d-dgm.cpp -o curvedgm.out -I../../include $(OPTIONS)

curveh:
	$(GPP) curvedmeshgen2dh.cpp -c $(OPTIONS) -I$(PRMM)/include -I$(SUPERLU_INC_DIR) -fopenmp
	$(GPP) curvedmeshgen2dh.o $(SUPERLU_LIB_DIR)lib$(SUPERLU_LIB).a -fopenmp -o curveh.out
	
curvee:
	$(GPP) linelast-curvedmeshgen2d.cpp -o curvee.out -std=c++11 -ggdb -I$(PRMM)/include
test:
	$(GPP) test.cpp -c $(OPTIONS) -I$(PRMM)/include -I$(SUPERLU_INC_DIR) -fopenmp
	$(GPP) test.o $(SUPERLU_LIB_DIR)lib$(SUPERLU_LIB).a -fopenmp -o test.out
