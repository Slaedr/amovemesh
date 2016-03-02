include ../make.inc

splinec:
	$(GPP) calc-bmotion.cpp -o cb.out -std=c++11 -I$(PRMM)/include -ggdb

crbf:
	$(GPP) gencurved.cpp -o crbf.out -std=c++11 -I$(PRMM)/include -ggdb

convert:
	$(GPP) convertDomnGmsh2.cpp -o conv.out -std=c++11 -I$(PRMM)/include -ggdb

curve:
	$(GPP) curvedmeshgen2d.cpp -o curve.out -std=c++11 -ggdb -I$(PRMM)/include

curvedgm:
	$(GPP) curvedmeshgen2d-dgm.cpp -o curvedgm.out -I../../include $(OPTIONS)

curveh:
	$(GPP) $(OPTIONS) -I$(PRMM)/include -fopenmp curvedmeshgen2dh.cpp -o curveh.out
	
curvee:
	$(GPP) linelast-curvedmeshgen2d.cpp -o curvee.out -std=c++11 -ggdb -I$(PRMM)/include

