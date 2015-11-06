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

curveh:
	$(GPP) curvedmeshgen2dh.cpp -o curveh.out -std=c++11 -ggdb -I$(PRMM)/include -fopenmp
	
curvee:
	$(GPP) linelast-curvedmeshgen2d.cpp -o curvee.out -std=c++11 -ggdb -I$(PRMM)/include -fopenmp
test:
	$(GPP) test.cpp -o test.out -std=c++11 -I$(PRMM)/include -ggdb
