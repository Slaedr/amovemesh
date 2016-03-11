include ../make.inc
skewmetric2d:
	$(GPP) skewmetric2d.cpp -o skewmetric2d.out -std=c++11 -ggdb -I$(PRMM)/include

shapemetric2d:
	$(GPP) shapemetric2d.cpp -o shapemetric2d.out -std=c++11 -ggdb -I$(PRMM)/include

convert:
	$(GPP) convertformat.cpp -o convert.out $(OPTIONS) -I../include

q2t:
	$(GPP) convertQuadToTri.cpp -o q2t.out -std=c++11 -ggdb -I$(PRMM)/include

quadratize:
	$(GPP) quadratize2d.cpp -o quadratize2d.out $(OPTIONS) -I$(PRMM)/include

jac3d:
	$(GPP) jac3d.cpp -o jac3d.out $(OPTIONS) -I../../include