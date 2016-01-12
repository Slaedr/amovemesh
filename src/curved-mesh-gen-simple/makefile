OPTIONS = -std=c++11 -ggdb

curvedrbf:
	g++-5 gencurved-rbf.cpp -o crbf.out -std=c++11 -I/media/aditya/Data/include -ggdb

curveddgm:
	$(GPP) gencurved-dgm.cpp -o cdgm.out -I../../include $(OPTIONS)

curveddgrbf:
	g++-5 gencurved-dgrbf.cpp -o cdgrbf.out -std=c++11 -I/media/aditya/Data/include -ggdb

compare:
	g++-5 compare.cpp -o comp.out -std=c++11 -I/media/aditya/Data/include -ggdb

crbf3d:
	g++-5 gencurved-rbf-3d.cpp -o crbf3d.out -std=c++11 -I/media/aditya/Data/include -ggdb -fopenmp

getq:
	g++-5 getquadmesh.cpp -o getq.out -std=c++11 -I/media/aditya/Data/include -ggdb	
