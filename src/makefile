include make.inc

all: amatrix2.o alinalg.o

amatrix2.o: amatrix2.cpp amatrix2.hpp
	$(GPP) $(OPTIONS) -c -o $(PREFIX)/$@ amatrix2.cpp

amatrix2.hpp:

alinalg.o: alinalg.cpp alinalg.hpp
	$(GPP) $(OPTIONS) -c -o $(PREFIX)/$@ alinalg.cpp

alinalg.hpp:

clean:
	rm $(PREFIX)/*.o
