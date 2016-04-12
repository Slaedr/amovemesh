/** \file convertformatGmshDomn3d.cpp
 * \brief Converts mesh from Gmsh format to rDGFlo format.
 * \author Aditya Kashi
 *
 * The control file has the format
 * 
 * 	-input-mesh
 * 	../../output/curved-mesh-gen/bump3d-vcoarse_curved.msh
 *	-output-mesh
 *	../../output/curved-mesh-gen/bump2d-vcoarse_curved.domn
 *	-number-of-Farfield-markers
 *	1
 *	-Farfield-markers
 *	5
 *	-number-of-symmetry-markers
 *	3
 *	-Symmetry-markers
 *	32 33 4
 *	-number-of-wall-markers
 *	1
 *	-Wall-markers
 *	1
 */

#ifndef __AMESH3D_H
#include <amesh3d.hpp>
#endif

using namespace amc;
using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 2) {
		cout << "Please provide control file name." << endl;
		return -1;
	}

	string confile(argv[1]), inmeshfile, outmeshfile, dum;
	vector<int> ffn, sn, wn;
	int ffsz, ssz, wsz;

	ifstream conf(confile);
	conf >> dum; conf >> inmeshfile;
	conf >> dum; conf >> outmeshfile;
	conf >> dum; conf >> ffsz;
	conf >> dum;
	ffn.resize(ffsz);
	for(int i = 0; i < ffsz; i++)
		conf >> ffn[i];
	conf >> dum; conf >> ssz;
	conf >> dum;
	sn.resize(ssz);
	for(int i = 0; i < ssz; i++)
		conf >> sn[i];
	conf >> dum; conf >> wsz;
	conf >> dum;
	wn.resize(wsz);
	for(int i = 0; i < wsz; i++)
		conf >> wn[i];
	conf.close();

	UMesh m;
	m.readGmsh2(inmeshfile,3);
	m.writeDomn(outmeshfile,ffn,sn,wn);

	return 0;
}
