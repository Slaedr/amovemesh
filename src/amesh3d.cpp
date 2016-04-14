/// @brief Implementation of data structure and setup for 3D unstructured mesh.
/// @author Aditya Kashi
/// @date August 20, 2015

#include "amesh3d.hpp"

namespace amc {

/** No-arg constructor. */
UMesh::UMesh() {alloc_jacobians = false; psup = nullptr; elsed = nullptr;}

UMesh::UMesh(const UMesh& other)
{
	npoin = other.npoin;
	nelem = other.nelem;
	nface = other.nface;
	ndim = other.ndim;
	nnode = other.nnode;
	naface = other.naface;
	nbface = other.nbface;
	nfael = other.nfael;
	nnofa = other.nnofa;
	nedel = other.nedel;
	nedfa = other.nedfa;
	nnoded = other.nnoded;
	nedge = other.nedge;
	nbedge = other.nbedge;
	nbtag = other.nbtag;
	ndtag = other.ndtag;
	coords = other.coords;
	inpoel = other.inpoel;
	bface = other.bface;
	flag_bpoin = other.flag_bpoin;
	vol_regions = other.vol_regions;
	esup = other.esup;
	esup_p = other.esup_p;
	if(other.psup != nullptr)
	{
		psup = new std::vector<int>[npoin];
		for(int i = 0; i < npoin; i++)
			psup[i] = other.psup[i];
	}
	if(other.elsed != nullptr)
	{
		elsed = new std::vector<int>[nedge];
		for(int i = 0; i < nedge; i++)
			elsed[i] = other.elsed[i];
	}
	esuel = other.esuel;
	edgepo = other.edgepo;
	elsed = other.elsed;
	intfac = other.intfac;
	alloc_jacobians = other.alloc_jacobians;
	jacobians = other.jacobians;
}

UMesh::~UMesh()
{
	if(psup != nullptr)
		delete [] psup;
	if(elsed != nullptr)
		delete [] elsed;
}

// Reads mesh from Gmsh 2 format file. For quadratic meshes, mapping has to be applied for node-ordering.
void UMesh::readGmsh2(std::string mfile, int dimensions)
{
	std::cout << "UMesh3d: readGmsh2(): Reading mesh file...\n";
	int dum; double dummy; std::string dums; char ch;
	ndim = dimensions;

	std::ifstream infile(mfile);
	for(int i = 0; i < 4; i++)		//skip 4 lines
		do
			ch = infile.get();
		while(ch != '\n');

	infile >> npoin;
	std::cout << "UMesh3d: readGmsh2(): No. of points = " << npoin << std::endl;
	coords.setup(npoin,ndim);

	// read coords of points
	for(int i = 0; i < npoin; i++)
	{
		infile >> dum;
		for(int j = 0; j < ndim; j++)
			infile >> coords(i,j);
		if(ndim < 3) infile >> dummy;
	}
	infile >> dums;
	infile >> dums;
	//std::cout << "UMesh3d: readGmsh2(): coords read." << std::endl;

	int elmtype, nbtags, ntags, nskipped= 0;
	amc_int nelm, i;
	ndtag = 0; nbtag = 0;
	infile >> nelm;
	amat::Matrix<amc_int> elms(nelm,40);
	nface = 0; nelem = 0;
	std::cout << "UMesh3d: readGmsh2(): Total number of elms is " << nelm << std::endl;
	amat::Matrix<amc_int> temper(27,1);		// to temporarily store nodes in Gmsh order

	for(i = 0; i < nelm; i++)
	{
		infile >> dum;
		infile >> elmtype;
		// assumption here is that elm-type is same for all edges, same for all boundary faces and same for all elements
		switch(elmtype)
		{
			case(2): // linear triangle face
				nnofa = 3;
				nnoded = 2;
				nedfa = 3;
				infile >> nbtags;
				if(nbtags > nbtag) nbtag = nbtags;
				for(int j = 0; j < nbtags; j++)
					infile >> elms(i,j+nnofa);		// get tags
				for(int j = 0; j < nnofa; j++)
					infile >> elms(i,j);			// get node numbers
				nface++;
				break;
			case(3): // linear quad (4-node) face
				nnofa = 4;
				nnoded = 2;
				nedfa = 4;
				infile >> nbtags;
				if(nbtags > nbtag) nbtag = nbtags;
				for(int j = 0; j < nbtags; j++)
					infile >> elms(i,j+nnofa);		// get tags
				for(int j = 0; j < nnofa; j++)
					infile >> elms(i,j);			// get node numbers
				nface++;
				break;
			case(9): // quadratic triangle face
				nnofa = 6;
				nnoded = 3;
				nedfa = 3;
				infile >> nbtags;
				if(nbtags > nbtag) nbtag = nbtags;
				for(int j = 0; j < nbtags; j++)
					infile >> elms(i,j+nnofa);		// get tags
				for(int j = 0; j < nnofa; j++)
					infile >> elms(i,j);			// get node numbers
				nface++;
				break;
			case(10): // quadratic quad (9-node) face
				nnofa = 9;
				nnoded = 3;
				nedfa = 4;
				infile >> nbtags;
				if(nbtags > nbtag) nbtag = nbtags;
				for(int j = 0; j < nbtags; j++)
					infile >> elms(i,j+nnofa);		// get tags
				for(int j = 0; j < nnofa; j++)
					infile >> elms(i,j);			// get node numbers
				nface++;
				break;
			case(4): // linear tet
				nnode = 4;
				nfael = 4;
				nnofa = 3;
				nnoded = 2;
				nedel = 6;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnode);		// get tags
				for(int j = 0; j < nnode; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(5):	// linear hex
				nnode = 8;
				nfael = 6;
				nnofa = 4;
				nnoded = 2;
				nedel = 12;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnode);		// get tags
				for(int j = 0; j < nnode; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(11):	// quadratic tet
				nnode = 10;
				nfael = 4;
				nnofa = 6;
				nnoded = 3;
				nedel = 6;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnode);		// get tags
				for(int j = 0; j < nnode; j++)
					infile >> elms(i,j);			// get node numbers
				nelem++;
				break;
			case(12):	// quadratic hex (27 nodes)
				nnode = 27;
				nfael = 6;
				nnofa = 9;
				nnoded = 3;
				nedel = 12;
				infile >> ntags;
				if(ntags > ndtag) ndtag = ntags;
				for(int j = 0; j < ntags; j++)
					infile >> elms(i,j+nnode);		// get tags
				for(int j = 0; j < nnode; j++)
					infile >> temper(j);			// get node numbers
				
				for(int j = 0; j <= 8; j++)		// first 8 nodes are same
					elms(i,j) = temper(j);
				
				// add mapping from gmsh format to RDGFlo format
				elms(i,9) = temper(11);
				elms(i,10) = temper(13);
				elms(i,11) = temper(9);
				elms(i,12) = temper(10);
				elms(i,13) = temper(12);
				elms(i,14) = temper(14);
				elms(i,15) = temper(15);
				elms(i,16) = temper(16);
				elms(i,17) = temper(18);
				elms(i,18) = temper(19);
				elms(i,19) = temper(17);
				elms(i,20) = temper(20);
				elms(i,21) = temper(21);
				elms(i,22) = temper(23);
				elms(i,23) = temper(24);
				elms(i,24) = temper(22);
				elms(i,25) = temper(25);
				elms(i,26) = temper(26);
				
				nelem++;
				break;
				
			default:
				std::cout << "! UMesh3d: readGmsh2(): Element type not recognized. Skipping." << std::endl;
				do
					ch = infile.get();
				while(ch != '\n');
				nskipped++;
				break;
		}
	}
	//std::cout << "UMesh3d: readGmsh2(): Done reading elms" << std::endl;

	if(nface > 0)
		bface.setup(nface, nnofa+nbtag);
	else std::cout << "UMesh3d: readGmsh2(): NOTE: There is no data to populate bface!" << std::endl;

	inpoel.setup(nelem, nnode);
	vol_regions.setup(nelem, ndtag);

	// write into inpoel and bface
	// the first nface rows to be read are boundary faces
	for(i = 0; i < nface; i++)
	{
		for(int j = 0; j < nnofa; j++)
			bface(i,j) = elms(i+nskipped,j)-1;			// -1 to correct for the fact that our numbering starts from zero
		for(int j = nnofa; j < nnofa+nbtag; j++)
			bface(i,j) = elms(i+nskipped,j);
	}
	for(i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnode; j++)
			inpoel(i,j) = elms(i+nface+nskipped,j)-1;
		for(int j = 0; j < ndtag; j++)
			vol_regions(i,j) = elms(i+nface+nskipped,j+nnode);
	}
	infile.close();

	std::cout << "UMesh3d: readGmsh2(): Setting flag_bpoin..." << std::endl;

	// set flag_bpoin
	flag_bpoin.setup(npoin,1);
	flag_bpoin.zeros();
	for(int i = 0; i < nface; i++)
		for(int j = 0; j < nnofa; j++)
			flag_bpoin(bface(i,j)) = 1;

	std::cout << "UMesh3d: readGmsh2(): Done. No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ",\n number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << ", number of nodes per edge " << nnoded << "." << std::endl;
}

void UMesh::readDomn(std::string mfile)
{
	std::ifstream fin(mfile);
	if(!fin) {
		std::cout << "UMesh: readDomn(): could not open file " << mfile << "!" << std::endl;
		return;
	}

	std::string line;
	char dumc = 'a';
	int i,j,idim, dumi;

	std::getline(fin, line);
	fin >> nelem >> npoin >> nnode >> nedel >> nfael >> nnofa >> nedfa >> nnoded;
	std::cout << "Started reading" << std::endl;
	while(dumc != '\n')
		dumc = fin.get();
	dumc = 'a';
	std::getline(fin, line);

	ndim = 3;
	std::cout << "UMesh: readDomn(): nelem = " << nelem << ", npoin = " << npoin << ", nnode = " << nnode << ", nedel = " << nedel << ", nfael = " << nfael << ", nnofa = " << nnofa << std::endl;

	coords.setup(npoin,ndim);

	if(nnode == 0)
	{
		// handle hybrid mesh
		std::cout << "Error!" << std::endl;
		return;
	}
	else
	{
		std::cout << "Reading inpoel" << std::endl;
		inpoel.setup(nelem,nnode);
		for(i = 0; i < nelem; i++)
		{
			fin >> dumi;
			for(j = 0; j < nnode; j++)
			{
				fin >> inpoel(i,j);
				inpoel(i,j) -= 1;
			}
			std::getline(fin,line);
		}

		std::getline(fin,line);

		for(i = 0; i < npoin; i++)
		{
			fin >> dumi;
			for(j = 0; j < ndim; j++)
				fin >> coords(i,j);
		}
	}
	fin.close();

	nbtag = 2; ndtag = 2;
	vol_regions.setup(nelem,ndtag);
	vol_regions.zeros();

	// Get elements surrounding points
	esup_p.setup(npoin+1,1);
	esup_p.zeros();

	amat::Matrix<amc_int> lpoin(npoin,1);
	lpoin.zeros();
	
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnode; j++)
		{
			esup_p(inpoel(i,j)+1,0) += 1;	// inpoel(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
		}
	}

	// Now make the members of esup_p cumulative
	for(int i = 1; i < npoin+1; i++)
		esup_p(i,0) += esup_p(i-1,0);
	// Now populate esup
	esup.setup(esup_p(npoin,0),1);
	esup.zeros();
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnode; j++)
		{
			int ipoin = inpoel(i,j);
			esup(esup_p(ipoin,0),0) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
			esup_p(ipoin,0) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
		}
	}
	
	//But now esup_p holds increased values - each member increased by the number elements surrounding the corresponding point.
	// So correct this.
	for(int i = npoin; i >= 1; i--)
		esup_p(i,0) = esup_p(i-1,0);
	esup_p(0,0) = 0;

	// Elements surrounding points is now done.
	
	//  Elements surrounding elements
	std::cout << "UMesh3d: readDomn(): Elements surrounding elements...\n";

	esuel.setup(nelem, nfael);
	for(int ii = 0; ii < nelem; ii++)
		for(int jj = 0; jj < nfael; jj++)
			esuel(ii,jj) = -1;

	lpofa.setup(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in [0:nnofa], i in [0:nfael])

	if(nnode == 4)								// if tet
		for(int i = 0; i < nfael; i++)
		{
			for(int j = 0; j < nnofa; j++)
			{
				//lpofa(i,j) = perm(0,nnode-1,i,j);
				lpofa(i,j) = (i+j)%nnode;
			}
		}

	if(nnode == 8)								// if hex
	{
		lpofa(0,0) = 0; lpofa(0,1) = 3; lpofa(0,2) = 2; lpofa(0,3) = 1;
		lpofa(1,0) = 0; lpofa(1,1) = 1; lpofa(1,2) = 5; lpofa(1,3) = 4;
		lpofa(2,0) = 0; lpofa(2,1) = 4; lpofa(2,2) = 7; lpofa(2,3) = 3;
		lpofa(3,0) = 1; lpofa(3,1) = 2; lpofa(3,2) = 6; lpofa(3,3) = 5;
		lpofa(4,0) = 2; lpofa(4,1) = 3; lpofa(4,2) = 7; lpofa(4,3) = 6;
		lpofa(5,0) = 4; lpofa(5,1) = 5; lpofa(5,2) = 6; lpofa(5,3) = 7;
	}

	amat::Matrix<int> lhelp(nnofa,1);
	lhelp.zeros();
	lpoin.zeros();

	for(int ielem = 0; ielem < nelem; ielem++)
	{
		for(int ifael = 0; ifael < nfael; ifael++)
		{
			for(int i = 0; i < nnofa; i++)
			{
				lhelp(i) = inpoel(ielem, lpofa(ifael,i));		// lhelp stores global node nos. of current face of current element
				lpoin(lhelp(i)) = 1;
			}
			int ipoin = lhelp(0);								// ipoin is the global node number of first point of this face
			for(int istor = esup_p(ipoin); istor < esup_p(ipoin+1); istor++)
			{
				int jelem = esup(istor);
				if(jelem != ielem)
				{
					for(int jfael = 0; jfael < nfael; jfael++)
					{
						//Assume that no. of nodes in face ifael is same as that in face jfael
						int icoun = 0;
						for(int jnofa = 0; jnofa < nnofa; jnofa++)
						{
							int jpoin = inpoel(jelem, lpofa(jfael,jnofa));
							if(lpoin(jpoin)==1) icoun++;
						}
						if(icoun == nnofa)		// if all nnofa points of face ifael are found in face jfael, they both are the same face
						{
							esuel(ielem,ifael) = jelem;
							esuel(jelem,jfael) = ielem;
						}
					}
				}
			}
			for(int i = 0; i < nnofa; i++)
				lpoin(lhelp(i)) = 0;
		}
	}

	/** Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac.
	 * The orientation of the face is such that the face points towards the element with larger index.
	 * NOTE: After the following portion, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.
	 * \sa compute_topological
	 */

	std::cout << "UMesh2d: readDomn(): Computing intfac..." << std::endl;
	nbface = naface = 0;

	// first run: calculate nbface
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int je = esuel(ie,in);
			if(je == -1)
				nbface++;
		}
	}
	std::cout << "UMesh2d: compute_topological(): Number of boundary faces = " << nbface << std::endl;
	// calculate number of internal faces
	naface = nbface;
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int je = esuel(ie,in);
			if(je > ie && je < nelem)
				naface++;
		}
	}
	std::cout << "UMesh2d: compute_topological(): Number of all faces = " << naface << std::endl;

	//allocate intfac
	intfac.setup(naface,nnofa+2);

	//reset face totals
	nbface = naface = 0;

	//second run: populate intfac
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int* inp = new int[nnofa];
			for(int inofa = 0; inofa < nnofa; inofa++)
				inp[inofa] = lpofa.get(in,inofa);

			int je = esuel.get(ie,in);
			if(je == -1)
			{
				esuel(ie,in) = nelem+nbface;
				intfac(nbface,0) = ie;
				intfac(nbface,1) = nelem+nbface;
				for(int j = 0; j < nnofa; j++)
				{
					intfac(nbface,2+j) = inpoel(ie,inp[j]);
				}

				nbface++;
			}

			delete [] inp;
		}
	}

	naface = nbface;
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int* inp = new int[nnofa];
			for(int inofa = 0; inofa < nnofa; inofa++)
				inp[inofa] = lpofa.get(in,inofa);

			int je = esuel(ie,in);
			if(je > ie && je < nelem)
			{
				intfac(naface,0) = ie;
				intfac(naface,1) = je;
				for(int j = 0; j < nnofa; j++)
					intfac(naface,2+j) = inpoel(ie,inp[j]);
				naface++;
			}

			delete [] inp;
		}
	}

	flag_bpoin.setup(nbface,1);
	flag_bpoin.zeros();

	// compute bface and set flag_bpoin using intfac
	nface = nbface;
	bface.setup(nface, nnofa+nbtag);
	for(i = 0; i < nbface; i++)
	{
		for(j = 0; j < nnofa; j++)
		{
			bface(i,j) = intfac(i,2+j);
			flag_bpoin(bface.get(i,j)) = 1;
		}
		for(j = nnofa; j < nnofa+nbtag; j++)
			bface(i,j) = 0;
	}
}

void UMesh::printmeshstats()
{
	std::cout << "UMesh3d: No. of points: " << npoin << ", number of elements: " << nelem << ", number of boundary faces " << nface << ", number of nodes per element: " << nnode << ", number of nodes per face: " << nnofa << ", number of faces per element: " << nfael << std::endl;
}

/// Changes node ordering. Use only for quadratic hexahedral mesh!!
/** Assuming inpoel contains rDGFlo format initially, m_inpoel contains Gmsh ordering on return
 */
void UMesh::mapinpoelRDGFloToGmsh()
{
	int temp;

	for(int ielem = 0; ielem < nelem; ielem++)
	{
		for(int inode = 0; inode <= 8; inode++)
			m_inpoel(ielem,inode) = inpoel(ielem,inode);

		m_inpoel(ielem,11) = inpoel(ielem,9);
		m_inpoel(ielem,13) = inpoel(ielem,10);
		m_inpoel(ielem,9) = inpoel(ielem,11);
		m_inpoel(ielem,10) = inpoel(ielem,12);
		m_inpoel(ielem,12) = inpoel(ielem,13);
		m_inpoel(ielem,14) = inpoel(ielem,14);
		m_inpoel(ielem,15) = inpoel(ielem,15);
		m_inpoel(ielem,16) = inpoel(ielem,16);
		m_inpoel(ielem,18) = inpoel(ielem,17);
		m_inpoel(ielem,19) = inpoel(ielem,18);
		m_inpoel(ielem,17) = inpoel(ielem,19);
		m_inpoel(ielem,20) = inpoel(ielem,20);
		m_inpoel(ielem,21) = inpoel(ielem,21);
		m_inpoel(ielem,23) = inpoel(ielem,22);
		m_inpoel(ielem,24) = inpoel(ielem,23);
		m_inpoel(ielem,22) = inpoel(ielem,24);
		m_inpoel(ielem,25) = inpoel(ielem,25);
		m_inpoel(ielem,26) = inpoel(ielem,26);
	}
}

// Changes node ordering from Gmsh to rDGFlo (inverse of mapinpoelRDGFloToGmsh ). Use only for quadratic hexahedral mesh!!
/* Assuming inpoel contains Gmsh format initially, m_inpoel contains rDGFlo ordering on return
 */
/*void UMesh::mapinpoelGmshToRDGFlo()
{
	int temp;

	for(int ielem = 0; ielem < nelem; ielem++)
	{
		for(int inode = 0; inode <= 8; inode++)
			m_inpoel(ielem,inode) = inpoel(ielem,inode);

		m_inpoel(ielem,9); = inpoel(ielem,11);
		m_inpoel(ielem,10) = inpoel(ielem,13);        
		m_inpoel(ielem,11) = inpoel(ielem,9) ;        
		m_inpoel(ielem,12) = inpoel(ielem,10);        
		m_inpoel(ielem,13) = inpoel(ielem,12);        
		m_inpoel(ielem,14) = inpoel(ielem,14);        
		m_inpoel(ielem,15) = inpoel(ielem,15);        
		m_inpoel(ielem,16) = inpoel(ielem,16);        
		m_inpoel(ielem,17) = inpoel(ielem,18);        
		m_inpoel(ielem,18) = inpoel(ielem,19);        
		m_inpoel(ielem,19) = inpoel(ielem,17);        
		m_inpoel(ielem,20) = inpoel(ielem,20);        
		m_inpoel(ielem,21) = inpoel(ielem,21);        
		m_inpoel(ielem,22) = inpoel(ielem,23);        
		m_inpoel(ielem,23) = inpoel(ielem,24);        
		m_inpoel(ielem,24) = inpoel(ielem,22);        
		m_inpoel(ielem,25) = inpoel(ielem,25);        
		m_inpoel(ielem,26) = inpoel(ielem,26);        
	}
}*/
	
/*void UMesh::computeIntfacNumberFromBfaceNumber()
{
	intfacFromBface.resize(nface);
	amc_int iface, jface;
	int inofa, jnofa;

	for(iface = 0; iface < nface; iface++)
	{
		for(jface = 0; jface < nbface; jface++)
		{
			std::vector<bool> foundi(nnofa,false);

			for(inofa = 0; inofa < nnofa; inofa++)
				for(jnofa = 0; jnofa < nnofa; jnofa++)
					if(intfac.get(jface,jnofa+2) == bface.get(iface,inofa))
						foundi[inofa] = true;

			bool found = true;
			for(inofa = 0; inofa < nnofa; inofa++)
				if(foundi[inofa] == false)
					found = false;

			if(found)
				intfacFromBface[iface] = jface;
		}
	}
}*/

void UMesh::findLocalFaceLocalPointConnectivityLinearElements(amat::Matrix<int>& lpofal)
{
	if(nnode == 27 || nnode == 8)
	{
		lpofal(0,0) = 0; lpofal(1,0) = 0;
		lpofal(0,1) = 1; lpofal(1,1) = 1;
		lpofal(0,2) = 2; lpofal(1,2) = 5;
		lpofal(0,3) = 3; lpofal(1,3) = 4;

		lpofal(2,0) = 0; lpofal(3,0) = 6;
		lpofal(2,1) = 4; lpofal(3,1) = 5;
		lpofal(2,2) = 7; lpofal(3,2) = 1;
		lpofal(2,3) = 3; lpofal(3,3) = 2;

		lpofal(4,0) = 6; lpofal(5,0) = 6;
		lpofal(4,1) = 2; lpofal(5,1) = 7;
		lpofal(4,2) = 3; lpofal(5,2) = 4;
		lpofal(4,3) = 7; lpofal(5,3) = 5;
	}
	else if(nnode == 10 || nnode == 4)
	{
		for(int ifael = 0; ifael < 4; ifael++)
			for(int i = 0; i < 3; i++)
				lpofal(ifael,i) = (ifael+i+1)%4;
	}
}

void UMesh::findBfaceHostCell()
{
	bfaceHostCell.setup(nface,1);
	amc_int iface, ielem;
	int inofa, inode, ifael, ncoun, i;

	amat::Matrix<int> lpofal;
	if(nnode == 8 || nnode == 27)
		lpofal.setup(nfael,4);
	else if(nnode == 4 || nnode == 10)
		lpofal.setup(nfael,3);
	findLocalFaceLocalPointConnectivityLinearElements(lpofal);

	for(iface = 0; iface < nface; iface++)
	{
		//std::cout << " " << iface << std::flush;
		bool set = false;
		for(ielem = 0; ielem < nelem; ielem++)
		{
			for(ifael = 0; ifael < nfael; ifael++)
			{
				ncoun = 0;
				// now for each bface node, check if it matches any node of this fael of this element
				for(inofa = 0; inofa < nnofa; inofa++)
				{
					if(nnode == 8 || nnode == 27)
						for(i = 0; i < 4; i++)
						{
							if(bface.get(iface,inofa) == inpoel.get(ielem, lpofal.get(ifael,i)) )
							{
								ncoun++;
								break;
							}
						}
					else if(nnode == 4 || nnode == 10)
						for(i = 0; i < 3; i++)
						{
							if(bface.get(iface,inofa) == inpoel.get(ielem, lpofal.get(ifael,i)) )
							{
								ncoun++;
								break;
							}
						}
				}
				if((nnode == 8 || nnode == 27) && ncoun == 4)
				{
					bfaceHostCell(iface) = ielem;
					set = true;
					break;
				}
				if((nnode == 4 || nnode == 10) && ncoun == 3)
				{
					bfaceHostCell(iface) = ielem;
					set = true;
					break;
				}
			}

			// if element has been found for this face, break from ielem loop
			if(set) break;
		}
	}
	//std::cout << std::endl;
}

/* Writes mesh to Gmsh2 file format. */
void UMesh::writeGmsh2(std::string mfile)
{
	std::cout << "UMesh2d: writeGmsh2(): writing mesh to file " << mfile << std::endl;

	m_inpoel.setup(nelem,nnode);
	
	if(nnode == 27)
		mapinpoelRDGFloToGmsh();
	else
		m_inpoel = inpoel;

	// decide element type first, based on nfael/nnode and nnofa
	int elm_type = 4;
	int face_type = 2;
	if(nnode == 8)
		elm_type = 5;
	else if(nnode == 10)
		elm_type = 11;
	else if(nnode == 27)
		elm_type = 12;

	if(nnofa == 4) face_type = 3;
	else if(nnofa == 6) face_type = 9;
	else if(nnofa == 9) face_type = 10;

	std::ofstream outf(mfile);
	outf << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
	outf << "$Nodes\n" << npoin << '\n';
	outf << std::setprecision(MESHDATA_DOUBLE_PRECISION);
	for(int ip = 0; ip < npoin; ip++)
	{
		outf << ip+1;
		for(int j = 0; j < ndim; j++)
			outf << " " << coords.get(ip,j);
		for(int j = 3-ndim; j > 0; j--)
			outf << " " << 0.0;
		outf << '\n';
	}
	outf << "$EndNodes\n";

	// "elements"
	outf << "$Elements\n" << nelem+nface << '\n';
	// boundary faces first
	for(int iface = 0; iface < nface; iface++)
	{
		outf << iface+1 << " " << face_type << " " << nbtag;
		for(int i = nnofa; i < nnofa+nbtag; i++)
			outf << " " << bface.get(iface,i);			// write tags
		for(int i = 0; i < nnofa; i++)
			outf << " " << bface.get(iface,i)+1;		// write nodes
		outf << '\n';
	}

	// std::cout << "elements\n";
	for(int iel = 0; iel < nelem; iel++)
	{
		outf << nface+iel+1 << " " << elm_type << " " << ndtag;
		for(int i = 0; i < ndtag; i++)
			outf << " " << vol_regions.get(iel,i);
		for(int i = 0; i < nnode; i++)
			outf << " " << m_inpoel.get(iel,i)+1;
		outf << '\n';
	}
	outf << "$EndElements\n";

	outf.close();
}

void UMesh::writeDomn(std::string mfile, std::vector<int> farfieldnumber, std::vector<int> symmetrynumber, std::vector<int> wallnumber)
{
	std::cout << "UMesh: writeDomn(): Writing domn file to " << mfile << std::endl;
	m_inpoel.setup(nelem,nnode);
	
	std::ofstream fout(mfile);
	fout << "npoin ntetr npyra npris nhexa ntria nquad time\n";
	if(nnode == 8 || nnode == 27)
		fout << npoin << " 0 0 0 " << nelem << " 0 " << nface << " 0.0\n";
	else if(nnode == 6 || nnode == 10)
		fout << npoin << " " << nelem << " 0 0 0 " << nface << " 0.0\n";

	fout << "element connectivity\n";
	for(int ielem = 0; ielem < nelem; ielem++)
	{
		fout << ielem+1;
		for(int inode = 0; inode < nnode; inode++)
			fout << ' ' << inpoel.get(ielem,inode)+1;
		fout << '\n';
	}

	std::cout << "UMesh: writeDomn: finding bface host cells..." << std::flush;
	findBfaceHostCell();
	std::cout << "Done." << std::endl;

	fout << "face connectivity\n";
	int marker, ffmarker;
	for(int iface = 0; iface < nface; iface++)
	{
		for(int inum = 0; inum < farfieldnumber.size(); inum++)	
			if(bface.get(iface,nnofa) == farfieldnumber[inum])
			{
				marker = 4;
				ffmarker = 1;
			}
		for(int inum = 0; inum < symmetrynumber.size(); inum++)
			if(bface.get(iface,nnofa) == symmetrynumber[inum])
			{
				marker = 3;
				ffmarker = 0;
			}
		for(int inum = 0; inum < wallnumber.size(); inum++)
			if(bface.get(iface,nnofa) == wallnumber[inum])
			{
				marker = 0;
				ffmarker = 0;
			}

		fout << iface+1 << " " << marker << " " << ffmarker;
		fout << " 0 0 " << bfaceHostCell.get(iface)+1 << " ";
		for(int inofa = 0; inofa < nnofa; inofa++)
			fout << " " << bface.get(iface,inofa)+1;
		// we need nine numbers, so if nnofa < 9, we fill in dummy zeros
		for(int i = 9-nnofa; i > 0; i--)
			fout << " 0";
		fout << '\n';
	}

	fout << "coordinates of nodes\n";
	fout << std::setprecision(MESHDATA_DOUBLE_PRECISION);
	for(int ipoin = 0; ipoin < npoin; ipoin++)
	{
		fout << ipoin+1;
		for(int idim = 0; idim < ndim; idim++)
			fout << " " << coords.get(ipoin,idim);
		fout << '\n';
	}
	fout.close();
}

/// Computes jacobians for linear elements
/** Currently only for tetrahedra
 */
void UMesh::compute_jacobians()
{
	if(!alloc_jacobians)
	{
		jacobians.setup(nelem,1);
		alloc_jacobians = true;
	}
	jacobians.zeros();
	
	if(nnode == 4)
	{
		amc_int ielem; int inode;
		amc_real x1, y1, z1, x21, y21, z21, x31, y31, z31, x41, y41, z41;
		
		for(ielem = 0; ielem < nelem; ielem++)
		{
			x1 = coords(inpoel(ielem,0),0); y1 = coords(inpoel(ielem,0),1); z1 = coords(inpoel(ielem,0),2);

			x21 = coords(inpoel(ielem,1),0) - x1;
			y21 = coords(inpoel(ielem,1),1) - y1;
			z21 = coords(inpoel(ielem,1),2) - z1;

			x31 = coords(inpoel(ielem,2),0) - x1;
			y31 = coords(inpoel(ielem,2),1) - y1;
			z31 = coords(inpoel(ielem,2),2) - z1;

			x41 = coords(inpoel(ielem,3),0) - x1;
			y41 = coords(inpoel(ielem,3),1) - y1;
			z41 = coords(inpoel(ielem,3),2) - z1;

			jacobians(ielem) = x21*(y31*z41-z31*y41) + x31*(z21*y41-y21*z41) + x41*(y21*z31-z21*y31);
		}
	}
	else
	{
		std::cout << "UMesh: compute_jacobians(): ! Not implemented for this mesh type!!" << std::endl;
	}
}

/** \brief Computes various connectivity data structures for the mesh.
 *
 * These include
 * - Elements surrounding points (esup and esup_p)
 * - Points surrounding points (psup)
 * - Elements surrounding elements (esuel)
 * - Elements surrounding edge (elsed)
 * - Edge data structure (edgepo)
 * - Face data structure (intfac)
 * 
 * \note NOTE: Currently only works for linear mesh - and psup works only for tetrahedral or hexahedral linear mesh
 */
void UMesh::compute_topological()
{
	std::cout << "UMesh: compute_topological(): Calculating and storing topological information..." << std::endl;
	amc_int ied, i, iel, jel;
	
	//1. Elements surrounding points
	esup_p.setup(npoin+1,1);
	esup_p.zeros();
	
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnode; j++)
		{
			esup_p(inpoel(i,j)+1,0) += 1;	// inpoel(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
		}
	}

	// Now make the members of esup_p cumulative
	for(int i = 1; i < npoin+1; i++)
		esup_p(i,0) += esup_p(i-1,0);
	// Now populate esup
	esup.setup(esup_p(npoin,0),1);
	esup.zeros();
	for(int i = 0; i < nelem; i++)
	{
		for(int j = 0; j < nnode; j++)
		{
			int ipoin = inpoel(i,j);
			esup(esup_p(ipoin,0),0) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
			esup_p(ipoin,0) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
		}
	}
	
	//But now esup_p holds increased values - each member increased by the number elements surrounding the corresponding point.
	// So correct this.
	for(int i = npoin; i >= 1; i--)
		esup_p(i,0) = esup_p(i-1,0);
	esup_p(0,0) = 0;

	// Elements surrounding points is now done.

	//2. Points surrounding points - works only for tets and hexes!!
	std::cout << "UMesh2d: compute_topological(): Points surrounding points\n";
	if(psup != nullptr)
		delete [] psup;
	psup = new std::vector<int>[npoin];
	for(int i = 0; i < npoin; i++)
		psup[i].reserve(10);

	amat::Matrix<int> lpoin(npoin,1);  // The ith member indicates the global point number of which the ith point is a surrounding point

	for(int i = 0; i < npoin; i++)
		lpoin(i,0) = -1;	// initialize this std::vector to -1

	for(int ip = 0; ip < npoin; ip++)
	{
		lpoin(ip) = ip;
		// iterate over elements surrounding point
		for(int i = esup_p(ip); i < esup_p(ip+1); i++)
		{
			int ielem = esup(i);
			int inode;

			// find local node number of ip in ielem -- needed for anything except tetrahedral mesh
			for(int jnode = 0; jnode < nnode; jnode++)
				if(inpoel(ielem,jnode) == ip) inode = jnode;

			std::vector<bool> nbd(nnode);		// contains true if that local node number is connected to inode
			for(int j = 0; j < nnode; j++)
				nbd[j] = false;

			if(nnode == 4)					// for a tet, all local nodes are connected to any given node
			{
				for(int inb = 0; inb < nbd.size(); inb++)
					nbd[inb] = true;
			}
			if(nnode == 8)					// for a hex, hard-code based on local point numbering
			{
				if(inode==0) { nbd[1] = true; nbd[3] = true; nbd[4] = true; }
				else if(inode==1) { nbd[0] = true; nbd[2] = true; nbd[5] = true; }
				else if(inode==2) { nbd[1] = nbd[3] = nbd[6] = true; }
				else if(inode==3) { nbd[2] = nbd[0] = nbd[7] = true; }
				else if(inode==4) { nbd[0] = nbd[5] = nbd[7] = true; }
				else if(inode==5) { nbd[1] = nbd[4] = nbd[6] = true; }
				else if(inode==6) { nbd[2] = nbd[5] = nbd[7] = true; }
				else if(inode==7) { nbd[3] = nbd[6] = nbd[4] = true; }
				else std::cout << "! UMesh3d: compute_topological(): Error in psup!" << std::endl;
			}

			for(int jnode = 0; jnode < nnode; jnode++)
			{
				int jpoin = inpoel(ielem, jnode);
				if(lpoin(jpoin) != ip && nbd[jnode] == true)
				{
					psup[ip].push_back(jpoin);
					lpoin(jpoin) = ip;
				}
			}
		}
		for(int i = 0; i < npoin; i++)
			lpoin(i,0) = -1;
	}
	//Points surrounding points is done.

	// 3. calculate number of edges using psup
	std::cout << "UMesh: compute_topological(): calculate number of edges using psup" << std::endl;
	nedge = 0; nbedge = 0;
	lpoin.zeros();

	for(int ipoin = 0; ipoin < npoin; ipoin++)
	{
		lpoin(ipoin) = 1;
		for(int jp = 0; jp < psup[ipoin].size(); jp++)
		{
			int jpoin = psup[ipoin].at(jp);

			if(lpoin(jpoin) != 1)
				nedge++;
		}
	}

	for(int i = 0; i < npoin; i++)
		lpoin(i) = 0;

	std::cout << "UMesh3d: compute_topological(): Number of edges = " << nedge << std::endl;

	elsed = new std::vector<int>[nedge];
	for(int i = 0; i < nedge; i++)
		elsed[i].reserve(8);
	edgepo.setup(nedge, nnoded);

	// 4. get edgepo
	std::cout << "UMesh3d: compute_topological(): Calculating edgepo" << std::endl;

	// first, boundary edges
	nbedge = 0;

	for(int ipoin = 0; ipoin < npoin; ipoin++)
	{
		lpoin(ipoin) = 1;
		for(int jp = 0; jp < psup[ipoin].size(); jp++)
		{
			int jpoin = psup[ipoin].at(jp);
			if(lpoin.get(jpoin) != 1 && flag_bpoin.get(ipoin)==1 && flag_bpoin.get(jpoin)==1)
			{
				edgepo(nbedge,0) = ipoin;
				edgepo(nbedge,1) = jpoin;
				nbedge++;
			}
		}
	}

	std::cout << "UMesh3d: compute_topological(): Number of boundary edges = " << nbedge << std::endl;
	for(int i = 0; i < npoin; i++)
		lpoin(i) = 0;

	//std::cout << "UMesh3d: compute_topological(): Calculating edgepo - interior" << std::endl;
	nedge = nbedge;
	for(int ipoin = 0; ipoin < npoin; ipoin++)
	{
		lpoin(ipoin) = 1;
		for(int jp = 0; jp < psup[ipoin].size(); jp++)
		{
			int jpoin = psup[ipoin].at(jp);
			if(lpoin(jpoin) != 1 && !(flag_bpoin.get(ipoin)==1 && flag_bpoin.get(jpoin)==1) )
			{
				edgepo(nedge,0) = ipoin;
				edgepo(nedge,1) = jpoin;
				nedge++;
			}
		}
	}

	// 5. Get elsed (elements surrounding each edge) using esup
	std::cout << "UMesh3d: compute_topological(): Calculating elsed" << std::endl;
	amat::Matrix<int> lelem(nelem,1);
	amc_int* ip = new amc_int[nnoded];

	for( ied = 0; ied < nedge; ied++)
	{
		for( i = 0; i < nnoded; i++)
			ip[i] = edgepo(ied,i);

		lelem.zeros();
		for( iel = esup_p(ip[0]); iel < esup_p(ip[0]+1); iel++)
		{
			lelem(esup(iel)) = 1;
		}

		for( jel = esup_p(ip[1]+1)-1; jel >= esup_p(ip[1]); jel--)
		{
			if(lelem(esup(jel)) == 1)
				elsed[ied].push_back(esup(jel));
		}
	}

	delete [] ip;

	// 6. Elements surrounding elements
	std::cout << "UMesh3d: compute_topological(): Elements surrounding elements...\n";

	esuel.setup(nelem, nfael);
	for(int ii = 0; ii < nelem; ii++)
		for(int jj = 0; jj < nfael; jj++)
			esuel(ii,jj) = -1;

	lpofa.setup(nfael, nnofa);	// lpofa(i,j) holds local node number of jth node of ith face (j in [0:nnofa], i in [0:nfael])

	if(nnode == 4)								// if tet
		for(int i = 0; i < nfael; i++)
		{
			for(int j = 0; j < nnofa; j++)
			{
				//lpofa(i,j) = perm(0,nnode-1,i,j);
				lpofa(i,j) = (i+j)%nnode;
			}
		}

	if(nnode == 8)								// if hex
	{
		lpofa(0,0) = 0; lpofa(0,1) = 3; lpofa(0,2) = 2; lpofa(0,3) = 1;
		lpofa(1,0) = 0; lpofa(1,1) = 1; lpofa(1,2) = 5; lpofa(1,3) = 4;
		lpofa(2,0) = 0; lpofa(2,1) = 4; lpofa(2,2) = 7; lpofa(2,3) = 3;
		lpofa(3,0) = 1; lpofa(3,1) = 2; lpofa(3,2) = 6; lpofa(3,3) = 5;
		lpofa(4,0) = 2; lpofa(4,1) = 3; lpofa(4,2) = 7; lpofa(4,3) = 6;
		lpofa(5,0) = 4; lpofa(5,1) = 5; lpofa(5,2) = 6; lpofa(5,3) = 7;
	}

	amat::Matrix<int> lhelp(nnofa,1);
	lhelp.zeros();
	lpoin.zeros();

	for(int ielem = 0; ielem < nelem; ielem++)
	{
		for(int ifael = 0; ifael < nfael; ifael++)
		{
			for(int i = 0; i < nnofa; i++)
			{
				lhelp(i) = inpoel(ielem, lpofa(ifael,i));		// lhelp stores global node nos. of current face of current element
				lpoin(lhelp(i)) = 1;
			}
			int ipoin = lhelp(0);								// ipoin is the global node number of first point of this face
			for(int istor = esup_p(ipoin); istor < esup_p(ipoin+1); istor++)
			{
				int jelem = esup(istor);
				if(jelem != ielem)
				{
					for(int jfael = 0; jfael < nfael; jfael++)
					{
						//Assume that no. of nodes in face ifael is same as that in face jfael
						int icoun = 0;
						for(int jnofa = 0; jnofa < nnofa; jnofa++)
						{
							int jpoin = inpoel(jelem, lpofa(jfael,jnofa));
							if(lpoin(jpoin)==1) icoun++;
						}
						if(icoun == nnofa)		// if all nnofa points of face ifael are found in face jfael, they both are the same face
						{
							esuel(ielem,ifael) = jelem;
							esuel(jelem,jfael) = ielem;
						}
					}
				}
			}
			for(int i = 0; i < nnofa; i++)
				lpoin(lhelp(i)) = 0;
		}
	}

	/** Computes, for each face, the elements on either side, the starting node and the ending node of the face. This is stored in intfac.
	The orientation of the face is such that the face points towards the element with larger index.
	NOTE: After the following portion, esuel holds (nelem + face no.) for each ghost cell, instead of -1 as before.*/

	//std::cout << "UMesh3d: compute_topological(): Computing intfac..." << std::endl;
	nbface = naface = 0;

	// first run: calculate nbface
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int je = esuel(ie,in);
			if(je == -1)
				nbface++;
		}
	}
	std::cout << "UMesh2d: compute_topological(): Number of boundary faces = " << nbface << std::endl;
	// calculate number of internal faces
	naface = nbface;
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int je = esuel(ie,in);
			if(je > ie && je < nelem)
				naface++;
		}
	}
	std::cout << "UMesh2d: compute_topological(): Number of all faces = " << naface << std::endl;

	//allocate intfac
	intfac.setup(naface,nnofa+2);

	//reset face totals
	nbface = naface = 0;

	//second run: populate intfac
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int* inp = new int[nnofa];
			for(int inofa = 0; inofa < nnofa; inofa++)
				inp[inofa] = lpofa.get(in,inofa);

			int je = esuel.get(ie,in);
			if(je == -1)
			{
				esuel(ie,in) = nelem+nbface;
				intfac(nbface,0) = ie;
				intfac(nbface,1) = nelem+nbface;
				for(int j = 0; j < nnofa; j++)
				{
					intfac(nbface,2+j) = inpoel(ie,inp[j]);
				}

				nbface++;
			}

			delete [] inp;
		}
	}

	naface = nbface;
	for(int ie = 0; ie < nelem; ie++)
	{
		for(int in = 0; in < nfael; in++)
		{
			int* inp = new int[nnofa];
			for(int inofa = 0; inofa < nnofa; inofa++)
				inp[inofa] = lpofa.get(in,inofa);

			int je = esuel(ie,in);
			if(je > ie && je < nelem)
			{
				intfac(naface,0) = ie;
				intfac(naface,1) = je;
				for(int j = 0; j < nnofa; j++)
					intfac(naface,2+j) = inpoel(ie,inp[j]);
				naface++;
			}

			delete [] inp;
		}
	}
}

// Computes topological properties of the boundary (surface) mesh
/*
 * - Boundary points (bpoints and bpointsinv)
 * - Boundary faces surrounding boudnary point (bfsubp and bfsubp_p)
 * - Boundary faces surrounding boundary face (bfsubf)
 * - Boundary faces adjoining boundary face (intbedge)
 */
void UMesh::compute_boundary_topological()
{
	// boundary data structures
	std::cout << "UMesh: compute_boundary_topological(): Computing bpoints, bpointsinv, bfsubp and bfsubf..." << std::endl;

	amc_int i, j, ipoin, jpoin;
	int ifa, iface, jface, istor, ied, jed, icoun, inode, jnode, inoded, jnoded, iedfa;

	nbpoin = 0;
	for(i = 0; i < npoin; i++)
	{
		if(flag_bpoin.get(i) == 1)
			nbpoin++;
	}
	bpoints.setup(nbpoin,1);
	bpointsinv.setup(npoin,1);

	nbpoin = 0;
	for(i = 0; i < npoin; i++)
	{
		bpointsinv(i) = -1;
		if(flag_bpoin.get(i) == 1)
		{
			bpoints(nbpoin) = i;
			bpointsinv(i) = nbpoin;
			nbpoin++;
		}
	}

	std::cout << "UMesh3d: compute_topological(): The mesh has " << nbpoin << " boundary points." << std::endl;

	amat::Matrix<int> lhelp(nnoded,1);
	amat::Matrix<int> lpoin(nbpoin,1);

	// bfaces surrounding bpoint

	bfsubp_p.setup(nbpoin+1,1);
	bfsubp_p.zeros();

	for(i = 0; i < nface; i++)
	{
		for(j = 0; j < nnofa; j++)
		{
			bfsubp_p(bpointsinv.get(bface(i,j))+1) += 1;		// bface(i,j) + 1 : the + 1 is there because the storage corresponding to the first node begins at 0, not at 1
		}
	}
	// Now make the members of bfsubp_p cumulative
	for(i = 1; i < nbpoin+1; i++)
		bfsubp_p(i) += bfsubp_p(i-1);
	// Now populate bfsubp
	bfsubp.setup(bfsubp_p(nbpoin),1);
	bfsubp.zeros();
	for(i = 0; i < nface; i++)
	{
		for(j = 0; j < nnofa; j++)
		{
			ipoin = bpointsinv.get(bface.get(i,j));
			bfsubp(bfsubp_p(ipoin)) = i;		// now put that element no. in the space pointed to by esup_p(ipoin)
			bfsubp_p(ipoin) += 1;				// an element corresponding to ipoin has been found - increment esup_p for that point
		}
	}
	//But now bfsubp_p holds increased values - each member increased by the number elements surrounding the corresponding point, so correct this.
	for(i = nbpoin; i >= 1; i--)
		bfsubp_p(i) = bfsubp_p(i-1);
	bfsubp_p(0) = 0;

	// bpoints surrounding bpoint
	
	std::cout << "UMesh2d: compute_boundary_topological(): Points surrounding points\n";
	bpsubp_p.setup(nbpoin+1,1,amat::ROWMAJOR);
	bpsubp_p.zeros();
	bpsubp_p(0) = 0;
	
	// lpoin: the ith member indicates the boundary point number of which the ith point is a surrounding point
	for(i = 0; i < nbpoin; i++) 
		lpoin(i) = -1;	// initialize this vector to -1
	istor = 0;
	
	std::vector<bool> nbd(nnofa);		// contains true if that local node number is connected to inode

	// first pass: calculate storage needed for psup
	for(int ip = 0; ip < nbpoin; ip++)
	{
		lpoin(ip) = ip;		// the point ip itself is not counted as a surrounding point of ip
		// Loop over elements surrounding this point
		for(ifa = bfsubp_p(ip); ifa < bfsubp_p(ip+1); ifa++)
		{
			iface = bfsubp(ifa,0);

			// find local node number of ip in iface
			for(jnode = 0; jnode < nnofa; jnode++)
				if(bface.get(iface,jnode) == ip) inode = jnode;

			for(j = 0; j < nnofa; j++)
				nbd[j] = false;

			if(nnofa == 3)
				for(i = 0; i < nbd.size(); i++)
					nbd[i] = true;
			else if(nnofa == 4)
				for(jnode = 0; jnode < nnofa; jnode++)
				{
					if(jnode == (inode+1)%nnofa /*perm(0,nnode-1,inode,1)*/ || jnode == perm(0,nnofa-1, inode, -1))
						nbd[jnode] = true;
				}

			//loop over nodes of the face
			for(inode = 0; inode < nnofa; inode++)
			{
				//Get global index of this node
				jpoin = bface.get(iface, inode);
				if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
				{
					istor++;
					lpoin(jpoin,0) = ip;		// set this point as a surrounding point of ip
				}
			}
		}
		bpsubp_p(ip+1) = istor;
	}

	bpsubp.setup(istor,1);

	//second pass: populate bpsubp
	istor = 0;
	for(i = 0; i < nbpoin; i++) 
		lpoin(i,0) = -1;	// initialize lpoin to -1
	
	for(int ip = 0; ip < nbpoin; ip++)
	{
		lpoin(ip,0) = ip;		// the point ip itself is not counted as a surrounding point of ip
		// Loop over elements surrounding this point
		for(ifa = bfsubp_p(ip); ifa < bfsubp_p(ip+1); ifa++)
		{
			iface = bfsubp(ifa);		// element number

			// find local node number of ip in ielem
			for(jnode = 0; jnode < nnofa; jnode++)
				if(bface.get(iface,jnode) == ip) 
					inode = jnode;

			for(j = 0; j < nnofa; j++)
				nbd[j] = false;

			if(nnofa == 3)
				for(i = 0; i < nbd.size(); i++)
					nbd[i] = true;
			else if(nnofa == 4)
				for(jnode = 0; jnode < nnofa; jnode++)
				{
					if(jnode == (inode+1)%nnofa || jnode == perm(0,nnode-1, inode, -1))
						nbd[jnode] = true;
				}

			//loop over nodes of the face
			for(inode = 0; inode < nnofa; inode++)
			{
				//Get global index of this node
				jpoin = bface.get(iface, inode);
				if(lpoin(jpoin,0) != ip && nbd[inode])		// test of this point as already been counted as a surrounding point of ip
				{
					bpsubp(istor,0) = jpoin;
					istor++;
					lpoin(jpoin) = ip;		// set this point as a surrounding point of ip
				}
			}
		}
	}

	istor = 0;

	// bfaces surrounding bface
	
	int ii, jj;
	bfsubf.setup(nface, nedfa);
	for(ii = 0; ii < nface; ii++)
		for(jj = 0; jj < nedfa; jj++)
			esuel(ii,jj) = -1;

	amat::Matrix<int> lpofab(nedfa, nnoded);			// lpofab(i,j) holds local node number of jth node of ith edge (j in [0:nnoded], i in [0:nedfa])
	for(int i = 0; i < nedfa; i++)
	{
		for(int j = 0; j < nnoded; j++)
		{
			lpofab(i,j) = (i+j) % nnofa;
		}
	}
	//lpofab.mprint();
	lhelp.zeros();
	lpoin.zeros();


	for(int iface = 0; iface < nface; iface++)
	{
		for(ied = 0; ied < nedfa; ied++)
		{
			for(i = 0; i < nnoded; i++)
			{
				lhelp(i) = bpointsinv.get(bface.get(iface, lpofab(ied,i)));	// lhelp stores global node nos. of current edge of current face
				lpoin(lhelp(i)) = 1;
			}
			ipoin = lhelp(0);
			for(istor = bfsubp_p(ipoin); istor < bfsubp_p(ipoin+1); istor++)
			{
				jface = bfsubp(istor);
				if(jface != iface)
				{
					for(jed = 0; jed < nedfa; jed++)
					{
						//Assume that no. of nodes in face ifael is same as that in face jfael
						icoun = 0;
						for(jnoded = 0; jnoded < nnoded; jnoded++)
						{
							jpoin = bpointsinv.get(bface.get(jface, lpofab(jed,jnoded)));
							if(lpoin(jpoin)==1) icoun++;
						}
						if(icoun == nnoded)		// nnoded is 2
						{
							bfsubf(iface,ied) = jface;
							bfsubf(jface,jed) = iface;
						}
					}
				}
			}
			for(i = 0; i < nnoded; i++)
				lpoin(lhelp(i)) = 0;
		}
	}


	// lbpoed(i,j) stores the local (bface) number of the jth node of the ith edge of a bface
	amat::Matrix<int> lbpoed(nedfa,nnoded);

	for(int i = 0; i < nedfa; i++)
		for(int j = 0; j < nnoded; j++)
			lbpoed(i,j) = (i+j)%nnofa;

	intbedge.setup(nbedge,nnoded+2);

	//reset b edge totals
	nbedge = 0;

	int* inp = new int[nnoded];
	nbedge = 0;
	for(iface = 0; iface < nface; iface++)
	{
		for( iedfa = 0; iedfa < nedfa; iedfa++)
		{
			for(inoded = 0; inoded < nnoded; inoded++)
				inp[inoded] = lbpoed.get(iedfa,inoded);

			jface = bfsubf(iface,iedfa);
			if(jface > iface && jface < nface)
			{
				intbedge(nbedge,0) = iface;
				intbedge(nbedge,1) = jface;
				for(j = 0; j < nnoded; j++)
					intbedge(nbedge,2+j) = bface(iface,inp[j]);
				nbedge++;
			}

		}
	}
	delete [] inp;
	std::cout << "UMesh3d: compute_boundary_topological(): Number of boundary edges = " << nbedge << std::endl;

	// finally, re-arrange intbedge so that its order matches that of edgepo and elsed
	// We're doing this the naive way.
	
	std::vector<int> temp(2+nnoded);
	for(ied = 0; ied < nbedge; ied++)
	{
		for(jed = ied; jed < nbedge; jed++)
		{
			// here we're assuming nnoded == 2, ie, linear mesh
			if(edgepo.get(ied,0)==intbedge.get(jed,2) && edgepo.get(ied,1)==intbedge.get(jed,3))
			{
				for(i = 0; i < 4; i++)
				{
					temp[i] = intbedge.get(ied,i);
					intbedge(ied,i) = intbedge.get(jed,i);
					intbedge(jed,i) = temp[i];
				}
			}
			else if(edgepo.get(ied,0)==intbedge.get(jed,3) && edgepo.get(ied,1)==intbedge.get(jed,2))
			{
				for(i = 0; i < 4; i++)
					temp[i] = intbedge.get(ied,i);

				// invert local point ordering and left-right face ordering
				intbedge(ied,0) = intbedge.get(jed,1);
				intbedge(ied,1) = intbedge.get(jed,0);
				intbedge(ied,2) = intbedge.get(jed,3);
				intbedge(ied,3) = intbedge.get(jed,2);

				for(i = 0; i < 4; i++)
					intbedge(jed,i) = temp[i];
			}
		}
	}

	// Now bedge. Note that this stores global point numbers, not boundary point numbers
	/*int inoded, in, iface, jface;
	amc_int nek = 0;
	bedge.setup(nbedge, 2 + nnoded);
	// first copy nodes
	for(amc_int ied = 0; ied < nbedge; ied++)
		for(inoded = 0; inoded < nnoded; inoded++)
			bedge(ied,inoded) = edgepo.get(ied,inoded);
	// now get bfaces using bfaces surrounding bpoint (bfsubp)
	int* inp = new int[nnoded];
	for(iface = 0; iface < nface; iface++)
	{
		for(in = 0; in < nedfa; in++)
		{
			for(inoded = 0; inoded < nnoded; inoded++)
				inp[inoded] = lpofab.get(in,inoded);

			jface = bfsubf(iface,in);
			if(jface > iface && jface < nface)
			{
				bedge(nek,nnoded+0) = iface;
				bedge(nek,nnoded+1) = jface;
				for(j = 0; j < nnoded; j++)
					bedge(nek,j) = bface(iface,inp[j]);
				nek++;
			}

		}
	}
	delete [] inp;**/

	std::cout << "UMesh3d: compute_boundary_topological(): Done." << std::endl;
}


/// Creates a UMesh object by adding one extra node at each edge centre, and also, in case of a hex mesh, one extra node at each face centre and cell centre.
/** \note Only works for tetrahedral and hexahedral elements.
 */
UMesh UMesh::convertLinearToQuadratic()
{
	int degree = 2;
	UMesh q;

	std::cout << "UMesh3d: convertLinearToQuadratic(): Producing a quadratic mesh from linear mesh..." << std::endl;
	if(nnofa != 3 && nnofa != 4) {
		std::cout << "! UMesh2d: convertLinearToQuadratic(): Mesh is not linear or is not supported!!" << std::endl;
		return q; }

	q.ndim = ndim;
	//q.npoin = npoin+nedge+naface+nelem;
	q.nelem = nelem;
	q.nface = nface; q.naface = naface; q.nbface = nbface;
	q.nedfa = nedfa;
	q.nedel = nedel;
	q.nfael = nfael;
	q.nnoded = nnoded + degree-1;
	q.nedge = nedge;
	q.nbedge = nbedge;

	if(nnode == 8)												// for hex
	{
		q.npoin = npoin + (degree-1)*nedge + (degree-1)*(degree-1)*naface + (degree-1)*(degree-1)*(degree-1)*nelem;
		q.nnode = nnode + nedel*(degree-1) + nfael*(degree-1)*(degree-1) + 1;
		q.nnofa = nnofa + nedfa*(degree-1) + (degree-1)*(degree-1);
	}
	else if(nnode == 4)											// for tet
	{
		q.npoin = npoin + (degree-1)*nedge + (degree-1)*(degree-2)/2*naface + 0*nelem;		// the nelem part is only true for 2nd and 3rd degree elements
		q.nnode = nnode + nedel*(degree-1);
		q.nnofa = nnofa + nedfa*(degree-1);
	}
	else  {
		std::cout << "! UMesh3d: convertLinearToQuadratic(): nnode is neither 4 nor 8 - mesh is not supported!" << std::endl;
		return q; }

	q.nbtag = nbtag; q.ndtag = ndtag;

	q.coords.setup(q.npoin, q.ndim);
	q.inpoel.setup(q.nelem, q.nnode);
	q.bface.setup(q.nface, q.nnofa+q.nbtag);
	q.vol_regions.setup(q.nelem, q.ndtag);
	//std::cout<< "Nodes per edge for P2 mesh = " << q.nnoded << ", number of edges " << q.nedge << std::endl;
	q.edgepo.setup(q.nedge,q.nnoded);
	q.flag_bpoin.setup(q.npoin,1);

	int ipoin, inode, idim, i, j, inofa, jnofa, ifnode, elem, ielem, iface, ibface;
	
	// copy nodes, elements and bfaces

	for(ipoin = 0; ipoin < npoin; ipoin++)
	{
		for(idim = 0; idim < ndim; idim++)
			q.coords(ipoin, idim) = coords(ipoin,idim);
	}

	for(int iel = 0; iel < nelem; iel++)
	{
		for(inode = 0; inode < nnode; inode++)
			q.inpoel(iel,inode) = inpoel(iel,inode);
		for(int itag = 0; itag < ndtag; itag++)
			q.vol_regions(iel,itag) = vol_regions(iel,itag);
	}

	for(iface = 0; iface < nface; iface++)
	{
		for(j = 0; j < nnofa; j++)
			q.bface(iface,j) = bface(iface,j);
		for(j = nnofa; j < nnofa+nbtag; j++)
			q.bface(iface,q.nnofa-nnofa+j) = bface(iface,j);
	}

	// copy edgepo for low-order nodes
	for(amc_int iedge = 0; iedge < nedge; iedge++)
		for(int inoded = 0; inoded < nnoded; inoded++)
			q.edgepo(iedge,inoded) = edgepo.get(iedge,inoded);

	double* centre = new double[ndim];

	if(nnode == 8)						// for hex mesh, we need face centres and cell centres
	{
		// get the body-centre nodes
		std::cout << "UMesh3d: convertLinearToQuadratic(): Adding nodes at cell-centres" << std::endl;

		for(int iel = 0; iel < nelem; iel++)
		{
			for(i = 0; i < ndim; i++)
				centre[i] = 0;

			for(inode = 0; inode < nnode; inode++)
			{
				for(idim = 0; idim < ndim; idim++)
					centre[idim] += coords(inpoel(iel,inode),idim);
			}

			for(idim = 0; idim < ndim; idim++)
				centre[idim] /= nnode;

			//add centre node to q.coords and update q.inpoel
			for(idim = 0; idim < ndim; idim++)
				q.coords(npoin+iel,idim) = centre[idim];

			q.inpoel(iel,q.nnode-1) = npoin+iel;
		}

		// face centre nodes
		std::cout << "UMesh3d: convertLinearToQuadratic(): Adding nodes at face centres." << std::endl;

		for(iface = 0; iface < nbface; iface++)
		{
			for(i = 0; i < ndim; i++)
				centre[i] = 0;

			for(ifnode = 2; ifnode < 2+nnofa; ifnode++)
				for(idim = 0; idim < ndim; idim++)
					centre[idim] += coords(intfac(iface,ifnode), idim);
			for(idim = 0; idim < ndim; idim++)
				centre[idim] /= nnofa;

			// add node to coords
			for(idim = 0; idim < ndim; idim++)
				q.coords(npoin+nelem+iface, idim) = centre[idim];

			// search for the bface corresponding to this face
			int rbface = -1; bool finmatch;
			bool* match = new bool[nnofa];

			for(int ibface = 0; ibface < nface; ibface++)
			{
				for(inofa = 0; inofa < nnofa; inofa++)
					match[inofa] = false;

				finmatch = true;

				for(inofa = 0; inofa < nnofa; inofa++)
				{
					for(jnofa = 0; jnofa < nnofa; jnofa++)
					{
						if(intfac(iface,inofa+2) == bface(ibface,jnofa)) {	// if ith node of iface matches any node of ibface, flag true
							match[inofa] = true;
							break; }
					}
				}

				/*for(int inofa = 0; inofa < nnofa; inofa++)
					std::cout << " " << match[inofa];
				std::cout << std::endl;*/

				for(inofa = 0; inofa < nnofa; inofa++)
				{
					if(match[inofa] == false) {					// if any node of iface did not match some node of ibface, there's no match
						finmatch = false;
						break;
					}
				}

				if(finmatch == true) rbface = ibface;
			}

			delete [] match;

			if(rbface == -1) std::cout << "! UMesh3d: convertLinearToQuadratic(): Bface corresponding to face " << iface << " not found!" << std::endl;

			// now add node to bface
			q.bface(rbface, q.nnofa-1) = npoin + nelem + iface;

			// --- add point to inpoel ---

			int nodenum = -1;
			elem = intfac(iface,0);

			bool ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,0)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 20;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// if found, this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,4)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 21;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,4)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 24;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,5)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 22;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,6)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 23;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,4) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,7)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 25;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}
		}

		// next, internal faces
		for(iface = nbface; iface < naface; iface++)
		{
			for(i = 0; i < ndim; i++)
				centre[i] = 0;

			for(ifnode = 2; ifnode < 2+nnofa; ifnode++)
				for(idim = 0; idim < ndim; idim++)
					centre[idim] += coords(intfac(iface,ifnode), idim);
			for(int idim = 0; idim < ndim; idim++)
				centre[idim] /= nnofa;

			// add node to coords
			for(idim = 0; idim < ndim; idim++)
				q.coords(npoin+nelem+iface, idim) = centre[idim];

			// add point to inpoel
			int nodenum = -1;
			elem = intfac(iface,0);

			bool ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,0)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 20;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				//continue;							// if found, this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,4)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 21;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				//continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,4)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 24;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				//continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,5)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 22;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				//continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,6)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 23;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				//continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,4) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,7)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 25;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				//continue;							// this face is done
			}

			// other element
			nodenum = -1;
			elem = intfac(iface,1);

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,0)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 20;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// if found, this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,4)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 21;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,0) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,4)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 24;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,1) || intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,5)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 22;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;							// this face is done
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,2) || intfac(iface,2+inode)==inpoel(elem,3) || intfac(iface,2+inode)==inpoel(elem,7) || intfac(iface,2+inode)==inpoel(elem,6)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 23;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;
			}

			ematch = true;
			for(inode = 0; inode < nnofa; inode++)
			{
				if(!(intfac(iface,2+inode)==inpoel(elem,4) || intfac(iface,2+inode)==inpoel(elem,5) || intfac(iface,2+inode)==inpoel(elem,6) || intfac(iface,2+inode)==inpoel(elem,7)))
					{ ematch = false; break; }
			}
			if(ematch==true)
			{
				nodenum = 25;
				q.inpoel(elem,nodenum) = npoin + nelem + iface;
				continue;
			}
		}

		// --- next, add points at edge centres ---
		std::cout << "UMesh3d: convertLinearToQuadratic(): Adding points at edge centres for hexes" << std::endl;

		// first, boundary edges
		for(int ied = 0; ied < nbedge; ied++)
		{
			for(i = 0; i < ndim; i++)
				centre[i] = 0;

			for(ifnode = 0; ifnode < nnoded; ifnode++)
				for(idim = 0; idim < ndim; idim++)
					centre[idim] += coords(edgepo(ied,ifnode), idim);
			for(idim = 0; idim < ndim; idim++)
				centre[idim] /= nnoded;

			// add node to coords
			int cono = npoin+nelem+naface+ied;
			for(idim = 0; idim < ndim; idim++)
				q.coords(cono, idim) = centre[idim];

			// add to elements surrounding edge
			//std::cout << "add to elements surr edge" << std::endl;
			for(ielem = 0; ielem < elsed[ied].size(); ielem++)
			{
				elem = elsed[ied].at(ielem);

				if((edgepo.get(ied,0)==inpoel(elem,0)&&edgepo.get(ied,1)==inpoel(elem,1)) || (edgepo.get(ied,1)==inpoel(elem,0)&&edgepo.get(ied,0)==inpoel(elem,1)))
					q.inpoel(elem,8) = cono;
				else if((edgepo.get(ied,0)==inpoel(elem,1)&&edgepo.get(ied,1)==inpoel(elem,2)) || (edgepo.get(ied,1)==inpoel(elem,1)&&edgepo.get(ied,0)==inpoel(elem,2)))
					q.inpoel(elem,9) = cono;
				else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,3)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,3)))
					q.inpoel(elem,10) = cono;
				else if((edgepo(ied,0)==inpoel(elem,3)&&edgepo(ied,1)==inpoel(elem,0)) || (edgepo(ied,1)==inpoel(elem,3)&&edgepo(ied,0)==inpoel(elem,0)))
					q.inpoel(elem,11) = cono;
				else if((edgepo(ied,0)==inpoel(elem,0)&&edgepo(ied,1)==inpoel(elem,4)) || (edgepo(ied,1)==inpoel(elem,0)&&edgepo(ied,0)==inpoel(elem,4)))
					q.inpoel(elem,12) = cono;
				else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,5)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,5)))
					q.inpoel(elem,13) = cono;
				else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,6)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,6)))
					q.inpoel(elem,14) = cono;
				else if((edgepo(ied,0)==inpoel(elem,3)&&edgepo(ied,1)==inpoel(elem,7)) || (edgepo(ied,1)==inpoel(elem,3)&&edgepo(ied,0)==inpoel(elem,7)))
					q.inpoel(elem,15) = cono;
				else if((edgepo(ied,0)==inpoel(elem,4)&&edgepo(ied,1)==inpoel(elem,5)) || (edgepo(ied,1)==inpoel(elem,4)&&edgepo(ied,0)==inpoel(elem,5)))
					q.inpoel(elem,16) = cono;
				else if((edgepo(ied,0)==inpoel(elem,5)&&edgepo(ied,1)==inpoel(elem,6)) || (edgepo(ied,1)==inpoel(elem,5)&&edgepo(ied,0)==inpoel(elem,6)))
					q.inpoel(elem,17) = cono;
				else if((edgepo(ied,0)==inpoel(elem,6)&&edgepo(ied,1)==inpoel(elem,7)) || (edgepo(ied,1)==inpoel(elem,6)&&edgepo(ied,0)==inpoel(elem,7)))
					q.inpoel(elem,18) = cono;
				else if((edgepo(ied,0)==inpoel(elem,7)&&edgepo(ied,1)==inpoel(elem,4)) || (edgepo(ied,1)==inpoel(elem,7)&&edgepo(ied,0)==inpoel(elem,4)))
					q.inpoel(elem,19) = cono;
			}

			// add to edgepo
			q.edgepo(ied,nnoded) = cono;

			//std::cout << "find bface" << std::endl;
			// find bfaces that this edge belongs to
			std::vector<int> edfa;
			bool bmatch1, bmatch2;
			for(int ibface = 0; ibface < nface; ibface++)
			{
				bmatch1 = bmatch2 = false;
				for(inode = 0; inode < nnofa; inode++)
					if(edgepo(ied,0)==bface(ibface,inode))
						bmatch1 = true;
				for(inode = 0; inode < nnofa; inode++)
					if(edgepo(ied,1)==bface(ibface,inode))
						bmatch2 = true;

				if(bmatch1 && bmatch2) edfa.push_back(ibface);
			}

			//std::cout << "add new point to bfaces" << std::endl;
			// add new point to the bfaces that were found
			for(int ibf = 0; ibf < edfa.size(); ibf++)
			{
				ibface = edfa.at(ibf);

				if((edgepo.get(ied,0)==bface.get(ibface,0) && edgepo.get(ied,1)==bface.get(ibface,1)) || (edgepo.get(ied,1)==bface.get(ibface,0) && edgepo.get(ied,0)==bface.get(ibface,1)))
					q.bface(ibface,4) = cono;
				else if((edgepo(ied,0)==bface(ibface,1) && edgepo(ied,1)==bface(ibface,2)) || (edgepo(ied,1)==bface(ibface,1) && edgepo(ied,0)==bface(ibface,2)))
					q.bface(ibface,5) = cono;
				else if((edgepo(ied,0)==bface(ibface,2) && edgepo(ied,1)==bface(ibface,3)) || (edgepo(ied,1)==bface(ibface,2) && edgepo(ied,0)==bface(ibface,3)))
					q.bface(ibface,6) = cono;
				else if((edgepo(ied,0)==bface(ibface,3) && edgepo(ied,1)==bface(ibface,0)) || (edgepo(ied,1)==bface(ibface,3) && edgepo(ied,0)==bface(ibface,0)))
					q.bface(ibface,7) = cono;
			}
		}

		// internal edges
		//std::cout << "internal edges" << std::endl;
		for(int ied = nbedge; ied < nedge; ied++)
		{
			for(i = 0; i < ndim; i++)
				centre[i] = 0;

			for(ifnode = 0; ifnode < nnoded; ifnode++)
				for(idim = 0; idim < ndim; idim++)
					centre[idim] += coords(edgepo(ied,ifnode), idim);
			for(idim = 0; idim < ndim; idim++)
				centre[idim] /= nnoded;

			// add node to coords
			int cono = npoin+nelem+naface+ied;
			for(idim = 0; idim < ndim; idim++)
				q.coords(cono, idim) = centre[idim];

			// add to elements surrounding edge
			for(ielem = 0; ielem < elsed[ied].size(); ielem++)
			{
				elem = elsed[ied].at(ielem);

				if((edgepo(ied,0)==inpoel(elem,0)&&edgepo(ied,1)==inpoel(elem,1)) || (edgepo(ied,1)==inpoel(elem,0)&&edgepo(ied,0)==inpoel(elem,1)))
					q.inpoel(elem,8) = cono;
				else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,2)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,2)))
					q.inpoel(elem,9) = cono;
				else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,3)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,3)))
					q.inpoel(elem,10) = cono;
				else if((edgepo(ied,0)==inpoel(elem,3)&&edgepo(ied,1)==inpoel(elem,0)) || (edgepo(ied,1)==inpoel(elem,3)&&edgepo(ied,0)==inpoel(elem,0)))
					q.inpoel(elem,11) = cono;
				else if((edgepo(ied,0)==inpoel(elem,0)&&edgepo(ied,1)==inpoel(elem,4)) || (edgepo(ied,1)==inpoel(elem,0)&&edgepo(ied,0)==inpoel(elem,4)))
					q.inpoel(elem,12) = cono;
				else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,5)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,5)))
					q.inpoel(elem,13) = cono;
				else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,6)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,6)))
					q.inpoel(elem,14) = cono;
				else if((edgepo(ied,0)==inpoel(elem,3)&&edgepo(ied,1)==inpoel(elem,7)) || (edgepo(ied,1)==inpoel(elem,3)&&edgepo(ied,0)==inpoel(elem,7)))
					q.inpoel(elem,15) = cono;
				else if((edgepo(ied,0)==inpoel(elem,4)&&edgepo(ied,1)==inpoel(elem,5)) || (edgepo(ied,1)==inpoel(elem,4)&&edgepo(ied,0)==inpoel(elem,5)))
					q.inpoel(elem,16) = cono;
				else if((edgepo(ied,0)==inpoel(elem,5)&&edgepo(ied,1)==inpoel(elem,6)) || (edgepo(ied,1)==inpoel(elem,5)&&edgepo(ied,0)==inpoel(elem,6)))
					q.inpoel(elem,17) = cono;
				else if((edgepo(ied,0)==inpoel(elem,6)&&edgepo(ied,1)==inpoel(elem,7)) || (edgepo(ied,1)==inpoel(elem,6)&&edgepo(ied,0)==inpoel(elem,7)))
					q.inpoel(elem,18) = cono;
				else if((edgepo(ied,0)==inpoel(elem,7)&&edgepo(ied,1)==inpoel(elem,4)) || (edgepo(ied,1)==inpoel(elem,7)&&edgepo(ied,0)==inpoel(elem,4)))
					q.inpoel(elem,19) = cono;
			}
			
			// add to edgepo
			q.edgepo(ied,nnoded) = cono;
		}
		delete [] centre;
		return q;
	}		// end if for hex

	// Tets

	// first, boundary edges
	for(int ied = 0; ied < nbedge; ied++)
	{
		for(int i = 0; i < ndim; i++)
			centre[i] = 0;

		for(int ifnode = 0; ifnode < nnoded; ifnode++)
			for(int idim = 0; idim < ndim; idim++)
				centre[idim] += coords(edgepo(ied,ifnode), idim);
		for(int idim = 0; idim < ndim; idim++)
			centre[idim] /= nnoded;

		// add node to coords
		int cono = npoin+ied;
		for(int idim = 0; idim < ndim; idim++)
			q.coords(cono, idim) = centre[idim];

		// add to elements surrounding edge NOTE: ordering of nodes is taken from Gmsh docs
		for(int ielem = 0; ielem < elsed[ied].size(); ielem++)
		{
			int elem = elsed[ied][ielem];

			if((edgepo(ied,0)==inpoel(elem,0)&&edgepo(ied,1)==inpoel(elem,1)) || (edgepo(ied,1)==inpoel(elem,0)&&edgepo(ied,0)==inpoel(elem,1)))
				q.inpoel(elem,4) = cono;
			else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,2)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,2)))
				q.inpoel(elem,5) = cono;
			else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,3)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,3)))
				q.inpoel(elem,8) = cono;
			else if((edgepo(ied,0)==inpoel(elem,3)&&edgepo(ied,1)==inpoel(elem,0)) || (edgepo(ied,1)==inpoel(elem,3)&&edgepo(ied,0)==inpoel(elem,0)))
				q.inpoel(elem,7) = cono;
			else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,0)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,0)))
				q.inpoel(elem,6) = cono;
			else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,3)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,3)))
				q.inpoel(elem,9) = cono;
		}

		// add to edgepo
		q.edgepo(ied,nnoded) = cono;

		// find bfaces that this edge belongs to
		std::vector<int> edfa;
		bool bmatch1, bmatch2;
		for(int ibface = 0; ibface < nface; ibface++)
		{
			bmatch1 = bmatch2 = false;

			for(int inode = 0; inode < nnofa; inode++)
				if(edgepo(ied,0)==bface(ibface,inode))
					bmatch1 = true;
			for(int inode = 0; inode < nnofa; inode++)
				if(edgepo(ied,1)==bface(ibface,inode))
					bmatch2 = true;

			if(bmatch1 && bmatch2) edfa.push_back(ibface);
		}

		// add new point to found bfaces
		for(int ibf = 0; ibf < edfa.size(); ibf++)
		{
			int ibface = edfa.at(ibf);

			if((edgepo(ied,0)==bface(ibface,0) && edgepo(ied,1)==bface(ibface,1)) || (edgepo(ied,1)==bface(ibface,0) && edgepo(ied,0)==bface(ibface,1)))
				q.bface(ibface,3) = cono;
			else if((edgepo(ied,0)==bface(ibface,1) && edgepo(ied,1)==bface(ibface,2)) || (edgepo(ied,1)==bface(ibface,1) && edgepo(ied,0)==bface(ibface,2)))
				q.bface(ibface,4) = cono;
			else if((edgepo(ied,0)==bface(ibface,2) && edgepo(ied,1)==bface(ibface,0)) || (edgepo(ied,1)==bface(ibface,2) && edgepo(ied,0)==bface(ibface,0)))
				q.bface(ibface,5) = cono;
		}
	}

	// internal edges
	for(int ied = nbedge; ied < nedge; ied++)
	{
		for(int i = 0; i < ndim; i++)
			centre[i] = 0;

		for(int ifnode = 0; ifnode < nnoded; ifnode++)
			for(int idim = 0; idim < ndim; idim++)
				centre[idim] += coords(edgepo(ied,ifnode), idim);
		for(int idim = 0; idim < ndim; idim++)
			centre[idim] /= nnoded;

		// add node to coords
		int cono = npoin+ied;
		for(int idim = 0; idim < ndim; idim++)
			q.coords(cono, idim) = centre[idim];

		// add to elements surrounding edge NOTE: ordering of nodes is taken from Gmsh docs
		for(int ielem = 0; ielem < elsed[ied].size(); ielem++)
		{
			int elem = elsed[ied].at(ielem);

			if((edgepo(ied,0)==inpoel(elem,0)&&edgepo(ied,1)==inpoel(elem,1)) || (edgepo(ied,1)==inpoel(elem,0)&&edgepo(ied,0)==inpoel(elem,1)))
				q.inpoel(elem,4) = cono;
			else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,2)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,2)))
				q.inpoel(elem,5) = cono;
			else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,3)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,3)))
				q.inpoel(elem,8) = cono;
			else if((edgepo(ied,0)==inpoel(elem,3)&&edgepo(ied,1)==inpoel(elem,0)) || (edgepo(ied,1)==inpoel(elem,3)&&edgepo(ied,0)==inpoel(elem,0)))
				q.inpoel(elem,7) = cono;
			else if((edgepo(ied,0)==inpoel(elem,2)&&edgepo(ied,1)==inpoel(elem,0)) || (edgepo(ied,1)==inpoel(elem,2)&&edgepo(ied,0)==inpoel(elem,0)))
				q.inpoel(elem,6) = cono;
			else if((edgepo(ied,0)==inpoel(elem,1)&&edgepo(ied,1)==inpoel(elem,3)) || (edgepo(ied,1)==inpoel(elem,1)&&edgepo(ied,0)==inpoel(elem,3)))
				q.inpoel(elem,9) = cono;
		}
			
		// add to edgepo
		q.edgepo(ied,nnoded) = cono;
	}
	
	// set flag_bpoin
	q.flag_bpoin.zeros();
	for(int i = 0; i < q.nface; i++)
		for(int j = 0; j < q.nnofa; j++)
			q.flag_bpoin(q.bface(i,j)) = 1;

	// set nbpoin
	q.nbpoin = 0;
	for(int i = 0; i < q.npoin; i++)
		q.nbpoin += q.flag_bpoin.get(i);

	delete [] centre;
	std::cout << "UMesh3d: convertLinearToQuadratic(): Quadratic mesh produced." << std::endl;
	return q;
}

}	// end namespace acfd
