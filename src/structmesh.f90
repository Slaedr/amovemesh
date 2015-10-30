!* This program generates a mixed grid in a unit cubic domain, and generates
!* grid and boundary condition files in an unstructured format (.ugrid).
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* The domain is divided into two parts, upper and lower.
!* A tetrahedral grid is generated in the upper part and
!* a prismatic grid is generated in the lower part.
!* It resembles a typical high-Reynolds-number viscous grid.
!*
!* This is Version 3 (2011).
!*
!* This F90 program is written and made available for an educational purpose.
!*
!* This file may be updated in future.
!*
!* Katate Masatsuka, March 2011. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* Input : nx, ny, nz ! Number of nodes in each coordinate direction (nz>>nx=ny)
!*                 zm ! Z-coordinate separating tetra and prismatic layers.
!*
!* Output: mixgrid.ugrid     ! Grid file (in .ugrid format)
!*         mixgrid.mapbc     ! Boundary condition file
!*         mixgrid.lines_fmt ! Line information in the prismatic layer
!*         mixgrid_tecplot.dat          !  Volume grid file for viewing by Tecplot
!*         mixgrid_boundary_tecplot.dat !Boundary grid file for viewing by Tecplot
!*
!* Note: Information on UGRID format can be found at
!*       http://www.simcenter.msstate.edu/docs/solidmesh/ugridformat.html
!*
!* Note: mixgrid.ugrid and mixgrid.mapbc are the files to be read by a flow solver
!*       for running a simulation.
!*
!* Note: BC numbers in mixgrid.mapbc must be chosen from available numbers in the
!*       flow solver.
!*
!********************************************************************************
program mixgrid_cube

  implicit none
  integer , parameter :: dp = selected_real_kind(15) !Double precision

! Structured data
  real(dp), allocatable, dimension(:,:,:) :: x, y, z

! Mapping from (i,j,k) to the node number
  integer , allocatable, dimension(:,:,:) :: ijk2n

! Unstructured data: these data will be used to write output files.
  real(dp), allocatable, dimension(:,:) :: xyz       ! Coordinates
  integer , allocatable, dimension(:,:) :: prz, tet, tria, quad,bfact,hex ! Element connectivities

! Number of nodes in each coordinate direction
  integer :: nx, ny, nz, nx1,ny1,nz1,nbin
  integer :: i0,j0,k0,i1,j1,k1,i2,j2,k2
  integer :: iwa0,iwa1
! Number of nodes in z-direction: lower and upper.
  integer :: nm, nz_lower, nz_upper
! Number of prisms, tetrahedra, triangles, quadrilaterals
  integer :: nprz, ntet, ntria, nquad,nhexa
  integer :: ntetr
  integer :: ntet1,ntet2,ntet3,ntet4,ntet5,ntet6,nboun,ncount
! Number of nodes
  integer :: nnodes,nnode
! Local numbering of a hexahedron
  integer :: hex_loc(27)
! Grid spacing in each coordinate direction
  real(dp) :: dx, dy,dz
! Z-coordinate of the intersection between tetrahedra and prismatic parts.
  real(dp) :: zm
! Spacing in z for the upper and lower regions (uniform in each region).
  real(dp) :: dz_lower, dz_upper,pi
! Local variables
  integer  :: i, j, k, os,ie
  integer  :: inn(27)
  integer,allocatable, dimension(:) :: inpt
  real*8   :: x1(3,9),x2(3,27)
  real(dp) :: zs
! Line info
  integer  :: n_lines, n_total_points, max_points, min_points, i_count, inode

!*******************************************************************************
! Unit cube to be gridded.
!        ______________ 
!       /             /|      Boundary group numbers are assigned as below:
!      /             / |       1:  Boundary face on x = 0.0
! z=1 /____________ /  |       2:  Boundary face on x = 1.0
!     |  |         |   |       3:  Boundary face on y = 0.0
!     |  |         |   |       4:  Boundary face on y = 1.0
!     |  |         |   |       5:  Boundary face on z = 0.0
!     |  |y=1______|___|       6:  Boundary face on z = 1.0
!     |  /         |  /
!     | /          | /
!     |/___________|/
! x=y=z=0        x=1
!
!
!  Here is the side-view of the grid.
!   Isotropic tetrahedra grid in the region: z > zm (node index k =  1, nm)
!              Prismatic grid in the region: z < zm (node index k = nm, nz)
!
!                                                        Index in z:
!                    z ^
!                      |  1 ____________________________ k=nz=(nz_upper-1)+nz_lower
!                      |   !                            !
!                      |   !                            !
!                      |   !                            !
!                      |   |                            |
!                      |   | Isotropic Tetrahedral grid |
!                      |   |                            |
!                      |   |                            |
!                      |   |                            |
!                      |   |                            |
!                      | zm|____________________________|k = nm = nz_lower
!                      |   |                            |
!                      |   |                            |
!                      |   |          Prismatic         |
!                      |   |                            |
!                      |   |____________________________|
!                      |  0                              1
!                      -------------------------------------> x
!
!
!*******************************************************************************
! Number of points in each direction.
! Set these values such that nz >> nx = ny.
! Also specify the vertical limit of the prismatic region by zm.

  open(2,file='bump3dgrid.p3dfmt')
  read(2,*) nbin
  read(2,*) nx,ny,nz
  allocate( x(nx,ny,nz), y(nx,ny,nz), z(nx,ny,nz) )
  read(2,*) (((x(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
             (((y(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
             (((z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  close(2)

!...For 3d bump
  iwa0 = 13
  iwa1 = 33

  iwa0 = 25
  iwa1 = 65

  iwa0 = 49
  iwa1 = 129

 ! iwa0 = 97
 ! iwa1 = 257

 ! iwa0 = 385
 ! iwa1 = 1025

  !iwa0 = 81
  !iwa1 = 129

!
!...get the mid point
!

! do i = 1,nx
!  do j = 1,ny
!  do k = 1,(nz-1)/2
  
!     x(i,j,2*k) = 0.5d0*(x(i,j,2*k-1)+x(i,j,2*k+1))
!     y(i,j,2*k) = 0.5d0*(y(i,j,2*k-1)+y(i,j,2*k+1))
!     z(i,j,2*k) = 0.5d0*(z(i,j,2*k-1)+z(i,j,2*k+1))
   
! enddo
! enddo
! enddo
! Allocate structured arrays

!
!(nx-1),(ny-1),(nz-1)should be able to divided by 2.
!

!The curved cell nodes number
  nx1 = (nx-1)/2 + 1
  ny1 = (ny-1)/2 + 1
  nz1 = (nz-1)/2 + 1

print*,nx,ny,nz
!
! Read from the existing code
!

! For smooth transition from prismatic layer to the tetrahedral region,
! you may want to apply some stretching in z-direction in the prismatic
! layer. I'll leave it to you.

! Construct a node array

! (1)Construct a map first.

    allocate(ijk2n(nx,ny,nz)) ! (i,j,k) to node map

     nnodes = 0

    do i = 1, nx
     do j = 1, ny
      do k = 1, nz
        nnodes = nnodes + 1
        ijk2n(i,j,k) = nnodes
      end do
     end do
    end do

! (2)Construct a node array using ijk2n(:,:,:)
!     xyz(i,1) = x at node i.
!     xyz(i,2) = y at node i.
!     xyz(i,3) = z at node i.

    allocate(xyz(nnodes,3))

    allocate(inpt(nnodes))

    do i = 1, nx
     do j = 1, ny
      do k = 1, nz
       xyz(ijk2n(i,j,k), 1) = x(i,j,k)
       xyz(ijk2n(i,j,k), 2) = y(i,j,k)
       xyz(ijk2n(i,j,k), 3) = z(i,j,k)
      end do
     end do
    end do

! Nodes are now ordered from 1 to nnodes(=total # of nodes).
! The nodal coordinates are stored in a one-dimensional fashion:
! x_i = xyz(i,1), y_i = xyz(i,2), z_i = xyz(i,3), i=1,2,3,...,nnodes.

! So, we don't need x, y, z, in the (i,j,k) data format.
! Goodbye, x, y, z!

   deallocate(x,y,z)

! Now, construct tet element and tria boundary connectivity data:
!  Prz: 6 nodes that define the prism.
! Quad: 4 nodes that define the quad face.
!  Tet: 4 nodes that define the tetrahedron.
! Tria: 3 nodes that define the tria face.

! Allocate the arrays with known dimensions:
!  - Hexahedron is divided into 6 tetrahedra
!  - Quad faces are divided into 2 triangles.

  nhexa  =(nx1-1)*(ny1-1)*(nz1-1)
  nquad = 2*(nx1-1)*(ny1-1) + 2*(ny1-1)*(nz1-1) + 2*(nz1-1)*(nx1-1)

  allocate(  hex(nhexa,27), quad(nquad,12) , bfact(nquad,11))

! Reset ntet and ntria, nhex and nquad.

  nquad = 0
  nhexa = 0

! Loop over (i,j,k) and construct connectivity data.

   i_loop : do i = 1, nx1-1
    j_loop : do j = 1, ny1-1
     k_loop : do k = 1, nz1-1

! At each step, we look at a hexahedron as shown below.
! We work with local numbering because it is convenient.
!
!  Local node numbering, {1,2,3,4,5,6,7,8}, is defined as follows.
!
!                         (i,j+1,k+1) 6----------------------7 (i+1,j+1,k+1)
!                                    /.                     /|
!                                   / .                    / |
!                                  /  .                   /  |
!                                 /   .                  /   |
!                                /    .                 /    |
!                               /     .                /     |
!                              /      .               /      |
!                             /       .              /       |
!                            /        .             /        |
!                           /         .            /         |
!                          /          .           /          |
!               (i,j,k+1) 2----------------------3 (i+1,j,k+1)
!                         !           .          !           !
!                         !           .          !           !
!                         !           .          !           !
!                         |           5..........|...........8
!                         | (i,j+1,k).           |          /(i+1,j+1,k)
!                         |         .            |         /
!                         |        .             |        /
!                         |       .              |       /
!                         |      .               |      /
!                         |     .                |     /
!                         |    .                 |    /
!                         |   .                  |   /
!                         |  .                   |  /
!                         | .                    | /
!                         |.                     |/
!                         1----------------------4
!                   (i,j,k)                     (i+1,j,k)
!

! Store the node numbers into hex_loc(1:8), so that
! we can look at the hex as defined by local numbers {1,2,3,4,5,6,7,8}.
!
     i0 = (i-1)*2 + 1
     j0 = (j-1)*2 + 1
     k0 = (k-1)*2 + 1
 
     i1 = i0 + 1
     j1 = j0 + 1
     k1 = k0 + 1    
  
     i2 = (i)*2 + 1
     j2 = (j)*2 + 1
     k2 = (k)*2 + 1

!    Lower 4 vertices

     hex_loc(1) = ijk2n(i0  , j0  , k0  )
     hex_loc(2) = ijk2n(i2  , j0  , k0  )
     hex_loc(3) = ijk2n(i2  , j2  , k0  )
     hex_loc(4) = ijk2n(i0  , j2  , k0  )

!    Upper 4 vertices
     hex_loc(5) = ijk2n(i0  , j0  , k2)
     hex_loc(6) = ijk2n(i2  , j0  , k2)
     hex_loc(7) = ijk2n(i2, j2, k2)
     hex_loc(8) = ijk2n(i0  , j2, k2)
!
!
!
!    Lower 4 vertices

     hex_loc(9) = ijk2n(i1  , j0  , k0  )
     hex_loc(10) = ijk2n(i2  , j1  , k0  )
     hex_loc(11) = ijk2n(i1  , j2  , k0  )
     hex_loc(12) = ijk2n(i0  , j1  , k0  )

     hex_loc(13) = ijk2n(i0  , j0  , k1)
     hex_loc(14) = ijk2n(i2  , j0  , k1)
     hex_loc(15) = ijk2n(i2  , j2  , k1)
     hex_loc(16) = ijk2n(i0  , j2  , k1)

     hex_loc(17) = ijk2n(i1  , j0  , k2  )
     hex_loc(18) = ijk2n(i2  , j1  , k2  )
     hex_loc(19) = ijk2n(i1  , j2  , k2  )
     hex_loc(20) = ijk2n(i0  , j1  , k2  )

     hex_loc(21) = ijk2n(i1  , j1  , k0)
     hex_loc(22) = ijk2n(i1  , j0  , k1)
     hex_loc(23) = ijk2n(i2  , j1  , k1)
     hex_loc(24) = ijk2n(i1  , j2, k1)

     hex_loc(25) = ijk2n(i0  , j1  , k1  )
     hex_loc(26) = ijk2n(i1  , j1  , k2  )
     hex_loc(27) = ijk2n(i1  , j1  , k1  )
! From here, we can think of the hex defined by the
! numbers {1,2,3,4,5,6,7,8}.

!--- Hexahedron defined by the nodes 1265-4378

!     Hexa #1
      nhexa = nhexa + 1
      hex(nhexa,1:27) = hex_loc(1:27)

!--- Boundary faces if the hexa is adjacent to the boundary:
!    Note: Nodes are ordered so that the boundary face points INWARD.
!
!    1. Left boundary face (x=0)
      if     (i == 1) then
        
      nquad = nquad + 1
      quad(nquad,1) = hex_loc(1)
      quad(nquad,2) = hex_loc(4)
      quad(nquad,3) = hex_loc(8)
      quad(nquad,4) = hex_loc(5)
      quad(nquad,5) = hex_loc(12)
      quad(nquad,6) = hex_loc(16)
      quad(nquad,7) = hex_loc(20)
      quad(nquad,8) = hex_loc(13)
      quad(nquad,9) = hex_loc(25)
      quad(nquad,10) = 4 !<------ Boundary group number
      quad(nquad,11) = nhexa
      quad(nquad,12) = 1
!    2. Right boundary face (x=1)
      elseif (i == nx1-1) then
      nquad = nquad + 1
      quad(nquad,1) = hex_loc(3)
      quad(nquad,2) = hex_loc(2)
      quad(nquad,3) = hex_loc(6)
      quad(nquad,4) = hex_loc(7)
      quad(nquad,5) = hex_loc(10)
      quad(nquad,6) = hex_loc(14)
      quad(nquad,7) = hex_loc(18)
      quad(nquad,8) = hex_loc(15)
      quad(nquad,9) = hex_loc(23)
      quad(nquad,10) = 4  !<------ Boundary group number
      quad(nquad,11) = nhexa
       quad(nquad,12) = 1
      endif

!    3. Front boundary face (y=0)
      if     (j == 1) then
      nquad = nquad + 1
      quad(nquad,1) = hex_loc(2)
      quad(nquad,2) = hex_loc(1)
      quad(nquad,3) = hex_loc(5)
      quad(nquad,4) = hex_loc(6)
      quad(nquad,5) = hex_loc(9)
      quad(nquad,6) = hex_loc(13)
      quad(nquad,7) = hex_loc(17)
      quad(nquad,8) = hex_loc(14)
      quad(nquad,9) = hex_loc(22)
      quad(nquad,10) = 3  !<------ Boundary group number
      quad(nquad,11) = nhexa
      quad(nquad,12) = 0
!    4. Back boundary face (y=1)
      nquad = nquad + 1
      quad(nquad,1) = hex_loc(4)
      quad(nquad,2) = hex_loc(3)
      quad(nquad,3) = hex_loc(7)
      quad(nquad,4) = hex_loc(8)
      quad(nquad,5) = hex_loc(11)
      quad(nquad,6) = hex_loc(15)
      quad(nquad,7) = hex_loc(19)
      quad(nquad,8) = hex_loc(16)
      quad(nquad,9) = hex_loc(24)
      quad(nquad,10) = 3  !<------ Boundary group number
      quad(nquad,11) = nhexa
      quad(nquad,12) = 0
      endif

!    5. Bottom boundary face (z=0)
      if     (k == 1) then
      nquad = nquad + 1
      quad(nquad,1) = hex_loc(1)
      quad(nquad,2) = hex_loc(2)
      quad(nquad,3) = hex_loc(3)
      quad(nquad,4) = hex_loc(4)
      quad(nquad,5) = hex_loc(9)
      quad(nquad,6) = hex_loc(10)
      quad(nquad,7) = hex_loc(11)
      quad(nquad,8) = hex_loc(12)
      quad(nquad,9) = hex_loc(21)
      quad(nquad,10) = 3
      quad(nquad,11) = nhexa

       if(i0.ge.iwa0.and.i2.le.iwa1)then
         print*,'iwa0,iwa1',i,j,k
         quad(nquad,10) = 0  !<------ Boundary group number
        endif
      quad(nquad,12) = 0
!    6. Top boundary face (z=1)
      elseif (k == nz1-1) then
       nquad = nquad + 1
      quad(nquad,1) = hex_loc(5)
      quad(nquad,2) = hex_loc(8)
      quad(nquad,3) = hex_loc(7)
      quad(nquad,4) = hex_loc(6)
      quad(nquad,5) = hex_loc(20)
      quad(nquad,6) = hex_loc(19)
      quad(nquad,7) = hex_loc(18)
      quad(nquad,8) = hex_loc(17)
      quad(nquad,9) = hex_loc(26)
      quad(nquad,10) = 3  !<------ Boundary group number
      quad(nquad,11) = nhexa
      quad(nquad,12) = 0
    
      endif


    end do k_loop
   end do j_loop
  end do i_loop

  write(*,*)
  write(*,*) "            Nodes = ", nnodes
  write(*,*) "           Hexas = ",  nhexa
  write(*,*) "   Quadrilaterals = ", nquad


!*******************************************************************************
! Find the periodic boundary face host.
!*******************************************************************************
 close(3)

 print*,'good1'
    print*,x2(1:3,27)

 print*,'good'
!
!...
!
  inpt = 0
 
!*******************************************************************************
! Write a grid file for rdgflo.
!******************************************************************************

  nnode = nnodes
  do i=1,nquad
   bfact(i,1)=quad(i,1)
   bfact(i,2)=quad(i,4)
   bfact(i,3)=quad(i,3)
   bfact(i,4)=quad(i,2)

   bfact(i,5)=quad(i,8)
   bfact(i,6)=quad(i,7)
   bfact(i,7)=quad(i,6)
   bfact(i,8)=quad(i,5)
   bfact(i,9)=quad(i,9)
  
   if(i==1721)then
    print*,'i==3567'
    print*,quad(i,1:9)
   endif
  enddo

   open(10,file='Taylorhexa.dat')
   write(10,*)'npoin,ntetr,npyra,npris,nhexa,ntria,nquad,time'
   write(10,'(7i9,e12.5e2)') nnode,0,0,0,nhexa,0,nquad,0.0  

   write(10,*)'element connectivity'

  do i=1,nhexa
!    write(10,'(i9,27i8)')i,hex(i,1:27)
    write(10,'(i9,27i8)')i,hex(i,1:27)
  enddo

    write(10,*)' face connectivity'
  do i=1,nquad

    write(10,'(15i8)')i,quad(i,10),quad(i,12),0,0,quad(i,11),bfact(i,1:9) ! 10 is periodic boundary
!    write(10,'(15i8)')i,quad(i,10),quad(i,12),0,0,quad(i,11),bfact(i,1:4),0,0,0,0,0 ! 10 is periodic boundary
 enddo

write(10,*)' coordinates of the nodes'

  do i=1,nnode
    write(10,'(i8,3e20.10e2)') i,(xyz(i,j),j=1,3)
  end do
  

  close(10)


  inn(1:9) = (/1,2,3,4,5,6,7,8,9 /)
  inn(10:18) = (/12,13,10,14,11,15,16,17,20 /)
  inn(19:20) = (/18,19 /)
  inn(21:27) = (/21,22,25,23,24,26,27/)


  open(10,file='bump3d.msh')
  write(10,'(11a)')'$MeshFormat'
  write(10,'(7a)') '2.2 0 8'
  write(10,'(14a)')'$EndMeshFormat'
  write(10,'(6a)')'$Nodes'
 
  write(10,*) nnode

  do i=1,nnode
   write(10,'(i8,3e20.10e2)')i,xyz(i,1:3)
  enddo

  write(10,'(9a)')'$EndNodes'
  write(10,'(9a)')'$Elements'


  write(10,*)  nquad+nhexa
 do i=1,nquad
  write(10,'(13i8)')i,10,2,17,quad(i,10)+1, bfact(i,1:9)
 enddo


 do i=1,nhexa
  write(10,'(32i8)')i+nquad,12,2,21,2, hex(i,inn(1:27))
 enddo
 
  write(10,'(12a)')'$EndElements'
  close(10)


 stop

 end program mixgrid_cube