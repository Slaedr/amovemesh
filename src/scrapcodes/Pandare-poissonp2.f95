! read P2-mesh rdgflo format---------------------------------------------------------
subroutine readfilep2(filename,btype,inpoeltri,inpoelfc,coord)
!-------------------------------------------------------------------------------------
integer                                 ::cntr,cntr2,ndimn,dump,i,n
integer,dimension(:),allocatable        ::btype
integer,dimension(:,:),allocatable      ::inpoeltri,inpoelfc
real*8::time
real*8,dimension(:,:),allocatable::coord
character(len=50),intent(in)            ::filename

open(11,file='./meshes/'//filename,action='read')
read(11,'(1I6)') n
do i = 1,n
   read(11,*) !dump descriptive lines
end do
read(11,*)
read(11,'(2I6)') ndimn,dump
read(11,*)
read(11,*) nelem,npoin,npoinp2,nbface,time
read(11,*) ! dumping the title

!----- reading connectivities -------------------------------------------------------!

allocate(inpoeltri(nnodep2,nelem))
do i=1,nelem
   read(11,*)dump,inpoeltri(1,i),inpoeltri(2,i),inpoeltri(3,i), &
             inpoeltri(4,i),inpoeltri(5,i),inpoeltri(6,i)
end do

read(11,*) ! dumping the title

!----- reading coordinates of grid-points -------------------------------------------!

allocate(coord(2,npoinp2))
do i=1,npoinp2
   read(11,*)dump,coord(1,i),coord(2,i)
end do

read(11,*) ! dumping the title

do i=1,npoinp2
   read(11,*) ! dumping the "initial value" data
end do

read(11,*) ! dumping the title

!----- reading boundary/face-cell connectivities ------------------------------------!

allocate(btype(nbface),inpoelfc(2,nbface))
do i=1,nbface
   read(11,*)dump,inpoelfc(1,i),inpoelfc(2,i),btype(i)
end do

close(11)

print*,"Mesh stats.:"
print*,"-------------------------------------------------------"
write(*,'(T2,4A10)')"npoin:","npoinp2:","nbface:","nelem:"
write(*,'(T2,4I10)')npoin,npoinp2,nbface,nelem
print*,"-------------------------------------------------------"

end subroutine readfilep2
!--------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
!----- lhs matrix generation for phi: computes the laplacian operator:
!-------------------------------------------------------------------------------------------

subroutine laplacev(geoel,a)
real*8 ,intent(in)::geoel(5)
integer  ::i,j,p,q
real*8   ::dbdx(3),dbdy(3),jacd,a(6,6)

a = 0.d0
!--- derivatives of lambda(1,2,3):
   dbdx(1) = geoel(1)
   dbdx(2) = geoel(2)
   dbdx(3) = -(dbdx(1) + dbdx(2))
   dbdy(1) = geoel(3)
   dbdy(2) = geoel(4)
   dbdy(3) = -(dbdy(1) + dbdy(2))
   jacd = geoel(5)

!--- A(i,i) for all i<4:
   a(1,1) = (dbdx(1)*dbdx(1) + dbdy(1)*dbdy(1))*jacd*0.5d0
   a(2,2) = (dbdx(2)*dbdx(2) + dbdy(2)*dbdy(2))*jacd*0.5d0
   a(3,3) = (dbdx(3)*dbdx(3) + dbdy(3)*dbdy(3))*jacd*0.5d0

!--- A(i,j) for all i,j<4 .and. i/=j:
   do i = 1,3
   do j = 1,3
      if (i.ne.j) then
         a(i,j) = -(dbdx(i)*dbdx(j) + dbdy(i)*dbdy(j))*jacd/6.d0
      end if
   end do !j
   end do !i

!--- A(i,i) for all i>3:
   a(4,4) = (dbdx(1)*dbdx(1) + dbdy(1)*dbdy(1) &
           + dbdx(2)*dbdx(2) + dbdy(2)*dbdy(2) &
           + dbdx(1)*dbdx(2) + dbdy(1)*dbdy(2)) * (4.d0*jacd/3.d0)

   a(5,5) = (dbdx(2)*dbdx(2) + dbdy(2)*dbdy(2) &
           + dbdx(3)*dbdx(3) + dbdy(3)*dbdy(3) &
           + dbdx(2)*dbdx(3) + dbdy(2)*dbdy(3)) * (4.d0*jacd/3.d0)

   a(6,6) = (dbdx(3)*dbdx(3) + dbdy(3)*dbdy(3) &
           + dbdx(1)*dbdx(1) + dbdy(1)*dbdy(1) &
           + dbdx(3)*dbdx(1) + dbdy(3)*dbdy(1)) * (4.d0*jacd/3.d0)

!--- A(i,j) for all i,j>3 .and. i/=j:
   a(4,5) = (dbdx(1)*dbdx(2) + dbdy(1)*dbdy(2) &
           + dbdx(2)*dbdx(2) + dbdy(2)*dbdy(2) &
           + dbdx(2)*dbdx(3) + dbdy(2)*dbdy(3)) * (2.d0*jacd/3.d0) &
           +(dbdx(1)*dbdx(3) + dbdy(1)*dbdy(3)) * (4.d0*jacd/3.d0)

   a(4,6) = (dbdx(1)*dbdx(3) + dbdy(1)*dbdy(3) &
           + dbdx(1)*dbdx(1) + dbdy(1)*dbdy(1) &
           + dbdx(2)*dbdx(1) + dbdy(2)*dbdy(1)) * (2.d0*jacd/3.d0) &
           +(dbdx(2)*dbdx(3) + dbdy(2)*dbdy(3)) * (4.d0*jacd/3.d0)

   a(5,6) = (dbdx(2)*dbdx(3) + dbdy(2)*dbdy(3) &
           + dbdx(3)*dbdx(3) + dbdy(3)*dbdy(3) &
           + dbdx(3)*dbdx(1) + dbdy(3)*dbdy(1)) * (2.d0*jacd/3.d0) &
           +(dbdx(2)*dbdx(1) + dbdy(2)*dbdy(1)) * (4.d0*jacd/3.d0)

   a(5,4) = a(4,5)
   a(6,4) = a(4,6)
   a(6,5) = a(5,6)

!--- A(i,j) for all i<4 .and. j>3 .and. i/=j  i.e. the rest of the terms:

   a(1,4) = (dbdx(1)*dbdx(2) + dbdy(1)*dbdy(2))*2.d0*jacd/3.d0
   a(2,5) = (dbdx(2)*dbdx(3) + dbdy(2)*dbdy(3))*2.d0*jacd/3.d0
   a(3,6) = (dbdx(3)*dbdx(1) + dbdy(3)*dbdy(1))*2.d0*jacd/3.d0
   a(1,5) = 0.d0
   a(1,6) = (dbdx(1)*dbdx(3) + dbdy(1)*dbdy(3))*2.d0*jacd/3.d0
   a(2,4) = (dbdx(2)*dbdx(1) + dbdy(2)*dbdy(1))*2.d0*jacd/3.d0
   a(2,6) = 0.d0
   a(3,4) = 0.d0
   a(3,5) = (dbdx(3)*dbdx(2) + dbdy(3)*dbdy(2))*2.d0*jacd/3.d0

   a(4,1) = a(1,4)
   a(5,2) = a(2,5)
   a(6,3) = a(3,6)
   a(5,1) = a(1,5)
   a(6,1) = a(1,6)
   a(4,2) = a(2,4)
   a(6,2) = a(2,6)
   a(4,3) = a(3,4)
   a(5,3) = a(3,5)

!do i = 1,6
!   print*,a(i,1),a(i,2),a(i,3),a(i,4),a(i,5),a(i,6)
!end do
!       stop

end subroutine laplacev

!--------------------------------------------------------------------------------------

subroutine intfac_gen(inpoel,inpoelfc,esuel,intfac)
integer::inpoel(nnodep2,nelem),esuel(nfael,nelem),ielem,ifael,ifc,icoun, &
         lpofa(2,3),inpoelfc(2,nbface)
integer,allocatable::intfac(:,:)
allocate(intfac(5,nafac))
intfac(:,:) = 0
lpofa(1,1) = 2
lpofa(2,1) = 3
lpofa(1,2) = 3
lpofa(2,2) = 1
lpofa(1,3) = 1
lpofa(2,3) = 2
icoun = 1
do ielem = 1,nelem
   do ifael = 1,nfael
      if(ielem.lt.esuel(ifael,ielem)) then
        intfac(1,nbface+icoun) = ielem
        intfac(2,nbface+icoun) = esuel(ifael,ielem)
        intfac(3:4,nbface+icoun) = inpoel(lpofa(1:2,ifael),ielem)
        icoun = icoun +1
      end if
   end do
end do

do ifc = 1,nbface
   intfac(3:4,ifc) = inpoelfc(1:2,ifc)
end do

do ifc = 1,nbface
    do ielem = 1,nelem
        do ifael = 1,nfael
            if(esuel(ifael,ielem).eq.0) then
            if(((inpoel(1,ielem).eq.intfac(3,ifc)) .or. &
                (inpoel(2,ielem).eq.intfac(3,ifc)) .or. &
                (inpoel(3,ielem).eq.intfac(3,ifc))) .and. &
                ((inpoel(1,ielem).eq.intfac(4,ifc)) .or. &
                (inpoel(2,ielem).eq.intfac(4,ifc)) .or. &
                (inpoel(3,ielem).eq.intfac(4,ifc)))) then
                intfac(1,ifc) = ielem
                intfac(2,ifc) = nelem + ifc
            end if
            end if
        end do
    end do
end do
!print*,'intfac'
!do ifael=1,nafac
!write(*,'(5I7)') ifael,intfac(1,ifael),intfac(2,ifael),intfac(3,ifael),intfac(4,ifael)
!end do
end subroutine intfac_gen

!-----------------------------------------------------------------------------!

!----------------- on-the-fly construction of the global stiffness ------------

   q = 0.d0
   do ie = 1,nelem
      ip1 = inpoel(1,ie)
      ip2 = inpoel(2,ie)
      ip3 = inpoel(3,ie)
      ip4 = inpoel(4,ie)
      ip5 = inpoel(5,ie)
      ip6 = inpoel(6,ie)
      !-- Matrix vector product:
      do j = 1,nnodep2
         jp = inpoel(j,ie)
         q(jp) = q(jp) + (a(j,1,ie)*p0(ip1) &
                         +a(j,2,ie)*p0(ip2) &
                         +a(j,3,ie)*p0(ip3) &
                         +a(j,4,ie)*p0(ip4) &
                         +a(j,5,ie)*p0(ip5) &
                         +a(j,6,ie)*p0(ip6))
      end do !j
   end do !ie

!----------------------------------------------------------------------------------
