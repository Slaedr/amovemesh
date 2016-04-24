c
c-----------------------------------------------------------------------
c
c     LocateTriangle sub locates the triangle in which  the point p lies
c
c----------------------------------------------------------------------*/
c
      function findel(ndimn ,nnode, melem, mpoin, nelem ,
     &                intmat,elsuel,coord, iele0, x, y )
c
      implicit real*8 (a-h,o-z)
c
      integer  intmat(nnode,melem),elsuel(nnode,melem),iele0
      integer  findel
      real*8   coord(ndimn,mpoin), x, y
c
      real*8   shape(3)
c
      data tolow/ -0.005 /
      data tolhi/  1.005 /
      data   c00/  0.0   /
      data   c10/  1.0   /
c
      ielem = iele0
c
 999  continue 
c
      do 100 in = 1, nnode
      ip1 = intmat(in,ielem)
      ip2 = intmat(mod(in,nnode)+1,ielem)
      if( (coord(2,ip1)-y)*(coord(1,ip2)-x) .gt.
     &    (coord(1,ip1)-x)*(coord(2,ip2)-y) ) then		! (A) using cross product?
      ielem = elsuel(mod(in+1,nnode)+1,ielem)
      if(ielem .ne. 0 ) goto 999
      endif
 100  continue
      findel = ielem 
c
c...  find the host element the hard way, if we can find it
c
      if(ielem .eq. 0) then
      do ie = 1, nelem
      ip1   = intmat(1,ie)
      ip2   = intmat(2,ie)
      ip3   = intmat(3,ie)
c
      x1  = coord(1,ip1)
      y1  = coord(2,ip1)
c
      x21 = coord(1,ip2) - x1
      y21 = coord(2,ip2) - y1
      x31 = coord(1,ip3) - x1
      y31 = coord(2,ip3) - y1
c
      rj  = x21*y31 - x31*y21
      rj1 = c10/rj
c
c...  find the local coordinates
c
      xp1 = x - x1
      yp1 = y - y1
c
c...  local coordinates
c
      shape(2) = rj1*( y31*xp1 - x31*yp1)
      shape(3) = rj1*(-y21*xp1 + x21*yp1)
      shape(1) = c10 - shape(2) - shape(3)
c
c...  max/min of these shape-functions
c
      shmin = min(shape(1),shape(2),shape(3))
      shmax = max(shape(1),shape(2),shape(3))
c
c...  see if in the element
c
      if(shmin .ge. tolow .and. shmax .le. tolhi) then
      findel = ie
      return
      endif
c
      enddo
      endif
c
      return
      end
c
c------------------------------------------------------------------------*/
