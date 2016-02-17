c
c====================================================================
c...  Unit of subs for direct solver
c====================================================================
c
c...  LU decomposition stocked in 1D profile
c
      subroutine lupoin(npoin, npoi1 ,nstor, upper, lower, ldiag )
c
c		Aditya : ldiag is an array of size nrows+1. ldiag(0) = 0, and ldiag(i+1) = ldiag(i) + bandwidth(row i)
c
      implicit real*8 (a-h,o-z)
c
      integer    ldiag(npoi1)
      real*8     upper(nstor),lower(nstor)
      do 10 i = 1, npoin
		lbeg    = ldiag(i) + 1
		lend    = ldiag(i+1)
		lband   = lend - lbeg + 1
		mi      = i + 1 - lband
	c
	c...  decompose the non diagonal terms
	c
		istor   = ldiag(i+1) - i
		do 20 j = mi, i-1 
			if(j .eq. 1) then
				mj      = 1
				else
				mj      = j + 1 - ldiag(j+1) + ldiag(j)
			endif
			kij     = max(mi, mj )
			ij      = istor + j
			jj      = ldiag(j+1)
			do 30 k = kij, j-1 
				ik      = istor + k
				jk      = ldiag(j+1) - j + k
				lower(ij) = lower(ij) - lower(ik)*upper(jk)
				upper(ij) = upper(ij) - lower(jk)*upper(ik)
			30   continue
			lower(ij) = lower(ij)*upper(jj)
		20   continue
	c
	c...  diagonal terms
	c
		ii      = ldiag(i+1)
		do 40 k = mi, i-1
			ik      = istor + k
			upper(ii) = upper(ii) - lower(ik)*upper(ik)
		40   continue
	c
	c...  inverse diagonal terms
	c
		upper(ii) = 1.0/upper(ii)
	c    
 10   continue
      return
      end
c
c...   LU solution stocked in 1D profile
c
      subroutine repoin(npoin, npoi1, nstor, upper, lower, ldiag ,
     &                  rhs  )
c
      implicit real*8 (a-h,o-z)
c
      integer    ldiag(npoi1)
      real*8     upper(nstor),lower(nstor)
      real*8     rhs(npoin)
c
c...  backword elemination
c
      do 10 i = 1, npoin
      lbeg    = ldiag(i) + 1
      lend    = ldiag(i+1)
      lband   = lend - lbeg + 1
      mi      = i + 1 - lband
      istor   = ldiag(i+1) - i
      do 20 j = mi, i-1 
      ij      = istor + j
      rhs(i)  = rhs(i) - lower(ij)*rhs(j)
 20   continue
 10   continue
c
c...  forward substitution
c
      do 30 i = npoin, 1, -1
      if(i .eq. 1) then
      mi      = 1
      else
      mi      = i + 1 - ldiag(i+1) + ldiag(i)
      endif
      istor   = ldiag(i+1) - i
      ii      = ldiag(i+1)
      rhs(i)  = upper(ii)*rhs(i)
      do 40 j = mi, i-1
      ij      = istor + j 
      rhs(j)  = rhs(j) - upper(ij)*rhs(i)
 40   continue
 30   continue
      return
      end
c
c-------------------------------------------------------------
c
