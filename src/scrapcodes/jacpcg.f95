!------------------------------------------------------------
!----- Jacobi preconditioned conjugate gradient method for phi-Poisson equation: cG(P1)
subroutine jpcg_cgp1(n_iter,btype,fcprd,fsup1,fsup2,intfac,upper,diago,rhspo,eps,n_it,resnorm,phi)
integer,intent(in):: n_iter,fsup1(:),fsup2(npoin+1),intfac(5,nafac),btype(nbface),fcprd(nbface)
integer::k,i,j,ifc,npj,n_it
real*8,intent(in) ::upper(nafac),diago(npoinp2),eps
real*8 ::auu,resnorm,resnorm0,anu,adu,alph,bnu,bdu,beta
real*8 ::phi(npoinp2),phi2(npoinp2),r0(npoinp2),r1(npoinp2), &
         p0(npoinp2),z0(npoinp2),z1(npoinp2),q(npoinp2), &
         rhspo(npoinp2)

!-- initial residual r0
resnorm0 = 0.d0
do i = 1,npoinp2
   auu = 0.d0
   do j = fsup2(i)+1,fsup2(i+1)
      ifc = fsup1(j)
      if (i.eq.intfac(4,ifc)) then
         npj = intfac(3,ifc)
      else
         npj = intfac(4,ifc)
      end if
      auu = auu + (upper(ifc)*phi(npj))
   end do !j
   auu = auu + diago(i)*phi(i)
   r0(i) = rhspo(i) - auu
   resnorm0 = resnorm0 + (r0(i)*r0(i))
end do !i
resnorm0 = dsqrt(resnorm0)

!-- solution of Mz0 = r0; and p0 = r0
do i = 1,npoinp2
   z0(i) = r0(i)/diago(i)
   p0(i) = z0(i)
end do !i

!-- main iteration loop
do k = 1,n_iter

!--- BC's for Periodic conditions on phi:
   if (iperi.gt.0) then
      call peri_loop(intfac,btype,fcprd,diago,phi,rhspo)
   end if

   resnorm = 0.d0
!-- calculation of alpha_k
   !- numerator of alpha_k
   anu = 0.d0
   do i = 1,npoinp2
      anu = anu + r0(i)*z0(i)

   !- denominator of alpha_k
      q(i) = 0.d0 ! q = A*p0
      do j = fsup2(i)+1,fsup2(i+1)
         ifc = fsup1(j)
         if (i.eq.intfac(4,ifc)) then
            npj = intfac(3,ifc)
         else
            npj = intfac(4,ifc)
         end if
         q(i) = q(i) + (upper(ifc)*p0(npj))
      end do !j
      q(i) = q(i) + (diago(i)*p0(i))
   end do !i
   adu = 0.d0
   do i = 1,npoinp2
      adu = adu + p0(i)*q(i)
   end do !i
   alph = anu/adu

!-- computation of (k+1)^th iterate
   do i = 1,npoinp2
      phi2(i) = phi(i) + alph*p0(i)

!-- computation of residual
      r1(i) = r0(i) - alph*q(i)

!-- solution of Mz1 = r1
      z1(i) = r1(i)/diago(i)
   end do !i

!-- calculation of beta_k
   !- numerator of beta_k
   bnu = 0.d0
   do i = 1,npoinp2
      bnu = bnu + r1(i)*z1(i)
   end do !i
   !- denominator of beta_k
   bdu = anu
   beta = bnu/bdu

!-- modify search direction p_k
   do i = 1,npoinp2
      p0(i) = z1(i) + beta*p0(i)

!-- update phi,z and r
      phi(i) = phi2(i)
      r0(i) = r1(i)
      z0(i) = z1(i)
   end do !i

!-- residual norm calculation
   do i = 1,npoinp2
      resnorm = resnorm + (r1(i)*r1(i)) 
   end do
   resnorm = dsqrt(resnorm)
!   resnorm = resnorm/resnorm0
   n_it = k
!   print*,k,resnorm
   if (resnorm .LT. eps) exit
end do !k

end subroutine jpcg_cgp1



