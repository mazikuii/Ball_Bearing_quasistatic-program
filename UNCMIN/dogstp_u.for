      subroutine dogstp_u(nr,n,g,a,p,sx,rnwtln,dlt,nwtake,fstdog,
     &     ssd,v,cln,eta,sc,iprint,lounit,stepmx)

c*********************************************************************72
c
cc DOGSTP finds the next double dogleg stepsize.
c
c  Discussion:
c
c    find new step by double dogleg algorithm
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c g(n)         --> gradient at current iterate, g(x)
c a(n,n)       --> cholesky decomposition of hessian in
c                  lower part and diagonal
c p(n)         --> newton step
c sx(n)        --> diagonal scaling matrix for x
c rnwtln       --> newton step length
c dlt         <--> trust region radius
c nwtake      <--> boolean, =.true. if newton step taken
c fstdog      <--> boolean, =.true. if on first leg of dogleg
c ssd(n)      <--> workspace (cauchy step to the minimum of the
c                  quadratic model in the scaled steepest descent
c                  direction) (retain value between successive calls)
c v(n)        <--> workspace  (retain value between successive calls)
c cln         <--> cauchy length
c                  (retain value between successive calls)
c eta              (retain value between successive calls)
c sc(n)       <--  current step
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                   0=no output.
c stepmx       --> maximum allowable step size
c
c internal variables
c ------------------
c cln              length of cauchy step
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      dimension g(n),p(n)
      dimension sx(n)
      dimension sc(n),ssd(n),v(n)
      dimension a(nr,1)
      logical nwtake,fstdog
c
c can we take newton step
c
      if(rnwtln.gt.dlt) go to 20
c     if(rnwtln.le.dlt)
c     then
        nwtake=.true.
        do i=1,n
          sc(i)=p(i)
        end do
        dlt=rnwtln
        if(iprint.gt. 1 )write(lounit,1000)
        go to 140
c     else
c
c newton step too long
c cauchy step is on double dogleg curve
c
   20   nwtake=.false.
        if(.not.fstdog) go to 80
c       if(fstdog)
c       then
c
c         calculate double dogleg curve (ssd)
          fstdog=.false.

          alpha = 0.0D+00
          do i=1,n
            alpha=alpha + (g(i)*g(i))/(sx(i)*sx(i))
          end do

          beta = 0.0D+00
          do i=1,n
            tmp = 0.0D+00
            do j=i,n
              tmp=tmp + (a(j,i)*g(j))/(sx(j)*sx(j))
            end do
            beta=beta+tmp*tmp
          end do

          do i=1,n
            ssd(i)=-(alpha/beta)*g(i)/sx(i)
          end do

          cln=alpha*sqrt(alpha)/beta

          eta= 0.2D+00 + ( 0.8D+00 *alpha*alpha)
     &      /(-beta * r8vec_dot ( n, g, p ) )

          do i=1,n
            v(i)=eta*sx(i)*p(i) - ssd(i)
          end do

          if (dlt .eq. ( -1.00D+00 )) dlt = min(cln, stepmx)
          if(iprint.gt. 1)write(lounit,1030) alpha,beta,cln,eta
          if(iprint.gt. 1)write(lounit,1040)
          if(iprint.gt. 1)write(lounit,1090) (ssd(i),i=1,n)
          if(iprint.gt. 1)write(lounit,1050)
          if(iprint.gt. 1)write(lounit,1090) (v(i),i=1,n)
c       end if
   80   if(eta*rnwtln.gt.dlt) go to 100
c       if(eta*rnwtln .le. dlt)
c       then
c
c  take partial step in newton direction
c
          do i=1,n
            sc(i)=(dlt/rnwtln)*p(i)
          end do
          if(iprint.gt. 1 )write(lounit,1060)
          go to 140
c       else
  100     if(cln.lt.dlt) go to 120
c
c  if(cln.ge.dlt) then
c  take step in steepest descent direction
c
            do i=1,n
              sc(i)=(dlt/cln)*ssd(i)/sx(i)
            end do
            if(iprint.gt. 1 )write(lounit,1070)
            go to 140
c         else
c
c           calculate convex combination of ssd and eta*p
c           which has scaled length dlt
c
  120       dot1 = r8vec_dot ( n, v, ssd )
            dot2 = r8vec_dot ( n, v, v )
            alam=(-dot1+sqrt((dot1*dot1)-dot2*(cln*cln-dlt*dlt)))/dot2
            do i=1,n
              sc(i)=(ssd(i) + alam*v(i))/sx(i)
            end do
            if(iprint.gt. 1 )write(lounit,1080)
c         end if
c       end if
c     end if
  140 continue
      if(iprint.gt. 1)write(lounit,1010) fstdog,nwtake,rnwtln,dlt
      if(iprint.gt. 1)write(lounit,1020)
      if(iprint.gt. 1)write(lounit,1090) (sc(i),i=1,n)
      return
 1000 format(27h0dogstp    take newton step)
 1010 format(18h dogstp    fstdog=,l1/
     &       18h dogstp    nwtake=,l1/
     &       18h dogstp    rnwtln=,g16.8/
     &       18h dogstp    dlt   =,g16.8)
 1020 format(28h dogstp    current step (sc))
 1030 format(18h dogstp    alpha =,g16.8/
     &       18h dogstp    beta  =,g16.8/
     &       18h dogstp    cln   =,g16.8/
     &       18h dogstp    eta   =,g16.8)
 1040 format(28h dogstp    cauchy step (ssd))
 1050 format(12h dogstp    v)
 1060 format(48h0dogstp    take partial step in newton direction)
 1070 format(50h0dogstp    take step in steepest descent direction)
 1080 format(39h0dogstp    take convex combination step)
 1090 format(14h dogstp       ,5(g16.8,3x))
      end
