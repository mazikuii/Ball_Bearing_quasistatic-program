      subroutine optchk_u(n,x,typsiz,sx,fscale,gradtl,itnlim,ndigit
     &			,epsm,
     &     dlt,method,iexp,iagflg,iahflg,stepmx,msg,iprint,lounit)

c*********************************************************************72
c
cc OPTCHK checks the input for reasonableness.
c
c  Discussion:
c
c    check input for reasonableness
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
c n            --> dimension of problem
c x(n)         --> on entry, estimate to root of fcn
c typsiz(n)   <--> typical size of each component of x
c sx(n)       <--  diagonal scaling matrix for x
c fscale      <--> estimate of scale of objective function fcn
c gradtl       --> tolerance at which gradient considered close
c                  enough to zero to terminate algorithm
c itnlim      <--> maximum number of allowable iterations
c ndigit      <--> number of good digits in optimization function fcn
c epsm         --> machine epsilon
c dlt         <--> trust region radius
c method      <--> algorithm indicator
c iexp        <--> expense flag
c iagflg      <--> =1 if analytic gradient supplied
c iahflg      <--> =1 if analytic hessian supplied
c stepmx      <--> maximum step size
c msg         <--> message and error code
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c
      implicit double precision ( a-h, o-z )

      dimension x(n),typsiz(n),sx(n)
c
c  check output unit first
c
      if(lounit.le.0)iprint=0
c
c  check that parameters only take on acceptable values.
c  if not, set them to default values.
c
      if(method.lt. 1 .or. method.gt. 3) method=1
      if(iagflg.ne. 1 ) iagflg=0
      if(iahflg.ne. 1 ) iahflg=0
      if(iexp.ne.0 ) iexp=1
c
c  check dimension of problem
c
      if(n.le.0) go to 40
c
c  compute scale matrix
c
      do i=1,n
        if(typsiz(i).eq. 0.0D+00 ) typsiz(i) = 1.0D+00
        if(typsiz(i).lt. 0.0D+00 ) typsiz(i)=-typsiz(i)
        sx(i)=1.0D+00 / typsiz(i)
      end do
c
c  check maximum step size
c
      if ( stepmx .le. 0.0D+00 ) then

        stpsiz = 0.0D+00
        do i = 1, n
          stpsiz = stpsiz + x(i)*x(i)*sx(i)*sx(i)
        end do
        stpsiz = sqrt(stpsiz)
        stepmx = max ( 1.0D+03 * stpsiz, 1.0D+03 )

      end if
c
c  check function scale
c
      if ( fscale .eq. 0.0D+00 ) then
        fscale = 1.0D+00
      else if ( fscale .lt. 0.0D+00 ) then
        fscale = - fscale
      end if
c
c  check gradient tolerance
c
      if(gradtl.lt. 0.0D+00 ) go to 50
c
c  check iteration limit
c
      if(itnlim.le.0) go to 60
c
c  check number of digits of accuracy in function fcn
c
      if(ndigit.eq.0) go to 70
      if(ndigit.lt.0) ndigit=-log10(epsm)
c
c  check trust region radius
c
      if(dlt.le. 0.0D+00 ) dlt = -1.0D+00
      if (dlt .gt. stepmx) dlt = stepmx
      return
c
c  error exits
c
   40 if(iprint.gt.0)write(lounit,1000) n
      msg=-1
      go to 80
   50 if(iprint.gt.0)write(lounit,1010) gradtl
      msg=-3
      go to 80
   60 if(iprint.gt.0)write(lounit,1020) itnlim
      msg=-4
      go to 80
   70 if(iprint.gt.0)write(lounit,1030) ndigit
      msg=-5
      go to 80
   80 return

 1000 format(32h0optchk    illegal dimension, n=,i5)
 1010 format(38h0optchk    illegal tolerance.  gradtl=,g16.8)
 1020 format(44h0optchk    illegal iteration limit.  itnlim=,i5)
 1030 format(52h0optchk    minimization function has no good digits.,
     &        9h  ndigit=,i5)
      end
