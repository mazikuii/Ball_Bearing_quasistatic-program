      subroutine dfault_u(n,typsiz,fscale,method,iexp,msg,ndigit,
     &     itnlim,iagflg,iahflg,iprint,lounit,dlt,gradtl,stepmx,steptl)

c*********************************************************************72
c
cc DFAULT sets default values for each input variable to the minimization package.
c
c  Discussion:
c
c    set default values for each input variable to
c    minimization algorithm.
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
c typsiz(n)   <--  typical size for each component of x
c                  set to 1.
c fscale      <--  estimate of scale of minimization function
c                  set to 1.
c method      <--  algorithm to use to solve minimization problem
c                  set to 1
c iexp        <--  =0 if minimization function not expensive to evaluate
c                  set to 1
c msg         <--  message to inhibit certain automatic checks + output
c                  set to 0
c ndigit      <--  number of good digits in minimization function
c                  use negative number if function assumed to have
c                  maximum possible number of good digits.
c                  set to -1
c itnlim      <--  maximum number of allowable iterations
c                  set to 150.
c iagflg      <--  =0 if analytic gradient not supplied
c                  set to 0.
c iahflg      <--  =0 if analytic hessian not supplied
c                  set to 0.
c iprint       --> amount of output desired.
c                  set to 1
c lounit      <--  device to which to send output
c                  0=no output.
c                  set to 6.
c dlt         <--  trust region radius
c                  if unknown, set to -1.
c                  set to -1.
c gradtl      <--  tolerance at which gradient considered close enough
c                  to zero to terminate algorithm
c                  set to cube root of relative machine precision.
c stepmx      <--  value of zero to trip default maximum in optchk
c                  set to 0.0
c steptl      <--  tolerance at which successive iterates considered
c                  close enough to terminate algorithm
c                  set to square root of relative machine precision.
c
      implicit double precision ( a-h, o-z )

      double precision r8_epsilon
      dimension typsiz(n)
c
c set typical size of x and minimization function
c
      do i=1,n
        typsiz(i) = 1.0D+00
      end do
      fscale = 1.0D+00
c
c set tolerances
c
      dlt = -1.0D+00

      epsm = r8_epsilon ( )

      gradtl=epsm**( 1.0D+00 / 3.0D+00 )
      stepmx = 0.0D+00
      steptl=sqrt(epsm)
c
c set flags
c
      method=1
      iexp=1
      msg=0
      ndigit=-1
      itnlim=150
      iagflg=0
      iahflg=0
      iprint=1
      lounit=6

      return
      end
