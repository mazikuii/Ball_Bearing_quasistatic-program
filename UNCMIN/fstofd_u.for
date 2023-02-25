      subroutine fstofd_u(nr,m,n,xpls,fcn,fpls,a,sx,rnoise,fhat,icase)

c*********************************************************************72
c
cc FSTOFD finds forward difference approximation to the gradient.
c
c  Discussion:
c
c    find first order forward finite difference approximation "a" to the
c    first derivative of the function defined by the subprogram "fname"
c    evaluated at the new iterate "xpls".
c
c
c    for optimization use this routine to estimate:
c    1) the first derivative (gradient) of the optimization function "fcn
c       analytic user routine has been supplied;
c    2) the second derivative (hessian) of the optimization function
c       if no analytic user routine has been supplied for the hessian but
c       one has been supplied for the gradient ("fcn") and if the
c       optimization function is inexpensive to evaluate
c
c    _m=1 (optimization) algorithm estimates the gradient of the function
c         (fcn).   fcn(x) # f: r(n)-->r(1)
c    _m=n (systems) algorithm estimates the jacobian of the function
c         fcn(x) # f: r(n)-->r(n).
c    _m=n (optimization) algorithm estimates the hessian of the optimizatio
c         function, where the hessian is the first derivative of "fcn"
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
c m            --> number of rows in a
c n            --> number of columns in a; dimension of problem
c xpls(n)      --> new iterate:  x(k)
c fcn          --> name of subroutine to evaluate function
c fpls(m)      --> _m=1 (optimization) function value at new iterate:
c                       fcn(xpls)
c                  _m=n (optimization) value of first derivative
c                       (gradient) given by user function fcn
c                  _m=n (systems)  function value of associated
c                       minimization function
c a(nr,n)     <--  finite difference approximation (see note).  only
c                  lower triangular matrix and diagonal are returned
c sx(n)        --> diagonal scaling matrix for x
c rnoise       --> relative noise in fcn (f(x))
c fhat(m)      --> workspace
c icase        --> =1 optimization (gradient)
c                  =2 systems
c                  =3 optimization (hessian)
c
c internal variables
c ------------------
c stepsz - stepsize in the j-th variable direction
c
      implicit double precision ( a-h, o-z )

      dimension xpls(n),fpls(m)
      dimension fhat(m)
      dimension sx(n)
      dimension a(nr,1)
      external fcn
c
c find j-th column of a
c each column is derivative of f(fcn) with respect to xpls(j)
c
      do j=1,n
        stepsz=sqrt(rnoise)*max(abs(xpls(j)), 1.0D+00 /sx(j))
        xtmpj=xpls(j)
        xpls(j)=xtmpj+stepsz
        call fcn(n,xpls,fhat)
        xpls(j)=xtmpj
        do i=1,m
          a(i,j)=(fhat(i)-fpls(i))/stepsz
        end do
      end do

      if(icase.ne. 3) return
c
c if computing hessian, a must be symmetric
c
      if(n.eq. 1) return
      nm1=n-1
      do j=1,nm1
        jp1=j+1
        do i=jp1,m
          a(i,j)=(a(i,j)+a(j,i)) / 2.0D+00
        end do
      end do

      return
      end
