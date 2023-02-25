      subroutine sndofd_u(nr,n,xpls,fcn,fpls,a,sx,rnoise,stepsz,anbr)

c*********************************************************************72
c
cc SNDOFD finds a second order approximation to the Hessian.
c
c  Discussion:
c
c    find second order forward finite difference approximation "a"
c    to the second derivative (hessian) of the function defined by the subp
c    "fcn" evaluated at the new iterate "xpls"
c
c    for optimization use this routine to estimate
c    1) the second derivative (hessian) of the optimization function
c    if no analytical user function has been supplied for either
c    the gradient or the hessian and if the optimization function
c    "fcn" is inexpensive to evaluate.
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
c xpls(n)      --> new iterate:   x(k)
c fcn          --> name of subroutine to evaluate function
c fpls         --> function value at new iterate, f(xpls)
c a(n,n)      <--  finite difference approximation to hessian
c                  only lower triangular matrix and diagonal
c                  are returned
c sx(n)        --> diagonal scaling matrix for x
c rnoise       --> relative noise in fname (f(x))
c stepsz(n)    --> workspace (stepsize in i-th component direction)
c anbr(n)      --> workspace (neighbor in i-th direction)
c
      implicit double precision ( a-h, o-z )

      dimension xpls(n)
      dimension sx(n)
      dimension stepsz(n),anbr(n)
      dimension a(nr,1)
      external fcn
c
c  find i-th stepsize and evaluate neighbor in direction
c  of i-th unit vector.
c
      ov3 = 1.00D+00 / 3.00D+00
      do i=1,n
        stepsz(i)=rnoise**ov3 * max(abs(xpls(i)), 1.0D+00 /sx(i))
        xtmpi=xpls(i)
        xpls(i)=xtmpi+stepsz(i)
        call fcn(n,xpls,anbr(i))
        xpls(i)=xtmpi
      end do
c
c  calculate column i of a
c
      do i=1,n

        xtmpi=xpls(i)
        xpls(i)=xtmpi+ 2.0D+00 * stepsz(i)
        call fcn(n,xpls,fhat)
        a(i,i)=((fpls-anbr(i))+(fhat-anbr(i)))/(stepsz(i)*stepsz(i))
c
c  calculate sub-diagonal elements of column
c
        if(i.eq.n) go to 30
        xpls(i)=xtmpi+stepsz(i)
        ip1=i+1
        do j=ip1,n
          xtmpj=xpls(j)
          xpls(j)=xtmpj+stepsz(j)
          call fcn(n,xpls,fhat)
          a(j,i)=((fpls-anbr(i))+(fhat-anbr(j)))/(stepsz(i)*stepsz(j))
          xpls(j)=xtmpj
        end do
   30   xpls(i)=xtmpi

      end do

      return
      end
