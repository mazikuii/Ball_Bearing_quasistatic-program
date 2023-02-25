      subroutine fstocd_u(n,x,fcn,sx,rnoise,g)

c*********************************************************************72
c
cc FSTOCD finds central difference approximation to the gradient.
c
c  Discussion:
c
c    find central difference approximation g to the first derivative
c    (gradient) of the function defined by fcn at the point x.
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
c x            --> point at which gradient is to be approximated.
c fcn          --> name of subroutine to evaluate function.
c sx           --> diagonal scaling matrix for x.
c rnoise       --> relative noise in fcn (f(x)).
c g           <--  central difference approximation to gradient.
c
      implicit double precision ( a-h, o-z )

      dimension x(n)
      dimension sx(n)
      dimension g(n)
      external fcn
c
c find i th  stepsize, evaluate two neighbors in direction of i th
c unit vector, and evaluate i th  component of gradient.
c
      third = 1.0D+00 / 3.0D+00
      do i = 1, n
         stepi = rnoise**third * max(abs(x(i)), 1.0D+00 / sx(i))
         xtempi = x(i)
         x(i) = xtempi + stepi
         call fcn (n, x, fplus)
         x(i) = xtempi - stepi
         call fcn (n, x, fminus)
         x(i) = xtempi
         g(i) = (fplus - fminus) / ( 2.0D+00 * stepi)
      end do

      return
      end
