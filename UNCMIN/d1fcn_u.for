      subroutine d1fcn_u(n,x,g)

c*********************************************************************72
c
cc D1FCN is a dummy routine for the analytic gradient.
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
c dummy routine to prevent unsatisfied external diagnostic
c when specific analytic gradient function not supplied.
c
      implicit none

	integer n
      double precision :: x(n),g(n)

      write ( *, 1000 )
      write ( *, 1010 )
      stop
 1000 format(' d1fcn - dummy gradient routine called.')
 1010 format(' d1fcn - program forced to stop.')
      end
