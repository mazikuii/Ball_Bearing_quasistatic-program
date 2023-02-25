      subroutine lltslv_u(nr,n,a,x,b)

c*********************************************************************72
c
cc LLTSLV solves A*X=B when A has been Cholesky factored.
c
c  Discussion:
c
c    solve ax=b where a has the form l(l-transpose)
c    but only the lower triangular part, l, is stored.
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
c a(n,n)       --> matrix of form l(l-transpose).
c                  on return a is unchanged.
c x(n)        <--  solution vector
c b(n)         --> right-hand side vector
c
c note
c ----
c if b is not required by calling program, then
c b and x may share the same storage.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),b(n)
c
c forward solve, result in x
c
      call forslv_u(nr,n,a,x,b)
c
c back solve, result in x
c
      call bakslv_u(nr,n,a,x,x)
      return
      end
