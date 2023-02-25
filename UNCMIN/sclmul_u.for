      subroutine sclmul_u(n,s,v,z)

c*********************************************************************72
c
cc SCLMUL multiplies a vector by a scalar.
c
c  Discussion:
c
c    multiply vector by scalar
c    result vector may be operand vector
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
c n            --> dimension of vectors
c s            --> scalar
c v(n)         --> operand vector
c z(n)        <--  result vector
c
      implicit double precision ( a-h, o-z )

      dimension v(n),z(n)

      do i=1,n
        z(i)=s*v(i)
      end do

      return
      end
