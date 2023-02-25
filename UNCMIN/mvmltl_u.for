      subroutine mvmltl_u(nr,n,a,x,y)

c*********************************************************************72
c
cc MVMLTL computes Y = L * X, where L is lower triangular.
c
c  Discussion:
c
c    compute y=lx where l is a lower triangular matrix stored in a
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
c a(n,n)       --> lower triangular (n*n) matrix
c x(n)         --> operand vector
c y(n)        <--  result vector
c
c note
c ----
c x and y cannot share storage
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),y(n)

      do i=1,n
        sum = 0.0D+00
        do j=1,i
          sum=sum+a(i,j)*x(j)
        end do
        y(i)=sum
      end do

      return
      end
