      subroutine forslv_u(nr,n,a,x,b)

c*********************************************************************72
c
cc FORSLV solves A*X=B, where A is lower triangular.
c
c  Discussion:
c
c    solve  ax=b  where a is lower triangular matrix
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
c a(n,n)       --> lower triangular matrix (preserved)
c x(n)        <--  solution vector
c b(n)         --> right-hand side vector
c
c note
c ----
c if b is no longer required by calling routine,
c then vectors b and x may share the same storage.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),b(n)
c
c solve lx=b. (foreward solve)
c
      x(1)=b(1)/a(1,1)
      if(n.eq. 1) return
      do i=2,n
        sum = 0.0D+00
        im1=i-1
        do j=1,im1
          sum=sum+a(i,j)*x(j)
        end do
        x(i)=(b(i)-sum)/a(i,i)
      end do

      return
      end
