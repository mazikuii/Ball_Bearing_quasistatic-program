      subroutine qraux2_u(nr,n,r,i,a,b)

c*********************************************************************72
c
cc QRAUX2 premultiplies a matrix by a Jacobi rotation.
c
c  Discussion:
c
c    pre-multiply r by the jacobi rotation j(i,i+1,a,b)
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
c n            --> dimension of matrix
c r(n,n)      <--> upper hessenberg matrix
c i            --> index of row
c a            --> scalar
c b            --> scalar
c
      implicit double precision ( a-h, o-z )

      dimension r(nr,1)

      den=sqrt(a*a + b*b)
      c=a/den
      s=b/den
      do j=i,n
        y=r(i,j)
        z=r(i+1,j)
        r(i,j)=c*y - s*z
        r(i+1,j)=s*y + c*z
      end do

      return
      end
 