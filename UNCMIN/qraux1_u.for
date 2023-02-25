      subroutine qraux1_u(nr,n,r,i)

c*********************************************************************72
c
cc QRAUX1 interchanges two rows of an upper Hessenberg matrix.
c
c  Discussion:
c
c    interchange rows i,i+1 of the upper hessenberg matrix r,
c    columns i to n
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
c i            --> index of row to interchange (i.lt.n)
c
      implicit double precision ( a-h, o-z )

      dimension r(nr,1)

      do j=i,n
        tmp=r(i,j)
        r(i,j)=r(i+1,j)
        r(i+1,j)=tmp
      end do

      return
      end
