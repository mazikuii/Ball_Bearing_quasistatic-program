      subroutine hsnint_u(nr,n,a,sx,method)

c*********************************************************************72
c
cc HSNINT approximates the initial Hessian by secant updates.
c
c  Discussion:
c
c    provide initial hessian when using secant updates
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
c a(n,n)      <--  initial hessian (lower triangular matrix)
c sx(n)        --> diagonal scaling matrix for x
c method       --> algorithm to use to solve minimization problem
c                    =1,2 factored secant method used
c                    =3   unfactored secant method used
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),sx(n)

      do 20 j=1,n
        if(method.eq. 3) a(j,j)=sx(j)*sx(j)
        if(method.ne. 3) a(j,j)=sx(j)
        if(j.eq.n) go to 20
        jp1=j+1
        do i=jp1,n
          a(i,j) = 0.0D+00
        end do
   20 continue

      return
      end
