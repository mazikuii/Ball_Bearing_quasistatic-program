      subroutine choldc_u(nr,n,a,diagmx,tol,addmax)

c*********************************************************************72
c
cc CHOLDC finds the perturbed Cholesky decomposition of A+D.
c
c  Discussion:
c
c    find the perturbed l(l-transpose), written ll+, decomposition
c    of a+d, where d is a non-negative diagonal matrix added to a if
c    necessary to allow the cholesky decomposition to continue.
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
c a(n,n)      <--> on entry: matrix for which to find perturbed
c                       cholesky decomposition
c                  on exit:  contains l of ll+ decomposition
c                  in lower triangular part and diagonal of "a"
c diagmx       --> maximum diagonal element of "a"
c tol          --> tolerance
c addmax      <--  maximum amount implicitly added to diagonal of "a"
c                  in forming the cholesky decomposition of a+d
c internal variables
c ------------------
c aminl    smallest element allowed on diagonal of l
c amnlsq   =aminl**2
c offmax   maximum off-diagonal element in column of a
c
c
c description
c -----------
c the normal cholesky decomposition is performed.  however, if at any
c point the algorithm would attempt to set l(i,i)=sqrt(temp)
c with temp < tol*diagmx, then l(i,i) is set to sqrt(tol*diagmx)
c instead.  this is equivalent to adding tol*diagmx-temp to a(i,i)
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)

      addmax = 0.0D+00
      aminl=sqrt(diagmx*tol)
      amnlsq=aminl*aminl
c
c form column j of l
c
      do 100 j=1,n
c
c find diagonal elements of l
c
        sum = 0.0D+00
        if(j.eq. 1 ) go to 20
        jm1=j-1
        do k=1,jm1
          sum=sum + a(j,k)*a(j,k)
        end do
   20   temp=a(j,j)-sum
        if(temp.lt.amnlsq) go to 30
c       if(temp.ge.aminl**2)
c       then
          a(j,j)=sqrt(temp)
          go to 60
c       else
c
c find maximum off-diagonal element in column
c
   30     offmax = 0.0D+00
          if(j.eq.n) go to 50
          jp1=j+1
          do i=jp1,n
            if(abs(a(i,j)).gt.offmax) offmax=abs(a(i,j))
          end do
   50     if(offmax.le.amnlsq) offmax=amnlsq
c
c add to diagonal element  to allow cholesky decomposition to continue
c
          a(j,j)=sqrt(offmax)
          addmax=max(addmax,offmax-temp)
c       end if
c
c find i,j element of lower triangular matrix
c
   60   if(j.eq.n) go to 100
        jp1=j+1
        do i=jp1,n
          sum = 0.0D+00
          if(j.eq. 1 ) go to 80
          jm1=j-1
          do k=1,jm1
            sum=sum+a(i,k)*a(j,k)
          end do
   80     a(i,j)=(a(i,j)-sum)/a(j,j)
        end do
  100 continue
      return
      end
