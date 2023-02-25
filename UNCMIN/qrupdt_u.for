      subroutine qrupdt_u(nr,n,a,u,v)

c*********************************************************************72
c
cc QRUPDT updates a QR factorization.
c
c  Discussion:
c
c    find an orthogonal (n*n) matrix (q*) and an upper triangular (n*n)
c    matrix (r*) such that (q*)(r*)=r+u(v+)
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
c a(n,n)      <--> on input:  contains r
c                  on output: contains (r*)
c u(n)         --> vector
c v(n)         --> vector
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)
      dimension u(n),v(n)
c
c determine last non-zero in u(.)
c
      k=n
   10 continue

      if(u(k).ne. 0.0D+00 .or. k.eq. 1) go to 20
c     if(u(k).eq. 0.0D+00 .and. k.gt. 1)
c     then
        k=k-1
        go to 10
c     end if
c
c (k-1) jacobi rotations transform
c     r + u(v+) --> (r*) + (u(1)*e1)(v+)
c which is upper hessenberg
c
   20 if(k.le. 1) go to 50
        km1=k-1
        do 40 ii=1,km1
          i=km1-ii+1
          if(u(i).ne. 0.0D+00 ) go to 30
            call qraux1_u(nr,n,a,i)
            u(i)=u(i+1)
            go to 40
   30       call qraux2_u(nr,n,a,i,u(i),-u(i+1))
            u(i)=sqrt(u(i)*u(i) + u(i+1)*u(i+1))
   40   continue
c     end if
c
c  r <-- r + (u(1)*e1)(v+)
c
   50 do j=1,n
        a(1,j)=a(1,j) +u(1)*v(j)
      end do
c
c  (k-1) jacobi rotations transform upper hessenberg r
c  to upper triangular (r*)
c
      if(k.le. 1) go to 90
        km1=k-1
        do 80 i=1,km1
          if(a(i,i).ne. 0.0D+00 ) go to 70
            icopy=i
            call qraux1_u(nr,n,a,icopy)
            go to 80
   70       t1=a(i,i)
            t2=-a(i+1,i)
            icopy=i
            call qraux2_u(nr,n,a,icopy,t1,t2)
   80   continue
   90 return
      end
