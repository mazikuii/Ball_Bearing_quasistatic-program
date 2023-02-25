      subroutine chlhsn_u(nr,n,a,epsm,sx,udiag)

c*********************************************************************72
c
cc CHLHSN Cholesky decomposes the perturbed model Hessian matrix.
c
c  Discussion:
c
c    find the l(l-transpose), written ll+, decomposition of the perturbed
c    model hessian matrix a+mu*i(where mu.ne.0 and i is identity matrix)
c    which is safely positive definite.  if a is safely positive definite
c    upon entry, then mu=0.
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
c a(n,n)      <--> on entry; "a" is model hessian (only lower
c                  triangular part and diagonal stored)
c                  on exit:  a contains l of ll+ decomposition of
c                  perturbed model hessian in lower triangular
c                  part and diagonal and contains hessian in upper
c                  triangular part and udiag
c epsm         --> machine epsilon
c sx(n)        --> diagonal scaling matrix for x
c udiag(n)    <--  on exit: contains diagonal of hessian
c
c internal variables
c ------------------
c tol              tolerance
c diagmn           minimum element on diagonal of a
c diagmx           maximum element on diagonal of a
c offmax           maximum off-diagonal element of a
c offrow           sum of off-diagonal elements in a row of a
c evmin            minimum eigenvalue of a
c evmax            maximum eigenvalue of a
c
c description
c -----------
c 1. if "a" has any negative diagonal elements, then choose mu>0
c such that the diagonal of a:=a+mu*i is all positive
c with the ratio of its smallest to largest element on the
c order of sqrt(epsm).
c
c 2. "a" undergoes a perturbed cholesky decomposition which
c results in an ll+ decomposition of a+d, where d is a
c non-negative diagonal matrix which is implicitly added to
c "a" during the decomposition if "a" is not positive definite.
c "a" is retained and not changed during this process by
c copying l into the upper triangular part of "a" and the
c diagonal into udiag.  then the cholesky decomposition routine
c is called.  on return, addmax contains maximum element of d.
c
c 3. if addmax=0, "a" was positive definite going into step 2
c and return is made to calling program.  otherwise,
c the minimum number sdd which must be added to the
c diagonal of a to make it safely strictly diagonally dominant
c is calculated.  since a+addmax*i and a+sdd*i are safely
c positive definite, choose mu=min(addmax,sdd) and decompose
c a+mu*i to obtain l.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),sx(n),udiag(n)
c
c scale hessian
c pre- and post- multiply "a" by inv(sx)
c
      do j=1,n
        do i=j,n
          a(i,j)=a(i,j)/(sx(i)*sx(j))
        end do
      end do
c
c step1
c -----
c note:  if a different tolerance is desired throughout this
c algorithm, change tolerance here:
c
      tol=sqrt(epsm)

      diagmx=a(1,1)
      diagmn=a(1,1)
      if(n.eq. 1 ) go to 40
      do i=2,n
        if(a(i,i).lt.diagmn) diagmn=a(i,i)
        if(a(i,i).gt.diagmx) diagmx=a(i,i)
      end do

   40 posmax=max(diagmx, 0.0D+00 )
c
c diagmn .le. 0
c
      if(diagmn.gt.posmax*tol) go to 110
c     if(diagmn.le.posmax*tol)
c     then
        amu=tol*(posmax-diagmn)-diagmn
        if(amu.ne. 0.0D+00 ) go to 90
c       if(amu.eq. 0.0D+00 )
c       then
c
c find largest off-diagonal element of a
c
          offmax=0.0D+00
          if(n.eq. 1 ) go to 70
          do 60 i=2,n
            im1=i-1
            do j=1,im1
              if(abs(a(i,j)).gt.offmax) offmax=abs(a(i,j))
            end do
   60     continue
   70     amu=offmax
          if(amu.ne. 0.0D+00 ) go to 80
c         if(amu.eq. 0.0D+00 )
c         then
            amu = 1.00D+00
            go to 90
c         else
   80       amu=amu* ( 1.00D+00 + tol )
c         end if
c       end if
c
c a=a + mu*i
c
   90   do 100 i=1,n
          a(i,i)=a(i,i)+amu
  100   continue
        diagmx=diagmx+amu
c     end if
c
c step2
c -----
c copy lower triangular part of "a" to upper triangular part
c and diagonal of "a" to udiag
c
  110 continue
      do 130 j=1,n
        udiag(j)=a(j,j)
        if(j.eq.n) go to 130
        jp1=j+1
        do 120 i=jp1,n
          a(j,i)=a(i,j)
  120   continue
  130 continue
c
      call choldc_u(nr,n,a,diagmx,tol,addmax)
c
c
c step3
c -----
c if addmax=0, "a" was positive definite going into step 2,
c the ll+ decomposition has been done, and we return.
c otherwise, addmax>0.  perturb "a" so that it is safely
c diagonally dominant and find ll+ decomposition
c
      if(addmax.le. 0.0D+00 ) go to 220
c     if(addmax.gt. 0.0D+00 )
c     then
c
c restore original "a" (lower triangular part and diagonal)
c
        do 150 j=1,n
          a(j,j)=udiag(j)
          if(j.eq.n) go to 150
          jp1=j+1
          do 140 i=jp1,n
            a(i,j)=a(j,i)
  140     continue
  150   continue
c
c find sdd such that a+sdd*i is safely positive definite
c note:  evmin<0 since a is not positive definite;
c
        evmin = 0.0D+00
        evmax=a(1,1)
        do 200 i=1,n
          offrow = 0.0D+00
          if(i.eq. 1 ) go to 170
          im1=i-1
          do 160 j=1,im1
            offrow=offrow+abs(a(i,j))
  160     continue
  170     if(i.eq.n) go to 190
          ip1=i+1
          do 180 j=ip1,n
            offrow=offrow+abs(a(j,i))
  180     continue
  190     evmin=min(evmin,a(i,i)-offrow)
          evmax=max(evmax,a(i,i)+offrow)
  200   continue
        sdd=tol*(evmax-evmin)-evmin
c
c perturb "a" and decompose again
c
        amu=min(sdd,addmax)
        do i=1,n
          a(i,i)=a(i,i)+amu
          udiag(i)=a(i,i)
        end do
c
c "a" now guaranteed safely positive definite
c
        call choldc_u(nr,n,a, 0.0D+00, tol, addmax )
c     end if
c
c unscale hessian and cholesky decomposition matrix
c
  220 do 260 j=1,n
        do i=j,n
          a(i,j)=sx(i)*a(i,j)
        end do
        if(j.eq. 1 ) go to 250
        jm1=j-1
        do 240 i=1,jm1
          a(i,j)=sx(i)*sx(j)*a(i,j)
  240   continue
  250   udiag(j)=udiag(j)*sx(j)*sx(j)
  260 continue
      return
      end
