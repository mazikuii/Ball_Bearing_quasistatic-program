      subroutine secunf_u(nr,n,x,g,a,udiag,xpls,gpls,epsm,itncnt,
     &     rnf,iagflg,noupdt,s,y,t)

c*********************************************************************72
c
cc SECUNF updates the Hessian by the BFGS unfactored method.
c
c  Discussion:
c
c    update hessian by the bfgs unfactored method
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
c x(n)         --> old iterate, x(k-1)
c g(n)         --> gradient or approximate at old iterate
c a(n,n)      <--> on entry: approximate hessian at old iterate
c                    in upper triangular part (and udiag)
c                  on exit:  updated approx hessian at new iterate
c                    in lower triangular part and diagonal
c                  (lower triangular part of symmetric matrix)
c udiag        --> on entry: diagonal of hessian
c xpls(n)      --> new iterate, x(k)
c gpls(n)      --> gradient or approximate at new iterate
c epsm         --> machine epsilon
c itncnt       --> iteration count
c rnf          --> relative noise in optimization function fcn
c iagflg       --> =1 if analytic gradient supplied, =0 otherwise
c noupdt      <--> boolean: no update yet
c                  (retain value between successive calls)
c s(n)         --> workspace
c y(n)         --> workspace
c t(n)         --> workspace
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      double precision r8vec_norm_l2
      dimension x(n),g(n),xpls(n),gpls(n)
      dimension a(nr,1)
      dimension udiag(n)
      dimension s(n),y(n),t(n)
      logical noupdt,skpupd
c
c  copy hessian in upper triangular part and udiag to
c  lower triangular part and diagonal
c
      do 20 j=1,n
        a(j,j)=udiag(j)
        if(j.eq.n) go to 20
        jp1=j+1
        do i=jp1,n
          a(i,j)=a(j,i)
        end do
   20 continue

      if(itncnt.eq. 1) noupdt=.true.

      do i=1,n
        s(i)=xpls(i)-x(i)
        y(i)=gpls(i)-g(i)
      end do

      den1 = r8vec_dot ( n, s, y )
      snorm2 = r8vec_norm_l2 ( n, s )
      ynrm2 = r8vec_norm_l2 ( n, y )
      if(den1 .lt.sqrt(epsm)*snorm2*ynrm2) go to 110
c     if(den1 .ge.sqrt(epsm)*snorm2*ynrm2)
c     then
        call mvmlts_u(nr,n,a,s,t)
        den2 = r8vec_dot ( n, s, t )
        if(.not. noupdt) go to 60
c       if(noupdt)
c       then
c
c         h <-- ((s+)y/(s+)hs)h
c
          gam=den1/den2
          den2=gam*den2
          do j=1,n
            t(j)=gam*t(j)
            do i=j,n
              a(i,j)=gam*a(i,j)
            end do
          end do
          noupdt=.false.
c       end if
   60   skpupd=.true.
c
c  check update condition on row i
c
        do 70 i=1,n
          tol=rnf*max(abs(g(i)),abs(gpls(i)))
          if(iagflg.eq.0) tol=tol/sqrt(rnf)
          if(abs(y(i)-t(i)).lt.tol) go to 70
c         if(abs(y(i)-t(i)).ge.tol)
c         then
            skpupd=.false.
            go to 80
c         end if
   70   continue
   80   if(skpupd) go to 110
c
c  bfgs update
c
          do j=1,n
            do i=j,n
              a(i,j)=a(i,j)+y(i)*y(j)/den1-t(i)*t(j)/den2
            end do
          end do
c       end if
c     end if
  110 return
      end
