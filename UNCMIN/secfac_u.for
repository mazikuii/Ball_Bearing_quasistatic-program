      subroutine secfac_u(nr,n,x,g,a,xpls,gpls,epsm,itncnt,rnf,
     &     iagflg,noupdt,s,y,u,w)

c*********************************************************************72
c
cc SECFAC updates the Hessian matrix by the BFGS factored method.
c
c  Discussion:
c
c    update hessian by the bfgs factored method
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
c a(n,n)      <--> on entry: cholesky decomposition of hessian in
c                    lower part and diagonal.
c                  on exit:  updated cholesky decomposition of hessian
c                    in lower triangular part and diagonal
c xpls(n)      --> new iterate, x(k)
c gpls(n)      --> gradient or approximate at new iterate
c epsm         --> machine epsilon
c itncnt       --> iteration count
c rnf          --> relative noise in optimization function fcn
c iagflg       --> =1 if analytic gradient supplied, =0 itherwise
c noupdt      <--> boolean: no update yet
c                  (retain value between successive calls)
c s(n)         --> workspace
c y(n)         --> workspace
c u(n)         --> workspace
c w(n)         --> workspace
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      double precision r8vec_norm_l2
      dimension x(n),xpls(n),g(n),gpls(n)
      dimension a(nr,1)
      dimension s(n),y(n),u(n),w(n)
      logical noupdt,skpupd

      if(itncnt.eq. 1) noupdt=.true.
      do i=1,n
        s(i)=xpls(i)-x(i)
        y(i)=gpls(i)-g(i)
      end do

      den1 = r8vec_dot ( n, x, y )

      snorm2 = r8vec_norm_l2 ( n, s )
      ynrm2 = r8vec_norm_l2 ( n, y )

      if ( den1 .lt. sqrt(epsm)*snorm2*ynrm2) go to 160
        call mvmltu_u(nr,n,a,s,u)
        den2 = r8vec_dot ( n, u, u )
c
c  l <-- sqrt(den1/den2)*l
c
        alp=sqrt(den1/den2)
        if(.not.noupdt) go to 40
          do j=1,n
            u(j)=alp*u(j)
            do i=j,n
              a(i,j)=alp*a(i,j)
            end do
          end do
          noupdt=.false.
          den2=den1
          alp = 1.0D+00
   40   skpupd=.true.
c
c  w = l(l+)s = hs
c
        call mvmltl_u(nr,n,a,u,w)
        i=1
        if(iagflg.ne.0) go to 50
          reltol=sqrt(rnf)
          go to 60
   50     reltol=rnf
   60   if(i.gt.n .or. .not.skpupd) go to 80
          if(abs(y(i)-w(i)) .lt. reltol*max(abs(g(i)),abs(gpls(i))))
     &         go to 70
            skpupd=.false.
            go to 60
   70       i=i+1
            go to 60
   80   if(skpupd) go to 160
c
c  w=y-alp*l(l+)s
c
          do i=1,n
            w(i)=y(i)-alp*w(i)
          end do
c
c  alp=1/sqrt(den1*den2)
c
          alp=alp/den1
c
c  u=(l+)/sqrt(den1*den2) = (l+)s/sqrt((y+)s * (s+)l(l+)s)
c
          do i=1,n
            u(i)=alp*u(i)
          end do
c
c  copy l into upper triangular part.  zero l.
c
          if(n.eq. 1) go to 130
          do i=2,n
            im1=i-1
            do j=1,im1
              a(j,i)=a(i,j)
              a(i,j)= 0.0D+00
            end do
          end do
c
c  find q, (l+) such that  q(l+) = (l+) + u(w+)
c
  130     call qrupdt_u(nr,n,a,u,w)
c
c  upper triangular part and diagonal of a now contain updated
c  cholesky decomposition of hessian.  copy back to lower
c  triangular part.
c
          if(n.eq. 1) go to 160
          do i=2,n
            im1=i-1
            do j=1,im1
              a(i,j)=a(j,i)
            end do
          end do
c       end if
c     end if
  160 return
      end
