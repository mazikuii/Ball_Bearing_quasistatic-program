      subroutine hookdr_u(nr,n,x,f,g,a,udiag,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,amu,dltp,phi,phip0,
     &     sc,xplsp,wrk0,epsm,itncnt,iprint,lounit)

c*********************************************************************72
c
cc HOOKDR finds the next Newton iterate by the More-Hebdon method.
c
c  Discussion:
c
c    find a next newton iterate (xpls) by the more-hebdon method
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
c x(n)         --> old iterate x(k-1)
c f            --> function value at old iterate, f(x)
c g(n)         --> gradient at old iterate, g(x), or approximate
c a(n,n)       --> cholesky decomposition of hessian in lower
c                  triangular part and diagonal.
c                  hessian in upper triangular part and udiag.
c udiag(n)     --> diagonal of hessian in a(.,.)
c p(n)         --> newton step
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c fcn          --> name of subroutine to evaluate function
c sx(n)        --> diagonal scaling matrix for x
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c dlt         <--> trust region radius
c iretcd      <--  return code
c                    =0 satisfactory xpls found
c                    =1 failed to find satisfactory xpls sufficiently
c                       distinct from x
c mxtake      <--  boolean flag indicating step of maximum length used
c amu         <--> (retain value between successive calls)
c dltp        <--> (retain value between successive calls)
c phi         <--> (retain value between successive calls)
c phip0       <--> (retain value between successive calls)
c sc(n)        --> workspace
c xplsp(n)     --> workspace
c wrk0(n)      --> workspace
c epsm         --> machine epsilon
c itncnt       --> iteration count
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n),p(n),xpls(n),sx(n)
      dimension a(nr,1),udiag(n)
      dimension sc(n),xplsp(n),wrk0(n)
      logical mxtake,nwtake
      logical fstime
      external fcn

      iretcd=4
      fstime=.true.
      tmp = 0.0D+00
      do i=1,n
        tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
      end do
      rnwtln=sqrt(tmp)
      if(iprint.gt. 1)write(lounit,1030) rnwtln
c
      if(itncnt.gt. 1) go to 50
c     if(itncnt.eq. 1)
c     then
        amu=0.0D+00
c
c       if first iteration and trust region not provided by user,
c       compute initial trust region.
c
        if(dlt.ne. ( -1.0D+00 )) go to 50
c       if(dlt.eq. ( -1.0D+00 ))
c       then
          alpha = 0.0D+00
          do i=1,n
            alpha=alpha+(g(i)*g(i))/(sx(i)*sx(i))
          end do

          beta=0.0D+00
          do i=1,n
            tmp = 0.0D+00
            do j=i,n
              tmp=tmp + (a(j,i)*g(j))/(sx(j)*sx(j))
            end do
            beta=beta+tmp*tmp
          end do
          dlt=alpha*sqrt(alpha)/beta
          dlt = min(dlt, stepmx)
          if(iprint.gt. 1)write(lounit,1000)
          if(iprint.gt. 1)write(lounit,1010) alpha,beta,dlt
c       end if
c     end if
c
   50 continue
c
c find new step by more-hebdon algorithm
c
      call hookst_u(nr,n,g,a,udiag,p,sx,rnwtln,dlt,amu,
     &     dltp,phi,phip0,fstime,sc,nwtake,wrk0,epsm,iprint,lounit)
      dltp=dlt
c
c check new point and update trust region
c
      call tregup_u(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,
     &   dlt,iretcd,xplsp,fplsp,xpls,fpls,mxtake,iprint,lounit,3,udiag)
      if(iretcd.le. 1) return
      go to 50
 1000 format(43h hookdr    initial trust region not given. ,
     &       21h compute cauchy step.)
 1010 format(18h hookdr    alpha =,g16.8/
     &       18h hookdr    beta  =,g16.8/
     &       18h hookdr    dlt   =,g16.8)
 1020 format(28h hookdr    current step (sc))
 1030 format(18h0hookdr    rnwtln=,g16.8)
 1040 format(14h hookdr       ,5(g16.8,3x))
      end
