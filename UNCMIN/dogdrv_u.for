
      subroutine dogdrv_u(nr,n,x,f,g,a,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,sc,wrk1,wrk2,wrk3,iprint,lounit)

c*********************************************************************72
c
cc DOGDRV finds the next Newton iterate by the double dogleg method.
c
c  Discussion:
c
c    find a next newton iterate (xpls) by the double dogleg method
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
c g(n)         --> gradient  at old iterate, g(x), or approximate
c a(n,n)       --> cholesky decomposition of hessian
c                  in lower triangular part and diagonal
c p(n)         --> newton step
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c fcn          --> name of subroutine to evaluate function
c sx(n)        --> diagonal scaling matrix for x
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c dlt         <--> trust region radius
c                  (retain value between successive calls)
c iretcd      <--  return code
c                    =0 satisfactory xpls found
c                    =1 failed to find satisfactory xpls sufficiently
c                       distinct from x
c mxtake      <--  boolean flag indicating step of maximum length used
c sc(n)        --> workspace (current step)
c wrk1(n)      --> workspace (and place holding argument to tregup)
c wrk2(n)      --> workspace
c wrk3(n)      --> workspace
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),xpls(n),g(n),p(n)
      dimension sx(n)
      dimension sc(n),wrk1(n),wrk2(n),wrk3(n)
      dimension a(nr,1)
      logical fstdog,nwtake,mxtake
      external fcn

      iretcd=4
      fstdog=.true.
      tmp = 0.0D+00
      do i=1,n
        tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
      end do
      rnwtln=sqrt(tmp)
      if(iprint.gt. 1 )write(lounit,1030) rnwtln

   20 continue
c
c find new step by double dogleg algorithm
c
      call dogstp_u(nr,n,g,a,p,sx,rnwtln,dlt,nwtake,fstdog,
     &     wrk1,wrk2,cln,eta,sc,iprint,lounit,stepmx)
c
c check new point and update trust region
c
      call tregup_u(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,dlt,
     &     iretcd,wrk3,fplsp,xpls,fpls,mxtake,iprint,lounit,2,wrk1)

      if(iretcd.le. 1 ) return
      go to 20
 1000 format(42h dogdrv    initial trust region not given.,
     &       22h  compute cauchy step.)
 1010 format(18h dogdrv    alpha =,g16.8/
     &       18h dogdrv    beta  =,g16.8/
     &       18h dogdrv    dlt   =,g16.8/
     &       18h dogdrv    nwtake=,l1    )
 1020 format(28h dogdrv    current step (sc))
 1030 format(18h0dogdrv    rnwtln=,g16.8)
 1040 format(14h dogdrv       ,5(g16.8,3x))
      end
