      subroutine hookst_u(nr,n,g,a,udiag,p,sx,rnwtln,dlt,amu,
     &     dltp,phi,phip0,fstime,sc,nwtake,wrk0,epsm,iprint,lounit)

c*********************************************************************72
c
cc HOOKST finds the More-Hebdon stepsize.
c
c  Discussion:
c
c    find new step by more-hebdon algorithm
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
c g(n)         --> gradient at current iterate, g(x)
c a(n,n)       --> cholesky decomposition of hessian in
c                  lower triangular part and diagonal.
c                  hessian or approx in upper triangular part
c udiag(n)     --> diagonal of hessian in a(.,.)
c p(n)         --> newton step
c sx(n)        --> diagonal scaling matrix for n
c rnwtln       --> newton step length
c dlt         <--> trust region radius
c amu         <--> (retain value between successive calls)
c dltp         --> trust region radius at last exit from this routine
c phi         <--> (retain value between successive calls)
c phip0       <--> (retain value between successive calls)
c fstime      <--> boolean. =.true. if first entry to this routine
c                  during k-th iteration
c sc(n)       <--  current step
c nwtake      <--  boolean, =.true. if newton step taken
c wrk0         --> workspace
c epsm         --> machine epsilon
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_norm_l2
      dimension g(n),p(n),sx(n),sc(n),wrk0(n)
      dimension a(nr,1),udiag(n)
      logical nwtake,done
      logical fstime
c
c hi and alo are constants used in this routine.
c change here if other values are to be substituted.
c
      hi = 1.5D+00
      alo = 0.75D+00
      if(rnwtln.gt.hi*dlt) go to 20
c
c  take newton step
c
        nwtake=.true.
        do i=1,n
          sc(i)=p(i)
        end do
        dlt=min(dlt,rnwtln)
        amu = 0.0D+00
        if(iprint.gt. 1)write(lounit,1000)
        return
c
c  newton step not taken
c
   20   continue
        if(iprint.gt. 1)write(lounit,1010)
        nwtake=.false.
        if(amu .le. 0.0D+00 ) go to 30
          amu=amu- (phi+dltp) *((dltp-dlt)+phi)/(dlt*phip)
          if(iprint.gt. 1)write(lounit,1050) amu
   30   continue
        phi=rnwtln-dlt
        if(.not.fstime) go to 50
          do i=1,n
            wrk0(i)=sx(i)*sx(i)*p(i)
          end do
c
c  solve l*y = (sx**2)*p
c
          call forslv_u(nr,n,a,wrk0,wrk0)

          phip0 = - ( r8vec_norm_l2 ( n, wrk0 ) )**2 / rnwtln

          fstime=.false.
   50   phip=phip0
        amulo=-phi/phip
        amuup = 0.0D+00
        do i=1,n
          amuup=amuup+(g(i)*g(i))/(sx(i)*sx(i))
        end do
        amuup=sqrt(amuup)/dlt
        done=.false.
        if(iprint.gt. 1)write(lounit,1050) amu
        if(iprint.gt. 1)write(lounit,1080) phi
        if(iprint.gt. 1)write(lounit,1090) phip
        if(iprint.gt. 1)write(lounit,1060) amulo
        if(iprint.gt. 1)write(lounit,1070) amuup
c
c  test value of amu; generate next amu if necessary
c
   70   continue
        if(done) return
        if(iprint.gt. 1)write(lounit,1110)
        if(amu.ge.amulo .and. amu.le.amuup) go to 80
          amu=max(sqrt(amulo*amuup),amuup * 1.0D-03 )
          if(iprint.gt. 1)write(lounit,1050) amu
   80   continue
c
c  copy (h,udiag) to l
c  where h <-- h+amu*(sx**2) (do not actually change (h,udiag))
c
        do 100 j=1,n
          a(j,j)=udiag(j) + amu*sx(j)*sx(j)
          if(j.eq.n) go to 100
          jp1=j+1
          do i=jp1,n
            a(i,j)=a(j,i)
          end do
  100   continue
c
c       factor h=l(l+)
c
        call choldc_u(nr,n,a, 0.0D+00, sqrt(epsm),addmax)
c
c       solve h*p = l(l+)*sc = -g
c
        do i=1,n
          wrk0(i)=-g(i)
        end do

        call lltslv_u(nr,n,a,sc,wrk0)
        if(iprint.gt. 1)write(lounit,1040)
        if(iprint.gt. 1)write(lounit,1120) (sc(i),i=1,n)
c
c  reset h.  note since udiag has not been destroyed we need do
c  nothing here.  h is in the upper part and in udiag, still intact
c
        stepln = 0.0D+00
        do i=1,n
          stepln=stepln + sx(i)*sx(i)*sc(i)*sc(i)
        end do

        stepln=sqrt(stepln)
        phi=stepln-dlt
        do i=1,n
          wrk0(i)=sx(i)*sx(i)*sc(i)
        end do
        call forslv_u(nr,n,a,wrk0,wrk0)

        phip = - ( r8vec_norm_l2 ( n, wrk0 ) )**2 / stepln

        if(iprint.gt. 1)write(lounit,1100) dlt,stepln
        if(iprint.gt. 1)write(lounit,1080) phi
        if(iprint.gt. 1)write(lounit,1090) phip
        if((alo*dlt.gt.stepln .or. stepln.gt.hi*dlt) .and.
     &       (amuup-amulo.gt. 0.0D+00 )) go to 140
c
c  sc is acceptable hookstep
c
          if(iprint.gt. 1)write(lounit,1030)
          done=.true.
          go to 70
c       else
c
c  sc not acceptable hookstep.  select new amu
c
  140     continue
          if(iprint.gt. 1)write(lounit,1020)
          amulo=max(amulo,amu-(phi/phip))
          if(phi.lt. 0.0D+00 ) amuup=min(amuup,amu)
          amu=amu-(stepln*phi)/(dlt*phip)
          if(iprint.gt. 1)write(lounit,1050) amu
          if(iprint.gt. 1)write(lounit,1060) amulo
          if(iprint.gt. 1)write(lounit,1070) amuup
          go to 70
c       end if
c     end if
c
 1000 format(27h0hookst    take newton step)
 1010 format(32h0hookst    newton step not taken)
 1020 format(31h hookst    sc is not acceptable)
 1030 format(27h hookst    sc is acceptable)
 1040 format(28h hookst    current step (sc))
 1050 format(18h hookst    amu   =,g16.8)
 1060 format(18h hookst    amulo =,g16.8)
 1070 format(18h hookst    amuup =,g16.8)
 1080 format(18h hookst    phi   =,g16.8)
 1090 format(18h hookst    phip  =,g16.8)
 1100 format(18h hookst    dlt   =,g16.8/
     &       18h hookst    stepln=,g16.8)
 1110 format(23h0hookst    find new amu)
 1120 format(14h hookst       ,5(g16.8,3x))
      end
