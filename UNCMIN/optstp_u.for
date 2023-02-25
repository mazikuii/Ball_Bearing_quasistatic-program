      subroutine optstp_u(n,xpls,fpls,gpls,x,itncnt,icscmx,
     & itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,iprint,
     & lounit,msg)

c*********************************************************************72
c
cc OPTSTP checks the optimization stopping criteria.
c
c  Discussion:
c
c    find whether the algorithm should terminate, due to any
c    of the following:
c    1) problem solved within user tolerance
c    2) convergence within user tolerance
c    3) iteration limit reached
c    4) divergence or too restrictive maximum step (stepmx) suspected
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
c n            --> dimension of problem
c xpls(n)      --> new iterate x(k)
c fpls         --> function value at new iterate f(xpls)
c gpls(n)      --> gradient at new iterate, g(xpls), or approximate
c x(n)         --> old iterate x(k-1)
c itncnt       --> current iteration k
c icscmx      <--> number consecutive steps .ge. stepmx
c                  (retain value between successive calls)
c itrmcd      <--  termination code
c gradtl       --> tolerance at which relative gradient considered close
c                  enough to zero to terminate algorithm
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c sx(n)        --> diagonal scaling matrix for x
c fscale       --> estimate of scale of objective function
c itnlim       --> maximum number of allowable iterations
c iretcd       --> return code
c mxtake       --> boolean flag indicating step of maximum length used
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c msg          --> if msg includes a term 8, suppress output
c
      implicit double precision ( a-h, o-z )

      integer n,itncnt,icscmx,itrmcd,itnlim
      dimension sx(n)
      dimension xpls(n),gpls(n),x(n)
      logical mxtake

      itrmcd=0
c
c last global step failed to locate a point lower than x
      if(iretcd.ne. 1) go to 10
c     if(iretcd.eq. 1)
c     then
        jtrmcd=3
        go to 50
c     end if
   10 continue
c
c find direction in which relative gradient maximum.
c check whether within tolerance
c
      d=max(abs(fpls),fscale)
      rgx = 0.0D+00
      do i=1,n
        relgrd=abs(gpls(i))*max(abs(xpls(i)), 1.0D+00 /sx(i))/d
        rgx=max(rgx,relgrd)
      end do

      jtrmcd=1
      if(rgx.le.gradtl) go to 50

      if(itncnt.eq.0) return
c
c find direction in which relative stepsize maximum
c check whether within tolerance.
c
      rsx=0.0D+00
      do i=1,n
        relstp=abs(xpls(i)-x(i))/max(abs(xpls(i)), 1.0D+00 /sx(i))
        rsx=max(rsx,relstp)
      end do
      jtrmcd=2
      if(rsx.le.steptl) go to 50
c
c check iteration limit
c
      jtrmcd=4
      if(itncnt.ge.itnlim) go to 50
c
c check number of consecutive steps.ne.stepmx
c
      if(mxtake) go to 40
c     if(.not.mxtake)
c     then
        icscmx=0
        return
c     else
   40   continue
        if (mod(msg/8,2).eq. 0 .and. iprint.gt. 1)write(lounit,1000)
        icscmx=icscmx+1
        if(icscmx.lt. 5) return
        jtrmcd=5
c     end if
c
c
c print termination code
c
   50 itrmcd=jtrmcd
      if (mod(msg/8,2) .eq. 0) go to(60,70,80,90,100), itrmcd
      go to 110
   60 if(iprint.gt. 1)write(lounit,1010)
      go to 110
   70 if(iprint.gt. 1)write(lounit,1020)
      go to 110
   80 if(iprint.gt. 1)write(lounit,1030)
      go to 110
   90 if(iprint.gt. 1)write(lounit,1040)
      go to 110
  100 if(iprint.gt. 1)write(lounit,1050)

  110 return

 1000 format(48h0optstp    step of maximum length (stepmx) taken)
 1010 format(43h0optstp    relative gradient close to zero./
     &       48h optstp    current iterate is probably solution.)
 1020 format(48h0optstp    successive iterates within tolerance./
     &       48h optstp    current iterate is probably solution.)
 1030 format(52h0optstp    last global step failed to locate a point,
     &       14h lower than x./
     &       51h optstp    either x is an approximate local minimum,
     &       17h of the function,/
     &       50h optstp    the function is too non-linear for this,
     &       11h algorithm,/
     &       34h optstp    or steptl is too large.)
 1040 format(36h0optstp    iteration limit exceeded./
     &       28h optstp    algorithm failed.)
 1050 format(39h0optstp    maximum step size exceeded 5,
     &       19h consecutive times./
     &       50h optstp    either the function is unbounded below,/
     &       47h optstp    becomes asymptotic to a finite value,
     &       30h from above in some direction,/
     &       33h optstp    or stepmx is too small)
      end
