      subroutine Nedogstp(n,scalex,dn,dnlen,dlt,nwtake,vssdag,vlen,
     *                  ssd,v,ssdlen,eta,d,dtype,lambda)
      integer n
      double precision  scalex(*), dn(*), ssd(*), v(*), d(*)
      double precision  dnlen, dlt, vssdag, vlen, ssdlen, eta, lambda
      logical nwtake
      integer dtype

c-------------------------------------------------------------------------
c
c     Find a new step by the double dogleg algorithm
c     Internal routine for nwddlg
c
c     Arguments
c
c     In       n       Integer         size of problem
c     In       scalex  Real(*)         scalig factors x()
c     In       dn      Real(*)         current newton step
c     Out      dnlen   Real            length dn()
c     In       dlt     Real            current trust region radius
c     Inout    nwtake  Logical         true if newton step taken
c     In       vssdag  Real            (internal)
c     In       vlen    Real            (internal)
c     In       ssd     Real(*)         (internal) steepest descent direction
c     In       v       Real(*)         (internal) eta * dn() - ssd()
c     In       ssdlen  Real            (internal) length ssd
c     In       eta     Real            (internal) double dogleg parameter
c     Out      d       Real(*)         new step for x()
c     Out      dtype   Integer         steptype
c                                       1 steepest descent
c                                       2 full newton direction
c                                       3 partial newton step
c                                       4 combination of dn and ssd
c     Out      lambda  Real            weight of eta*dn() in d()
c                                      closer to 1 ==> more of eta*dn()
c
c-----------------------------------------------------------------------

      nwtake = .false.

      if(dnlen .le. dlt) then

c         Newton step smaller than trust radius ==> take it

          nwtake = .true.
          call dcopy(n, dn, 1, d, 1)
          dlt = dnlen
          dtype = 2

      elseif(eta*dnlen .le. dlt) then

c         take partial step in newton direction

          call dcopy(n, dn, 1, d, 1)
          call dscal(n, dlt / dnlen, d, 1)
          dtype = 3

      elseif(ssdlen .ge. dlt) then

c         take step in steepest descent direction

          call dcopy(n, ssd, 1, d, 1)
          call dscal(n, dlt / ssdlen, d, 1)
          call vunsc(n,d,scalex)
          dtype = 1

      else

c         calculate convex combination of ssd and eta*p
c         which has scaled length dlt

          lambda =(-vssdag+sqrt(vssdag**2-vlen*(ssdlen**2-dlt**2)))/vlen
          call dcopy(n, ssd, 1, d, 1)
          call daxpy(n, lambda, v, 1, d, 1)
          call vunsc(n,d,scalex)
          dtype = 4

      endif

      return
      end
