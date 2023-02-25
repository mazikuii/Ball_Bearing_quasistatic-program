      subroutine neddlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol,
     *                  mxtake,dlt,qtf,scalex,scaley,fvec,d,xprev,
     *                  ssd,v,wa,fprev,xp, fp, fpnorm,retcd,gcnt,
     *                  priter,iter)
	use neintface
	implicit none
      integer ldr, n, retcd, gcnt, priter, iter
      double precision  fcnorm, stepmx, xtol, dlt
      double precision  rjac(n,n), dn(n), g(n), xc(n), qtf(n)
      double precision  scalex(n),scaley(n), d(n)
      double precision  xprev(n)
      double precision  ssd(n), v(n), wa(n), fprev(n)
      logical mxtake
      external fvec
      double precision  xp(n), fp(n), fpnorm

c-------------------------------------------------------------------------
c
c     Find a next iterate xp by the double dogleg method
c
c     Arguments
c
c     In       n       Integer         size of problem: dimension x, f
c     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of Rjac
c     Inout    dn      Real(*)         newton direction
c     Inout    g       Real(*)         gradient at current point
c                                      trans(jac)*f()
c     In       xc      Real(*)         current iterate
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       stepmx  Real            maximum stepsize
c     In       xtol    Real            x-tolerance (stepsize)
c     Out      mxtake  Logical         true if maximum step taken
c     Inout    dlt     Real            on input: initial trust region radius
c                                                if -1 then set to something
c                                                reasonable
c                                      on output: final value
c                                      ! Do not modify between calls while
c                                        still iterating
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       scalex  Real(*)         scaling factors for x()
c     In       fvec    Name            name of subroutine to evaluate f(x)
c                                      ! must be declared external in caller
c     Wk       d       Real(*)         work vector
c     Wk       xprev   Real(*)         work vector
c     Wk       ssd     Real(*)         work vector
c     Wk       v       Real(*)         work vector
c     Wk       wa      Real(*)         work vector
c     Wk       fprev   Real(*)         work vector
c     Out      xp      Real(*)         new x()
c     Out      fp      Real(*)         new f(xp)
c     Out      fpnorm  Real            new .5*||f(xp)||**2
c     Out      retcd   Integer         return code
c                                       0  new satisfactory x() found
c                                       1  no  satisfactory x() found
c     Out      gcnt    Integer         number of steps taken
c     In       priter  Integer         print flag
c                                       -1 no intermediate printing
c                                       >0 yes for print of intermediate results
c     In       iter    Integer         current iteration (only used for above)
c
c     All vectors at least size n
c
c-------------------------------------------------------------------------

      integer i
      double precision  dnlen,ssdlen,alpha,beta,lambda,vlen,vssdag,fpred
      double precision  eta,gamma,fpnsav,oarg(7)
      double precision  dnrm2, ddot
      logical nwtake
      integer dtype

      integer idamax

      double precision Rone, Rtwo, Rten, Rhalf, Rp2, Rp8  
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
      parameter(Rp2 = Rtwo/Rten, Rp8 = Rone - Rp2)

	xprev=xp
	fprev=fp
	fcnorm=fpnorm

c     length scaled newton direction

      do i=1,n
         wa(i) = dn(i) * scalex(i)
      enddo

      dnlen = dnrm2(n, wa, 1)

c     steepest descent direction and length

      do i=1,n
         wa(i) = g(i) / scalex(i)
      enddo
      alpha = dnrm2(n,wa,1)**2

      do i=1,n
         d(i) = wa(i) / scalex(i)
      enddo
      call dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      call dcopy(n, wa, 1, ssd, 1)
      call dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*sqrt(alpha)/beta

c     set trust radius to ssdlen or dnlen if required

      if( dlt .eq. -Rone ) then
         dlt = min(ssdlen, stepmx)
      elseif( dlt .eq. -Rtwo ) then
         dlt = min(dnlen, stepmx)
      endif

c     calculate double dogleg parameter

      gamma = alpha*alpha/(-beta*ddot(n,g,1,dn,1))
      eta = Rp2 + Rp8*gamma

      do 10 i=1,n
          v(i) = eta*scalex(i)*dn(i) - ssd(i)
  10  continue

      vssdag = ddot(n,v,1,ssd,1)
      vlen   = dnrm2(n,v,1)**2

      retcd = 4
      gcnt  = 0

  20  continue

c       find new step by double dogleg algorithm

        call Nedogstp(n,scalex,dn,dnlen,dlt,nwtake,vssdag,vlen,
     *              ssd,v,ssdlen,eta,d,dtype,lambda)

c       compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

        call dcopy(n,d,1,wa,1)
        call dtrmv('U','N','N',n,rjac,ldr,wa,1)
        call daxpy(n, Rone, qtf,1,wa,1)
        fpred = Rhalf * dnrm2(n,wa,1)**2

c       evaluate function at xp = xc + d

        do 30 i = 1,n
           xp(i) = xc(i) + d(i)
  30    continue

	fpnorm=nefmin(xp,fvec,fp,scaley)
!        call nwfvec(xp,n,fvec,fp,fpnorm)
        gcnt = gcnt + 1

c       check whether the global step is acceptable

        oarg(2) = dlt
        call netrup(n,fcnorm,g,d,scalex,nwtake,stepmx,xtol,dlt,mxtake,
     *              fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm,wa)


      if(retcd .gt. 1) goto 20

      return
      end