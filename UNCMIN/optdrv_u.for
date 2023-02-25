      subroutine optdrv_u(nr,n,x,fcn,d1fcn_u,d2fcn_u,typsiz,fscale,
     &     method,iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,
     &     dlt,gradtl,stepmx,steptl,xpls,fpls,gpls,itrmcd,
     &     a,udiag,g,p,sx,wrk0,wrk1,wrk2,wrk3)

c*********************************************************************72
c
cc OPTDRV is a driver for the nonlinear optimization code.
c
c  Discussion:
c
c    driver for non-linear optimization problem
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
c x(n)         --> on entry: estimate to a root of fcn
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                            fcn: r(n) --> r(1)
c d1fcn        --> (optional) name of subroutine to evaluate gradient
c                  of fcn.  must be declared external in calling routine
c d2fcn        --> (optional) name of subroutine to evaluate hessian of
c                  of fcn.  must be declared external in calling routine
c typsiz(n)    --> typical size for each component of x
c fscale       --> estimate of scale of objective function
c method       --> algorithm to use to solve minimization problem
c                    =1 line search
c                    =2 double dogleg
c                    =3 more-hebdon
c iexp         --> =1 if optimization function fcn is expensive to
c                  evaluate, =0 otherwise.  if set then hessian will
c                  be evaluated by secant update instead of
c                  analytically or by finite differences
c msg         <--> on input:  (.gt.0) message to inhibit certain
c                    automatic checks
c                  on output: (.lt.0) error code; =0 no error
c ndigit       --> number of good digits in optimization function fcn
c itnlim       --> maximum number of allowable iterations
c iagflg       --> =1 if analytic gradient supplied
c iahflg       --> =1 if analytic hessian supplied
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c dlt          --> trust region radius
c gradtl       --> tolerance at which gradient considered close
c                  enough to zero to terminate algorithm
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c xpls(n)     <--> on exit:  xpls is local minimum
c fpls        <--> on exit:  function value at solution, xpls
c gpls(n)     <--> on exit:  gradient at solution xpls
c itrmcd      <--  termination code
c a(n,n)       --> workspace for hessian (or estimate)
c                  and its cholesky decomposition
c udiag(n)     --> workspace (for diagonal of hessian)
c g(n)         --> workspace (for gradient at current iterate)
c p(n)         --> workspace for step
c sx(n)        --> workspace (for diagonal scaling matrix)
c wrk0(n)      --> workspace
c wrk1(n)      --> workspace
c wrk2(n)      --> workspace
c wrk3(n)      --> workspace
c
c
c internal variables
c ------------------
c analtl           tolerance for comparison of estimated and
c                  analytical gradients and hessians
c epsm             machine epsilon
c f                function value: fcn(x)
c itncnt           current iteration, k
c rnf              relative noise in optimization function fcn.
c                       noise=10.**(-ndigit)
c
      implicit double precision ( a-h, o-z )

      double precision r8_epsilon
      dimension x(n),xpls(n),g(n),gpls(n),p(n)
      dimension typsiz(n),sx(n)
      dimension a(nr,1),udiag(n)
      dimension wrk0(n),wrk1(n),wrk2(n),wrk3(n)
      logical mxtake,noupdt
      external fcn,d1fcn_u,d2fcn_u
c
c initialization
c
      do i=1,n
        p(i) = 0.0D+00
      end do

      itncnt=0
      iretcd=-1

      epsm = r8_epsilon ( )

      call optchk_u(n,x,typsiz,sx,fscale,gradtl,itnlim,ndigit,epsm,
     &     dlt,method,iexp,iagflg,iahflg,stepmx,msg,iprint,lounit)
      if(msg.lt.0) return
      rnf=max( 10.0D+00 **(-ndigit),epsm)
      analtl=max( 1.0D-02, sqrt(rnf))

      if(mod(msg/8,2).eq. 1) go to 20
      if(iprint.gt. 1)write(lounit,1010)
      if(iprint.gt. 1)write(lounit,1000) (typsiz(i),i=1,n)
      if(iprint.gt. 1)write(lounit,1020)
      if(iprint.gt. 1)write(lounit,1000) (sx(i),i=1,n)
      if(iprint.gt. 1)write(lounit,1030) fscale
      if(iprint.gt. 1)write(lounit,1040)
     & ndigit,iagflg,iahflg,iexp,method,itnlim,epsm
      if(iprint.gt. 1)write(lounit,1050) stepmx,steptl,gradtl,dlt,rnf,
     & analtl
   20 continue
c
c evaluate fcn(x)
c
      call fcn(n,x,f)
c
c evaluate analytic or finite difference gradient and check analytic
c gradient, if requested.
c
      if (iagflg .eq. 1) go to 30
c     if (iagflg .eq. 0)
c     then
        call fstofd_u (1, 1, n, x, fcn, f, g, sx, rnf, wrk, 1)
        go to 40
c
   30 call d1fcn_u (n, x, g)
      if (mod(msg/2,2) .eq. 1) go to 40
c     if (mod(msg/2,2).eq.0)
c     then
        call grdchk_u(n,x,fcn,f,g,typsiz,sx,fscale,
     &  rnf,analtl,wrk1,msg,iprint,lounit)
        if (msg .lt. 0) return
   40 continue

      call optstp_u(n,x,f,g,wrk1,itncnt,icscmx,
     &            itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,
     &            iprint,lounit,msg)
      if(itrmcd.ne.0) go to 250

      if(iexp.ne. 1) go to 50
c
c if optimization function expensive to evaluate (iexp=1), then
c hessian will be obtained by secant updates.  get initial hessian.
c
      call hsnint_u(nr,n,a,sx,method)
      go to 90
   50 continue
c
c evaluate analytic or finite difference hessian and check analytic
c hessian if requested (only if user-supplied analytic hessian
c routine d2fcn fills only lower triangular part and diagonal of a).
c
      if (iahflg .eq. 1) go to 60
c     if (iahflg .eq. 0)
c     then
       if (iagflg .eq. 1) call fstofd_u (nr, n, n, x, d1fcn_u, g, a, sx,
     &      rnf, wrk1, 3)
         if (iagflg .ne. 1) call sndofd_u (nr, n, x, fcn, f, a, sx, rnf,
     &      wrk1, wrk2)
         go to 80
c
c     else
   60    if (mod(msg/4,2).eq.0) go to 70
c        if (mod(msg/4, 2) .eq. 1)
c        then
            call d2fcn_u (nr, n, x, a)
            go to 80
c
c        else
   70  call heschk_u (nr, n, x, fcn, d1fcn_u, d2fcn_u, f, g, a, typsiz,
     &         sx,rnf,analtl,iagflg,udiag,wrk1,wrk2,msg,iprint,lounit)
c
c           heschk evaluates d2fcn and checks it against the finite
c           difference hessian which it calculates by calling fstofd
c           (if iagflg .eq. 1) or sndofd (otherwise).
c
            if (msg .lt. 0) return
   80 continue

   90 if(mod(msg/8,2).eq.0 .and. iprint.gt.0)then
        call result_u(nr,n,x,f,g,a,p,itncnt,1)
       end if
c
c
c iteration
c 
  100 itncnt=itncnt+1
c
c find perturbed local model hessian and its ll+ decomposition
c (skip this step if line search or dogstep techniques being used with
c secant updates.  cholesky decomposition l already obtained from
c secfac.)
c
      if(iexp.eq. 1 .and. method.ne. 3) go to 120
  110   call chlhsn_u(nr,n,a,epsm,sx,udiag)
  120 continue
c
c solve for newton step:  ap=-g
c
      do i=1,n
        wrk1(i)=-g(i)
      end do

      call lltslv_u(nr,n,a,p,wrk1)
c
c decide whether to accept newton step  xpls=x + p
c or to choose xpls by a global strategy.
c
      if (iagflg .ne. 0 .or. method .eq. 1) go to 140
      dltsav = dlt
      if (method .eq. 2) go to 140
      amusav = amu
      dlpsav = dltp
      phisav = phi
      phpsav = phip0
  140 if(method.eq. 1)
     &     call lnsrch_u(n,x,f,g,p,xpls,fpls,fcn,mxtake,iretcd,
     &     stepmx,steptl,sx,iprint,lounit)
      if(method.eq. 2)
     &     call dogdrv_u(nr,n,x,f,g,a,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,wrk0,wrk1,wrk2,wrk3,iprint,lounit)
      if(method.eq. 3)
     &     call hookdr_u(nr,n,x,f,g,a,udiag,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,amu,dltp,phi,phip0,wrk0,
     &     wrk1,wrk2,epsm,itncnt,iprint,lounit)
c
c if could not find satisfactory step and forward difference
c gradient was used, retry using central difference gradient.
c
      if (iretcd .ne. 1 .or. iagflg .ne. 0) go to 150
c     if (iretcd .eq. 1 .and. iagflg .eq. 0)
c     then
c
c        set iagflg for central differences
c
         iagflg = -1
         if(iprint.gt. 1)write(lounit,1060) itncnt

         call fstocd_u (n, x, fcn, sx, rnf, g)
         if (method .eq. 1) go to 120
         dlt = dltsav
         if (method .eq. 2) go to 120
         amu = amusav
         dltp = dlpsav
         phi = phisav
         phip0 = phpsav
         go to 110
c     end if
c
c calculate step for output
c
  150 continue

      do i = 1, n
         p(i) = xpls(i) - x(i)
      end do
c
c calculate gradient at xpls
c
      if (iagflg .eq. (-1)) go to 170
      if (iagflg .eq. 0) go to 180
c
c analytic gradient
c
      call d1fcn_u (n, xpls, gpls)
      go to 190
c
c central difference gradient
c
  170 call fstocd_u (n, xpls, fcn, sx, rnf, gpls)
      go to 190
c
c forward difference gradient
c
  180 call fstofd_u (1, 1, n, xpls, fcn, fpls, gpls, sx, rnf, wrk, 1)
  190 continue
c
c check whether stopping criteria satisfied
c
      call optstp_u(n,xpls,fpls,gpls,x,itncnt,icscmx,
     &            itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,
     &            iprint,lounit,msg)
      if(itrmcd.ne.0) go to 240
c
c evaluate hessian at xpls
c
      if(iexp.eq.0) go to 200

      if(method.eq. 3)
     &     call secunf_u(nr,n,x,g,a,udiag,xpls,gpls,epsm,itncnt,rnf,
     &     iagflg,noupdt,wrk1,wrk2,wrk3)
      if(method.ne. 3)
     &     call secfac_u(nr,n,x,g,a,xpls,gpls,epsm,itncnt,rnf,iagflg,
     &     noupdt,wrk0,wrk1,wrk2,wrk3)
      go to 220
  200 if(iahflg.eq. 1) go to 210
      if(iagflg.eq. 1)
     &     call fstofd_u(nr,n,n,xpls,d1fcn_u,gpls,a,sx,rnf,wrk1,3)

      if(iagflg.ne. 1) then
        call sndofd_u(nr,n,xpls,fcn,fpls,a,sx,rnf,wrk1,wrk2)
      end if

      go to 220
  210 call d2fcn_u(nr,n,xpls,a)
  220 continue
      if(mod(msg/16,2).eq. 1 .and. iprint.gt.0)then
        call result_u(nr,n,xpls,fpls,gpls,a,p,itncnt,1)
      end if
c
c x <-- xpls  and  g <-- gpls  and  f <-- fpls
c
      f=fpls
      do i=1,n
        x(i)=xpls(i)
        g(i)=gpls(i)
      end do
      go to 100
c
c termination
c reset xpls,fpls,gpls,  if previous iterate solution
c
  240 if(itrmcd.ne. 3) go to 270
  250 continue
      fpls=f
      do i=1,n
        xpls(i)=x(i)
        gpls(i)=g(i)
      end do
c
c print results
c
  270 continue
      if(mod(msg/8,2).eq.0 .and. iprint.gt.0)then
        call result_u(nr,n,xpls,fpls,gpls,a,p,itncnt,0)
      end if
      msg=0
      return

 1000 format(14h optdrv       ,5(g16.8,3x))
 1010 format(20h0optdrv    typical x)
 1020 format(40h optdrv    diagonal scaling matrix for x)
 1030 format(22h optdrv    typical f =,g16.8)
 1040 format(40h0optdrv    number of good digits in fcn=,i5/
     &       27h optdrv    gradient flag  =,i5,18h   (=1 if gradient,
     &       ' supplied)'/
     &       27h optdrv    hessian flag   =,i5,'   (=1 if hessian',
     &       ' supplied)'/
     &       27h optdrv    expense flag   =,i5, 9h   (=1 if,
     &       ' function expensive)'/
     &       ' optdrv    method to use=',i5/
     &       '           1=line search'/
     &       '           2=double dogleg'/
     &       '           3=more-hebdon'/
     &       27h optdrv    iteration limit=,i5/
     &       27h optdrv    machine epsilon=,g16.8)
 1050 format(33h0optdrv    maximum step size    =,g16.8/
     &       33h optdrv    step tolerance       =,g16.8/
     &       33h optdrv    gradient tolerance   =,g16.8/
     &       33h optdrv    trust region radius  =,g16.8/
     &       33h optdrv    rel noise in fcn     =,g16.8/
     &       33h optdrv    anal-fd tolerance    =,g16.8)
 1060 format(52h optdrv    shift from forward to central differences,
     &   14h in iteration , i5)
      end
