      subroutine tregup_u(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,
     & dlt,iretcd,xplsp,fplsp,xpls,fpls,mxtake,iprint,lounit,method,
     & udiag)

c*********************************************************************72
c
cc TREGUP accepts the next iterate and updates the trust region.
c
c  Discussion:
c
c    decide whether to accept xpls=x+sc as the next iterate and update the
c    trust region dlt.
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
c a(n,n)       --> cholesky decomposition of hessian in
c                  lower triangular part and diagonal.
c                  hessian or approx in upper triangular part
c fcn          --> name of subroutine to evaluate function
c sc(n)        --> current step
c sx(n)        --> diagonal scaling matrix for x
c nwtake       --> boolean, =.true. if newton step taken
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c dlt         <--> trust region radius
c iretcd      <--> return code
c                    =0 xpls accepted as next iterate;
c                       dlt trust region for next iteration.
c                    =1 xpls unsatisfactory but accepted as next iterate
c                       because xpls-x .lt. smallest allowable
c                       step length.
c                    =2 f(xpls) too large.  continue current iteration
c                       with new reduced dlt.
c                    =3 f(xpls) sufficiently small, but quadratic model
c                       predicts f(xpls) sufficiently well to continue
c                       current iteration with new doubled dlt.
c xplsp(n)    <--> workspace (value needs to be retained between
c                  succesive calls of k-th global step)
c fplsp       <--> (retain value between successive calls)
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c mxtake      <--  boolean flag indicating step of maximum length used
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c method       --> algorithm to use to solve minimization problem
c                    =1 line search
c                    =2 double dogleg
c                    =3 more-hebdon
c udiag(n)     --> diagonal of hessian in a(.,.)
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      dimension x(n),xpls(n),g(n)
      dimension sx(n),sc(n),xplsp(n)
      dimension a(nr,1)
      logical nwtake,mxtake
      dimension udiag(n)
      external fcn

      mxtake=.false.
      do i=1,n
        xpls(i)=x(i)+sc(i)
      end do

      call fcn(n,xpls,fpls)
      dltf=fpls-f
      slp = r8vec_dot ( n, g, sc )
c
c  next statement added for case of compilers which do not optimize
c  evaluation of next "if" statement (in which case fplsp could be
c  undefined).
c
      if(iretcd.eq. 4) fplsp = 0.0D+00
      if(iprint.gt. 1)write(lounit,1100) iretcd,fpls,fplsp,dltf,slp
      if(iretcd.ne. 3 .or. (fpls.lt.fplsp.and.dltf.le. 1.0D-04*slp))
     &                                                     go to 30
c
c  reset xpls to xplsp and terminate global step
c
        iretcd=0

        do i=1,n
          xpls(i)=xplsp(i)
        end do

        fpls=fplsp
        dlt = 0.5D+00 * dlt
        if(iprint.gt. 1)write(lounit,1000)
        go to 180
c     else
c
c  fpls too large
c
   30   if(dltf.le. 1.0D-04 *slp) go to 80
          if(iprint.gt. 1)write(lounit,1010)
          rln = 0.0D+00
          do i=1,n
            rln=max(rln,abs(sc(i))/max(abs(xpls(i)), 1.0D+00 /sx(i)))
          end do
          if(iprint.gt. 1)write(lounit,1110) rln
          if(rln.ge.steptl) go to 50
c         if(rln.lt.steptl)
c         then
c
c  cannot find satisfactory xpls sufficiently distinct from x
c
            iretcd=1
            if(iprint.gt. 1)write(lounit,1030)
            go to 180
c         else
c
c  reduce trust region and continue global step
c
   50       iretcd=2
            dltmp=-slp*dlt/( 2.0D+00 *(dltf-slp))
            if(iprint.gt. 1)write(lounit,1120) dltmp
            if(dltmp.ge. 0.1D+00 * dlt) go to 60
              dlt= 0.1D+00 * dlt
              go to 70
   60         dlt=dltmp
   70       continue
            if(iprint.gt. 1)write(lounit,1040)
            go to 180
c
c  fpls sufficiently small
c
   80     continue
          if(iprint.gt. 1)write(lounit,1070)
          dltfp = 0.0D+00
          if (method .eq. 3) go to 110
          do i = 1, n
            temp = 0.0D+00
            do j = i, n
              temp = temp + (a(j, i)*sc(j))
            end do
            dltfp = dltfp + temp*temp
          end do
          go to 140
c
c         else
c
  110     do 130 i = 1, n
             dltfp = dltfp + udiag(i)*sc(i)*sc(i)
             if (i .eq. n) go to 130
             temp = 0
             ip1 = i + 1
             do j = ip1, n
                temp = temp + a(i, j)*sc(i)*sc(j)
             end do
             dltfp = dltfp + 2.0D+00 * temp
  130     continue
c
c         end if
c
  140     dltfp = slp + dltfp / 2.0D+00
          if(iprint.gt. 1)write(lounit,1130) dltfp,nwtake
          if(iretcd.eq. 2 .or. (abs(dltfp-dltf).gt. 0.1D+00 *abs(dltf))
     &         .or. nwtake .or. (dlt.gt. 0.99D+00 *stepmx)) go to 160
c         if(iretcd.ne. 2 .and. (abs(dltfp-dltf) .le. 0.1D+00 *abs(dltf))
c    &         .and. (.not.nwtake) .and. (dlt.le. 0.99D+00 *stepmx))
c         then
c
c  double trust region and continue global step
c
            iretcd=3
            do i=1,n
              xplsp(i)=xpls(i)
            end do
            fplsp=fpls
            dlt=min ( 2.0D+00 * dlt,stepmx)
            if(iprint.gt. 1)write(lounit,1080)
            go to 180
c         else
c
c  accept xpls as next iterate.  choose new trust region.
c
  160       continue
            if(iprint.gt. 1)write(lounit,1090)
            iretcd=0
            if(dlt.gt. 0.99D+00 *stepmx) mxtake=.true.
            if(dltf.lt. 0.1D+00 *dltfp) go to 170
c           if(dltf.ge. 0.1D+00 *dltfp)
c           then
c
c  decrease trust region for next iteration
c
              dlt= 0.5D+00 * dlt
              go to 180
c           else
c
c  check whether to increase trust region for next iteration
c
  170         if(dltf.le. 0.75D+00 *dltfp) 
     &          dlt=min ( 2.0D+00 *dlt,stepmx)

c           end if
c         end if
c       end if
c     end if
  180 continue
      if(iprint.gt. 1)write(lounit,1020)
      if(iprint.gt. 1)write(lounit,1050) iretcd,mxtake,dlt,fpls
      if(iprint.gt. 1)write(lounit,1060)
      if(iprint.gt. 1)write(lounit,1140) (xpls(i),i=1,n)
      return
 1000 format('0tregup    reset xpls to xplsp. terminate global step')
 1010 format(26h0tregup    fpls too large.)
 1020 format(38h0tregup    values after call to tregup)
 1030 format(54h tregup    cannot find satisfactory xpls distinct from,
     &       27h x.  terminate global step.)
 1040 format(53h tregup    reduce trust region. continue global step.)
 1050 format(21h tregup       iretcd=,i3/
     &       21h tregup       mxtake=,l1/
     &       21h tregup       dlt   =,g16.8/
     &       21h tregup       fpls  =,g16.8)
 1060 format(32h tregup       new iterate (xpls))
 1070 format(35h tregup    fpls sufficiently small.)
 1080 format(54h tregup    double trust region.  continue global step.)
 1090 format(' tregup    accept xpls.  choose new trust region.',
     &       ' terminate global step.')
 1100 format(18h tregup    iretcd=,i5/
     &       18h tregup    fpls  =,g16.8/
     &       18h tregup    fplsp =,g16.8/
     &       18h tregup    dltf  =,g16.8/
     &       18h tregup    slp   =,g16.8)
 1110 format(18h tregup    rln   =,g16.8)
 1120 format(18h tregup    dltmp =,g16.8)
 1130 format(18h tregup    dltfp =,g16.8/
     &       18h tregup    nwtake=,l1)
 1140 format(14h tregup       ,5(g16.8,3x))
      end
