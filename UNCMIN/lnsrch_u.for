      subroutine lnsrch_u(n,x,f,g,p,xpls,fpls,fcn,mxtake,
     &   iretcd,stepmx,steptl,sx,iprint,lounit)

c*********************************************************************72
c
cc LNSRCH finds the next Newton iterate by line search.
c
c  Discussion:
c
c    find a next newton iterate by line search.
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
c x(n)         --> old iterate:   x(k-1)
c f            --> function value at old iterate, f(x)
c g(n)         --> gradient at old iterate, g(x), or approximate
c p(n)         --> non-zero newton step
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c fcn          --> name of subroutine to evaluate function
c iretcd      <--  return code
c mxtake      <--  boolean flag indicating step of maximum length used
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c sx(n)        --> diagonal scaling matrix for x
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c
c internal variables
c ------------------
c sln              newton length
c rln              relative length of newton step
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      integer n,iretcd
      dimension sx(n)
      dimension x(n),g(n),p(n)
      dimension xpls(n)
      logical mxtake
      external fcn

      mxtake=.false.
      iretcd=2
      if(iprint.gt. 1)write(lounit,1040)
      if(iprint.gt. 1)write(lounit,1060) (p(i),i=1,n)

      tmp = 0.0D+00
      do i=1,n
        tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
      end do

      sln=sqrt(tmp)
      if(sln.le.stepmx) go to 20
c
c newton step longer than maximum allowed
c
        scl=stepmx/sln
        call sclmul_u(n,scl,p,p)
        sln=stepmx
        if(iprint.gt. 1)write(lounit,1050)
        if(iprint.gt. 1)write(lounit,1060) (p(i),i=1,n)
   20 continue

      slp = r8vec_dot ( n, g, p )

      rln = 0.0D+00
      do i=1,n
        rln=max(rln,abs(p(i))/max(abs(x(i)), 1.0D+00 /sx(i)))
      end do

      rmnlmb=steptl/rln
      almbda = 1.0D+00
      if(iprint.gt. 1)write(lounit,1020) sln,slp,rmnlmb
c
c loop
c check if new iterate satisfactory.  generate new lambda if necessary.
c
   40 continue
      if(iretcd.lt. 2) return
      do i=1,n
        xpls(i)=x(i) + almbda*p(i)
      end do

      call fcn(n,xpls,fpls)
      if(iprint.gt. 1)write(lounit,1000) almbda
      if(iprint.gt. 1)write(lounit,1010)
      if(iprint.gt. 1)write(lounit,1060) (xpls(i),i=1,n)
      if(iprint.gt. 1)write(lounit,1030) fpls
      if(fpls.gt. f+slp * 1.0D-04 * almbda) go to 60
c
c solution found
c
        iretcd=0
        if(almbda.eq. 1.0D+00 .and. sln.gt. 0.99D+00 * stepmx) then
          mxtake=.true.
        end if
        go to 40
c
c solution not (yet) found
c
   60   if(almbda .ge. rmnlmb) go to 70
c
c no satisfactory xpls found sufficiently distinct from x
c
          iretcd=1
          go to 40
c
c calculate new lambda
c
   70     if(almbda.ne. 1.0D+00 ) go to 80
c
c first backtrack: quadratic fit
c
            tlmbda=-slp/( 2.0D+00 *(fpls-f-slp))
            go to 110
c
c all subsequent backtracks: cubic fit
c
   80       t1=fpls-f-almbda*slp
            t2=pfpls-f-plmbda*slp
            t3 = 1.0D+00 /(almbda-plmbda)
            a=t3*(t1/(almbda*almbda) - t2/(plmbda*plmbda))
            b=t3*(t2*almbda/(plmbda*plmbda)
     &           - t1*plmbda/(almbda*almbda) )
            disc=b*b- 3.0D+00 * a*slp
            if(disc.le. b*b) go to 90
c
c only one positive critical point, must be minimum
c
              tlmbda=(-b+sign ( 1.0D+00, a)*sqrt(disc)) 
     &          / ( 3.0D+00 * a )
              go to 100
c
c both critical points positive, first is minimum
c
   90         tlmbda=(-b-sign ( 1.0D+00, a )*sqrt(disc)) 
     &          / ( 3.0D+00 * a )
  100       if(tlmbda.gt. 0.5D+00 * almbda) tlmbda= 0.5D+00 *almbda
  110     plmbda=almbda
          pfpls=fpls

          if(tlmbda.ge. almbda * 0.1D+00 ) go to 120
            almbda=almbda * 0.1D+00
            go to 130
  120       almbda=tlmbda
  130 go to 40
 
 1000 format(18h lnsrch    almbda=,g16.8)
 1010 format(29h lnsrch    new iterate (xpls))
 1020 format(' lnsrch    newton length =',g16.8/
     &       ' lnsrch    slp =          ',g16.8/
     &       ' lnsrch    rmnlmb   =     ',g16.8)
 1030 format(19h lnsrch    f(xpls)=,g16.8)
 1040 format(26h0lnsrch    newton step (p))
 1050 format(' lnsrch      reduced newton step')
 1060 format(14h lnsrch       ,5(g16.8,3x))
      end
