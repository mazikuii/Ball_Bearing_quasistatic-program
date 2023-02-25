      subroutine grdchk_u(n,x,fcn,f,g,typsiz,sx,fscale,rnf,
     &     analtl,wrk1,msg,iprint,lounit)

c*********************************************************************72
c
cc GRDCHK compares the analytic gradient against a finite difference estimate.
c
c  Discussion:
c
c    check analytic gradient against estimated gradient
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
c x(n)         --> estimate to a root of fcn
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                       fcn:  r(n) --> r(1)
c f            --> function value:  fcn(x)
c g(n)         --> gradient:  g(x)
c typsiz(n)    --> typical size for each component of x
c sx(n)        --> diagonal scaling matrix:  sx(i)=1./typsiz(i)
c fscale       --> estimate of scale of objective function fcn
c rnf          --> relative noise in optimization function fcn
c analtl       --> tolerance for comparison of estimated and
c                  analytical gradients
c wrk1(n)      --> workspace
c msg         <--  message or error code
c                    on output: =-21, probable coding error of gradient
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n)
      dimension sx(n),typsiz(n)
      dimension wrk1(n)
      external fcn
c
c compute first order finite difference gradient and compare to
c analytic gradient.
c
      call fstofd_u(1,1,n,x,fcn,f,wrk1,sx,rnf,wrk,1)
c
      ker=0
      do i=1,n
        gs=max(abs(f),fscale)/max(abs(x(i)),typsiz(i))
        if(abs(g(i)-wrk1(i)).gt.max(abs(g(i)),gs)*analtl) then
		  ker=1
		  write(*,*)'grdchk suspects an error in gradient ',i
		  end if
        end do

      if(ker.ne.0)then
        if(iprint.gt. 0)write(lounit,1000)
        if(iprint.gt. 1)write(lounit,1010)
        if(iprint.gt. 1)write(lounit,1020) (i,g(i),wrk1(i),i=1,n)
        msg=-21
		end if

      if(ker.eq.0.and.iprint.gt. 1)write(lounit,1030)
      return
 1000 format(47h grdchk    probable error in coding of analytic,
     &       19h gradient function.)
 1010 format(16h grdchk     comp,12x,8hanalytic,12x,8hestimate)
 1020 format(11h grdchk    ,i5,3x,g16.8,3x,g16.8)
 1030 format(' grdchk    user gradient seems correct')
      end