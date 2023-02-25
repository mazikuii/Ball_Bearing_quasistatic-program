
      subroutine heschk_u(nr,n,x,fcn,d1fcn_u,d2fcn_u,f,g,a,typsiz,
     &     sx,rnf,analtl,iagflg,udiag,wrk1,wrk2,msg,iprint,lounit)

c*********************************************************************72
c
cc HESCHK compares the analytic Hessian against a finite difference estimate.
c
c  Discussion:
c
c    check analytic hessian against estimated hessian
c    (this may be done only if the user supplied analytic hessian
c    d2fcn fills only the lower triangular part and diagonal of a)
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
c x(n)         --> estimate to a root of fcn
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                       fcn:  r(n) --> r(1)
c d1fcn        --> name of subroutine to evaluate gradient of fcn.
c                  must be declared external in calling routine
c d2fcn        --> name of subroutine to evaluate hessian of fcn.
c                  must be declared external in calling routine
c f            --> function value:  fcn(x)
c g(n)        <--  gradient:  g(x)
c a(n,n)      <--  on exit:  hessian in lower triangular part and diag
c typsiz(n)    --> typical size for each component of x
c sx(n)        --> diagonal scaling matrix:  sx(i)=1./typsiz(i)
c rnf          --> relative noise in optimization function fcn
c analtl       --> tolerance for comparison of estimated and
c                  analytical gradients
c iagflg       --> =1 if analytic gradient supplied
c udiag(n)     --> workspace
c wrk1(n)      --> workspace
c wrk2(n)      --> workspace
c msg         <--> message or error code
c                    on input : if =1xx do not compare anal + est hess
c                    on output: =-22, probable coding error of hessian
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n),a(nr,1)
      dimension typsiz(n),sx(n)
      dimension udiag(n),wrk1(n),wrk2(n)
      external fcn
      external d1fcn_u
      external d2fcn_u
c
c compute finite difference approximation a to the hessian.
c
      if(iagflg.eq. 1) call fstofd_u(nr,n,n,x,d1fcn_u,g,a,sx,rnf,wrk1,3)
      if(iagflg.ne. 1) call sndofd_u(nr,n,x,fcn,f,a,sx,rnf,wrk1,wrk2)
c
      ker=0
c
c copy lower triangular part of "a" to upper triangular part
c and diagonal of "a" to udiag
c
      do 20 j=1,n
        udiag(j)=a(j,j)
        if(j.eq.n) go to 20
        jp1=j+1
        do i=jp1,n
          a(j,i)=a(i,j)
        end do
   20 continue
c
c compute analytic hessian and compare to finite difference
c approximation.
c
      call d2fcn_u(nr,n,x,a)
      do 40 j=1,n
        hs=max(abs(g(j)), 1.0D+00 )/max(abs(x(j)),typsiz(j))
        if(abs(a(j,j)-udiag(j)).gt.max(abs(udiag(j)),hs)*analtl)
     &       ker=1
        if(j.eq.n) go to 40
        jp1=j+1
        do i=jp1,n
          if(abs(a(i,j)-a(j,i)).gt.max(abs(a(i,j)),hs)*analtl) ker=1
        end do
   40 continue

      if(ker.eq.0) go to 80
        if(iprint.gt. 0)write(lounit,1000)
        if(iprint.gt. 1)write(lounit,1010)
        do i=1,n
          if(i.eq. 1) go to 60
          im1=i-1
          do j=1,im1
            if(iprint.gt. 1)write(lounit,1020) i,j,a(i,j),a(j,i)
          end do
   60     if(iprint.gt. 1)write(lounit,1020) i,i,a(i,i),udiag(i)
        end do
        msg=-22

   80 continue
      if(ker.eq.0.and.iprint.gt. 1)write(lounit,1030)
      return
 1000 format(47h heschk    probable error in coding of analytic,
     &       18h hessian function.)
 1010 format(21h heschk      row  col,14x,8hanalytic,14x,10h(estimate))
 1020 format(11h heschk    ,2i5,2x,g16.8,2x,1h(,g16.8,1h))
 1030 format('0heschk    user hessian seems correct')
      end
