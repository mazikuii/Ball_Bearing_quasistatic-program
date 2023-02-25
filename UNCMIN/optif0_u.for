      subroutine optif0_u(nr,n,x,fcn,xpls,fpls,gpls,itrmcd,a,wrk)

c*********************************************************************72
c
cc OPTIF0 provides the simplest interface to the optimization package.
c
c  Discussion:
c
c    provide simplest interface to minimization package.
c    user has no control over options.
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
c                  a positive integer specifying the row dimension
c                  of the matrices a and wrk in the user's calling
c                  program.  nr must be greater or equal to n.
c
c n            --> dimension of problem
c                  a positive integer specifying the order or dimension
c                  of the problem.  the program will abort if
c                  n.le.0.  the program is inefficient for n.eq.1.
c                  for n.eq.1, do not call this routine, but rather
c                  optif9.
c
c x(n)        <--> initial estimate of minimum on input.
c                  on return, x contains x(k-1), the next to last
c                  iterate.
c
c fcn          --> name of routine to evaluate minimization function.
c                  fcn must be declared external in the calling
c                  routine, and must have the form
c                  subroutine fcn(n,x,f) where x is a vector of
c                  length n.  the subroutine must not alter x or n.
c
c xpls(n)     <--  local minimum
c                  an n dimensional array containing the best
c                  approximation to the local minimum on return
c                  if the algorithm has converged.
c
c fpls        <--  function value at local minimum xpls
c
c gpls(n)     <--  gradient at local minimum xpls
c
c itrmcd      <--  termination code
c                  0=bad input (see msg).
c                  1=relative gradient is close to zero.  current
c                  iterate is probably solution.
c                  2=successive iterates within tolerance. current
c                  iterate is probably solution.
c                  3=last global step failed to locate a point
c                  lower than xpls.  either xpls is an
c                  approximate local minimum of the function,
c                  the function is too nonlinear for this program,
c                  or steptl is too large.
c                  4=iteration limit exceeded.
c                  5=maximum step size stepmx exceeded 5 consecutive
c                  times.  either the function is unbounded below,
c                  becomes asymptotic to a finite value from above
c                  in some direction, or stepmx is too small.
c
c a(n,n)       --> workspace
c                  a matrix used to store the hessian and its
c                  cholesky decomposition.  the user must declare
c                  this array of dimensions nr by nc where nr
c                  and nc are both at least n.  the hessian is
c                  not evaluated at the last iterate so on return
c                  the information in a is old.
c
c wrk(n,9)     --> workspace
c                  workspace required by the subroutine.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),xpls(n),gpls(n)
      dimension a(nr,1),wrk(nr,1)
      external fcn,d1fcn_u,d2fcn_u
c
c equivalence wrk(n,1) = udiag(n)
c             wrk(n,2) = g(n)
c             wrk(n,3) = p(n)
c             wrk(n,4) = typsiz(n)
c             wrk(n,5) = sx(n)
c             wrk(n,6) = wrk0(n)
c             wrk(n,7) = wrk1(n)
c             wrk(n,8) = wrk2(n)
c             wrk(n,9) = wrk3(n)
c
      call dfault_u(n,wrk(1,4),fscale,method,iexp,msg,ndigit,
     &     itnlim,iagflg,iahflg,iprint,lounit,dlt,gradtl,stepmx,steptl)

      call optdrv_u(nr,n,x,fcn,d1fcn_u,d2fcn_u,wrk(1,4),fscale,
     &     method,iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,
     &     dlt,gradtl,stepmx,steptl,
     &     xpls,fpls,gpls,itrmcd,
     &     a,wrk(1,1),wrk(1,2),wrk(1,3),wrk(1,5),wrk(1,6),
     &     wrk(1,7),wrk(1,8),wrk(1,9))

      return
      end
