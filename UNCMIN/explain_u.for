      subroutine explain_u ( itrmcd )

c*********************************************************************72
c
cc EXPLAIN prints an explanation for the UNCMIN termination code.
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
c    Input, integer ITRMCD, the UNCMIN termination code.
c
      implicit double precision ( a-h, o-z )

      integer itrmcd

      write(*,*)' '
      write(*,*)'EXPLAIN:'
      write(*,*)'  UNCMIN returned a termination code of ', itrmcd
      write(*,*)' '
      write(*,*)'This code has the following meaning:'
      write(*,*)' '

      if(itrmcd.eq.0)then

        write(*,*)'This termination code means UNCMIN received'
        write(*,*)'illegal or unacceptable input.'
      elseif(itrmcd.eq. 1)then
        write(*,*)'The gradient is relatively close to zero.'
        write(*,*)'The current iterate is probably a good solution.'
      elseif(itrmcd.eq. 2)then
        write(*,*)'Successive iterates are not changing much.'
        write(*,*)'The current iterate is probably a good solution.'
      elseif(itrmcd.eq. 3)then
        write(*,*)'The last global step could not find a point with '
        write(*,*)'lower function value than that achieved by XPLS.'
        write(*,*)' '
        write(*,*)'XPLS may be an approximate local minimum, or '
        write(*,*)' '
        write(*,*)'The function is too nonlinear for this program, or'
        write(*,*)'the stepsize tolerance is too large.'
      elseif(itrmcd.eq. 4)then
        write(*,*)'The iteration limit was exceeded.'
        write(*,*)'The current iterate does not satisfy the '
        write (*,*)'requirements.'
      elseif(itrmcd.eq. 5)then
        write(*,*)'maximum stepsize was exceeded 5 consecutive times.'
        write(*,*)' '
        write(*,*)'Either the function is unbounded below, or '
        write(*,*)'becomes asymptotic to a finite value from above in'
        write(*,*)'some direction, or'
        write(*,*)'the value of STEPMX is simply too small.'
      else
        write(*,*)'This is an unknown error code!'
      end if

      write(*,*)' '

      return
      end
