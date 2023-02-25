      subroutine result_u(nr,n,x,f,g,a,p,itncnt,iflg)

c*********************************************************************72
c
cc RESULT prints information from the optimization procedure.
c
c  Discussion:
c
c    print information
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
c x(n)         --> iterate x(k)
c f            --> function value at x(k)
c g(n)         --> gradient at x(k)
c a(n,n)       --> hessian at x(k)
c p(n)         --> step taken
c itncnt       --> iteration number k
c iflg         --> flag controlling info to print
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)
      dimension g(n)
      dimension p(n)
      dimension x(n)

      write(*,*)' '
      write(*,*)'RESULT - Iterate ',itncnt
      if(iflg.ne.0)then
        write(*,*)' '
        write(*,*)'Step'
        write(*,*)' '
        write(*,'(5g14.6)') (p(i),i=1,n)
      end if
      write(*,*)' '
      write(*,*)'X'
      write(*,*)' '
      write(*,'(5g14.6)') (x(i),i=1,n)
      write(*,*)'Function value = ',f
      write(*,*)' '
      write(*,*)'Gradient vector'
      write(*,*)' '
      write(*,'(5g14.6)') (g(i),i=1,n)
      if(iflg.ne.0)then
        write(*,*)' '
        write(*,*)'Hessian matrix'
        write(*,*)' '
        do i=1,n
          write(*,'(5g14.6)') (a(i,j),j=1,i),(a(k,i),k=i+1,n)
        end do
      end if

      return
      end
