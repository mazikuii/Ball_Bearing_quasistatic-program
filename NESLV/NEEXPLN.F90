subroutine neexpln ( itrmcd,consecmax)
	implicit none
  integer itrmcd,consecmax

  write(*,*)'NESLV: //////////////////////////////////////////////////'
  write(*,*)'NESLV: EXPLAIN:'
  write(*,200)'NESLV:   NESLV returned a termination code of', itrmcd
  write(*,*)'NESLV:  '
  write(*,*)'NESLV: This code has the following meaning:'
  write(*,*)'NESLV:  '

  if(itrmcd.eq.0)then
    write(*,*)'NESLV: This termination code means NESLV received'
    write(*,*)'NESLV: illegal or unacceptable input.'
  elseif(itrmcd.eq. 1)then
    write(*,*)'NESLV: The gradient is relatively close to zero.'
    write(*,*)'NESLV: The current iterate is probably a good solution.'
  elseif(itrmcd.eq. 2)then
    write(*,*)'NESLV: Successive iterates are not changing much.'
    write(*,*)'NESLV: The current iterate is probably a good solution.'
  elseif(itrmcd.eq. 3)then
    write(*,*)'NESLV: The last global step could not find a point with '
    write(*,*)'NESLV: lower function value than that achieved by XPLS.'
    write(*,*)'NESLV:  '
    write(*,*)'NESLV: XPLS may be an approximate local minimum, or '
    write(*,*)'NESLV:  '
    write(*,*)'NESLV: The function is too nonlinear for this program, or'
    write(*,*)'NESLV: the stepsize tolerance is too large.'
  elseif(itrmcd.eq. 4)then
    write(*,*)'NESLV: The iteration limit was exceeded.'
    write(*,*)'NESLV: The current iterate does not satisfy the '
    write (*,*)'NESLV: requirements.'
  elseif(itrmcd.eq. 5)then
    write(*,100)'NESLV: maximum stepsize was exceeded', consecmax,'consecutive times.'
    write(*,*)'NESLV:  '
    write(*,*)'NESLV: Either the function is unbounded below, or '
    write(*,*)'NESLV: becomes asymptotic to a finite value from above in'
    write(*,*)'NESLV: some direction, or'
    write(*,*)'NESLV: the value of STEPMX is simply too small.'
  else
    write(*,*)'NESLV: This is an unknown error code!'
  end if

  write(*,*)'NESLV: /////////////////////////////////////////////////'
	100 format(1x,a36,1x,i3,1x,a18)
	200 format(1x,a45,1x,i1)
endsubroutine neexpln
