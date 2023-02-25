subroutine neinchk(ftol,xtol,mtol,btol,etak,itnlimit,sx,sy,&
									 restart,maxtaken,&
									 lshmsd,global,method,termcode,&
									 itncount,itconsec,consecmax,nmt,&
									 macheps)
	use nrtype
	implicit none
	real(dp)::ftol,xtol,mtol,btol,etak
	integer::itnlimit
	real(dp)::sx(:),sy(:)
	logical::restart,maxtaken
	integer::lshmsd,global,method,termcode
	integer::itncount,itconsec,consecmax,nmt
	real(dp)::macheps
	
	!1st block
  if(size(sx(:)) .le. 0 .or. size(sy(:)) .le. 0) then
     termcode = -1
     return
  endif
  macheps = epsilon(0.0_dp)
  if(ftol .lt. 0.0_dp) ftol=macheps**(2.0_dp/3.0_dp)
  if(xtol .lt. 0.0_dp) xtol=macheps**(2.0_dp/3.0_dp)
  if(mtol .lt. 0.0_dp) mtol=macheps**(2.0_dp/3.0_dp)
  if(btol .lt. xtol) btol = xtol
	if(etak .lt. 0.0_dp .or. etak .gt. 1.0_dp) etak=.85_dp
  if(itnlimit .le. 0) itnlimit = 150
	where(sx .lt. 0.0_dp) sx=-sx
	where(sy .lt. 0.0_dp) sy=-sy
	
	!2nd block
	if(.not.restart) restart=.true.
	if(maxtaken) maxtaken=.false.

	!3rd block
	if(lshmsd<1 .or. lshmsd>4) lshmsd=1
  if(global .lt. 1 .or. global .gt. 4) global = 1
	if(method .lt.1 .or. method .gt. 2) method=2
  termcode = 0
	itncount=0
	itconsec=0
	if(consecmax<0 .or. consecmax>200) consecmax=5
	if(nmt<=0) nmt=1
  
	!4th block
endsubroutine neinchk