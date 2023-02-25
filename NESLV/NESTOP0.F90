subroutine nestop0(fvec,sy,ftol,itconsec,termcode)
	use nrtype
	implicit none

	real(dp),dimension(:)::fvec,sy
	real(dp)::ftol
	integer ::itconsec,termcode

	itconsec=0
	if(maxval(abs(fvec(:))*sy(:))<=0.01_dp*ftol) then
		termcode=1	!made
	else
		termcode=0	!failed
	endif
endsubroutine nestop0

