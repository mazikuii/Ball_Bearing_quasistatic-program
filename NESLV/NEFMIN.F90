function nefmin(x,func,fvec,sy)
	use nrtype
	implicit none
	real(dp),dimension(:), intent(in) :: x,sy
	real(dp),dimension(:)::fvec
	real(dp)::nefmin
	interface
		function func(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	func(size(x(:)))
		endfunction func	
	end interface

	fvec=func(x)
	nefmin=0.5_dp*dot_product(fvec*sy(:),fvec*sy(:))
!	nefmin=0.5_dp*vabs(fvec*sy(:))**2
!	nefmin=0.5_dp*ddot(size(x),fvec*sy(:),1,fvec*sy(:),1)
end function nefmin
