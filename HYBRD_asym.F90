subroutine hybrd_asym(x)
	use fsbr
	use nrtype; use nrutil
	use nr
	use intface,only:fun_asym
	implicit none
	real(dp), parameter :: tol=1.0e-16_dp
	integer(i4b), parameter :: n=2
	integer(i4b) :: info,lwa
	real(dp), dimension(n) :: fvec,x
	real(dp),allocatable  :: wa(:)


	lwa=(n*(3*n+13))/2+1
	allocate(wa(lwa))
	call hybrd1(fun_asym,n,x,fvec,tol,info,wa,lwa)
endsubroutine hybrd_asym

subroutine fun_asym(n,x,fvec,iflag)
	use nrtype; use nrutil
	use nr
	use neintface,only:slv2eqs_asym
	implicit none
	integer			::n,iflag
	real(dp), dimension(n) :: x,fvec

	call slv2eqs_asym(x,fvec,.false.)
end subroutine fun_asym
