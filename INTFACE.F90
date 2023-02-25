module intface
		interface
			function ei1st(x)
				use fsbr
				implicit none
				real(dp), dimension(:), intent(in) :: x
				real(dp), dimension(size(x)) :: ei1st
			end function ei1st
		end interface

		interface
			function ei2nd(x)
				use fsbr
				implicit none
				real(dp), dimension(:), intent(in) :: x
				real(dp), dimension(size(x)) :: ei2nd
			end function ei2nd
		end interface

		interface
			function firstinte(x)
				use fsbr
				implicit none
				real(dp), dimension(:), intent(in) :: x
				real(dp), dimension(size(x)) :: firstinte
			end function firstinte
		end interface

		interface
			function secondinte(x)
				use fsbr
				implicit none
				real(dp), dimension(:), intent(in) :: x
				real(dp), dimension(size(x)) :: secondinte
			end function secondinte
		end interface

		interface
			function ellip_new(rdiff)
				use fsbr
				implicit none
				real(dp) :: ellip_new,rdiff
			end function ellip_new
		end interface

		interface
			function pow(base,exp)
				use fsbr
				implicit none
				real(dp),intent(in)	:: base,exp
				real(dp)			:: pow
			end function pow		
		end interface

		interface
			function rtsafe(funcd,x1,x2,xacc)
			use nrtype; use nrutil, only : nrerror
			use fsbr
			implicit none
			real(dp), intent(in) :: x1,x2,xacc
			real(dp) :: rtsafe
			interface
				subroutine funcd(x,fval,fderiv)
				use nrtype
				use fsbr
				implicit none
				real(dp), intent(in) :: x
				real(dp), intent(out) :: fval,fderiv
				end subroutine funcd
			end interface
			end function rtsafe
		endinterface

		interface
			subroutine funcd(x,fval,fderiv)
				use fsbr
				use nrtype
				implicit none
				real(dp), intent(in) :: x
				real(dp), intent(out) :: fval,fderiv
				real(dp) ::fa,a,kn,alpha0
				integer	 ::z
			endsubroutine funcd
		endinterface

		interface
			function rad2deg(rad)
				use fsbr
				implicit none
				real(dp),intent(in)	:: rad
				real(dp)			:: rad2deg
			end function rad2deg
			function deg2rad(deg)
				use nrtype
				implicit none
				real(dp),intent(in)	:: deg
				real(dp)			:: deg2rad
			end function deg2rad
		endinterface

		interface
			function rev2rad(rev)
				use nrtype
				implicit none
				real(dp),intent(in)	:: rev
				real(dp)			:: rev2rad
			end function rev2rad		
		endinterface

		interface
			function rad2rev(rev)
				use nrtype
				implicit none
				real(dp),intent(in)	:: rev
				real(dp)			:: rad2rev
			end function rad2rev
		endinterface

		interface
			function fv(x)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(size(x)) :: fv
			end function fv
		end interface

		interface
			function funcv(x)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(size(x)) :: funcv
			end function funcv
		end interface

		interface
			subroutine frosen ( x, f, m, n )
				use fsbr
				use nrtype
				implicit none
				integer   ::m, n
				real(dp)	::x(n), f(m)
			end
		end interface

		interface
			subroutine jwood ( x, jac, maxm, m, n )
				integer            maxm, m, n
				double precision   x(n), jac(maxm,n)
			end
		end interface

		interface
			subroutine usrfun(x,fvec,fjac)
			use nrtype; use nrutil
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(:), intent(out) :: fvec
			real(dp), dimension(:,:), intent(out) :: fjac
			end subroutine usrfun
		end interface

		interface
			subroutine fun(n,x,fvec,iflag)
				use nrtype; use nrutil
				use nr
				implicit none
				integer			::n,iflag
				real(dp), dimension(n) :: x,fvec
			end subroutine fun
			subroutine FUN_asym(n,x,fvec,iflag)
				use nrtype; use nrutil
				use nr
				implicit none
				integer			::n,iflag
				real(dp), dimension(n) :: x,fvec
			end subroutine FUN_asym
		end interface

		interface
			subroutine fun10(x,f,n)
				use nrtype
				use nr
				implicit none
				integer			::n
				real(dp) :: x(n),f(n)
			end subroutine fun10
		end interface

		interface
			subroutine fun12(n,x,f)
				use nrtype
				use nr
				implicit none
				integer			::n
				real(dp) :: x(n),f(n)
			end subroutine fun12
		end interface

		interface
			subroutine dfun12(n, x, fjac,ldfjac)		
				use nrtype
				use nr
				implicit none
				integer			::n,ldfjac
				real(dp) :: x(n),fjac(n,n)
			end subroutine dfun12
		end interface

		interface
			subroutine lsjac(n, x, fjac)		
				use nrtype
				use nr
				implicit none
				integer			::n
				real(dp) :: x(n),fjac(n,n),fvec(n)
			end subroutine lsjac
		end interface

		interface
			subroutine f6eq(n, xcur, fcur, rpar, ipar, itrmf)
				use nrtype
				implicit none
				integer i, itrmf, j, j1, j2, n, ipar(*), nx, ny 
				real(dp)	:: cl, cr, h2l, xcur(n), fcur(n), rpar(*)
			endsubroutine f6eq
		end interface

		interface
			subroutine j6eq(n,xcur,fcur,ijob,v,z,rpar,ipar,itrmjv)
				implicit none
				integer n, ijob, ipar(*), itrmjv
				double precision xcur(n), fcur(n), v(n), z(n), rpar(*)
			endsubroutine j6eq
		end interface

		interface
			function func(p)
			use nrtype
			real(dp), dimension(:), intent(in) :: p
			real(dp) :: func
			end function func
		end interface

		interface
			function dfunc(p)
			use nrtype
			real(dp), dimension(:), intent(in) :: p
			real(dp), dimension(size(p)) :: dfunc
			end function dfunc
		end interface

		interface
			subroutine slv24_neslv()
				use	nrtype
				use neintface
				USE FSBR
				implicit none
					
			endsubroutine slv24_neslv
		endinterface
endmodule intface	