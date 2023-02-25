	SUBROUTINE trapzd(func,a,b,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_dp*del,del,it)))
		s=0.5_dp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd
