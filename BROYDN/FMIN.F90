MODULE fminln
	USE nrtype; USE nrutil, ONLY : nrerror
	REAL(DP), DIMENSION(:), POINTER :: fmin_fvecp
CONTAINS
!BL
	FUNCTION fmin(x)
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP) :: fmin
	INTERFACE
		FUNCTION funcv(x)
		USE FSBR
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	if (.not. associated(fmin_fvecp)) call &
		nrerror('fmin: problem with pointer for returned values')

	fmin_fvecp=funcv(x)
	fmin=0.5_Dp*dot_product(fmin_fvecp,fmin_fvecp)
	END FUNCTION fmin
END MODULE fminln
