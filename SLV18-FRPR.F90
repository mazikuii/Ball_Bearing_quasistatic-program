	SUBROUTINE SLV18_FRPRMN()
!	driver for routine frprmn
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=6
	REAL(DP), PARAMETER :: FTOL=1.0e-10_Dp
	INTEGER(I4B) :: iter
	REAL(DP) :: fret
	REAL(DP), DIMENSION(NDIM) :: p,F
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP) :: func
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
!BL
		FUNCTION dfunc(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE


	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************  BY FRPRMN METHOD  ****************'
	WRITE(*,*)'            ***************************************************'

	CALL X0VCT(P)
	CALL PRINTX0(P)					!µü´ú³õÖµ
	call frprmn(p,FTOL,iter,fret)

	write(*,'(1x,a,i3)') 'Iterations:',iter
	write(*,'(1x,a,e14.6)') 'Func. value at solution',fret
	CALL SLV6EQS(P,F,.false.)
	CALL PRINTXF(P,F)
	F=DFUNC(P)
	WRITE(*,*) 'CONVERGENCE GRADIENCE'
	CALL PRINTXF(P,F)

	ENDSUBROUTINE SLV18_FRPRMN

