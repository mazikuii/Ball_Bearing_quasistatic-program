	SUBROUTINE SLV10_DNEQNF()
	USE FSBR
	USE nrtype; USE nrutil
	USE nr
	USE IMSLF90
	USE INTFACE,ONLY:FUN10
	IMPLICIT NONE
	REAL(DP), PARAMETER :: ERRREL=1.0e-13_DP
	INTEGER(I4B), PARAMETER :: N=6,ITMAX=500
	REAL(DP),DIMENSION(N)		:: XGUESS,X
	REAL(DP)								::	FNORM

	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************     BY DNEQNF SUBROUTINE   ****************'
	WRITE(*,*)'         ***********************************************************'

	CALL X0VCT(XGUESS)
	CALL PRINTX0(XGUESS)
	CALL DNEQNF(FUN10,ERRREL,N,ITMAX,XGUESS,X,FNORM)

	CALL SLV6EQS(X,XGUESS,.false.)
	CALL SAVEXf(X,XGUESS)
	CALL PRINTXF(X,XGUESS)
!	CALL WRTTM()  !输出拖动力、拖动力矩
	ENDSUBROUTINE SLV10_DNEQNF

	SUBROUTINE FUN10(X,F,N)
		USE NRTYPE
		USE NR
		IMPLICIT NONE
		INTEGER			::N
		REAL(DP) :: X(N),F(N)

		CALL SLV6EQS(X,F,.false.)
	END SUBROUTINE FUN10
