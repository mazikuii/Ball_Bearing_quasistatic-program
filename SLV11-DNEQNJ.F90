	SUBROUTINE SLV11_DNEQNJ()
	USE FSBR
	USE nrtype; USE nrutil
	USE nr
	USE INTFACE,ONLY:FUN10,LSJAC
	IMPLICIT NONE
	REAL(DP), PARAMETER :: ERRREL=1.0e-13_DP
	INTEGER(I4B), PARAMETER :: N=6,ITMAX=100
	REAL(DP),DIMENSION(N)		:: XGUESS,X
	REAL(DP)								::	FNORM

	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************     BY DNEQNJ SUBROUTINE   ****************'
	WRITE(*,*)'         ***********************************************************'

	CALL X0VCT(XGUESS)
	CALL PRINTX0(XGUESS)
	CALL DNEQNJ(FUN10,LSJAC,ERRREL,N,ITMAX,XGUESS,X,FNORM)

	CALL SLV6EQS(X,XGUESS,.false.)
	CALL SAVEXf(X,XGUESS)
	WRITE(*,*)	"FNORM=",FNORM
	CALL PRINTXF(X,XGUESS)
!	CALL WRTTM()  !输出拖动力、拖动力矩
	ENDSUBROUTINE SLV11_DNEQNJ

	SUBROUTINE LSJAC(N, X, FJAC)		
		USE NRTYPE
		USE NR
		IMPLICIT NONE
		INTEGER			::N
		REAL(DP) :: X(N),FJAC(N,N),FVEC(N)

		CALL SLV6EQS(X,FVEC,.false.)
		CALL FDJAC(X,FVEC,FJAC)
	END SUBROUTINE LSJAC
