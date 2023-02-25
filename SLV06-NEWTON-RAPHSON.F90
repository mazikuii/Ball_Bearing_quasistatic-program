	SUBROUTINE SLV06_NEWTON_RAPHSON()
		USE FSBR
		USE nrtype; USE nrutil
		USE nr
		USE INTFACE,ONLY:USRFUN
		IMPLICIT NONE
		REAL(DP), PARAMETER :: TOLX=1.0e-15_DP,TOLF=1.0e-15_DP
		INTEGER(I4B), PARAMETER :: NTRIAL=50,N=6
		INTEGER(I4B) :: j
		REAL(DP), DIMENSION(N) :: fvec,x
		REAL(DP), DIMENSION(N,N) :: fjac

		WRITE(*,*)'         ***********************************************************'
		WRITE(*,*)'         ***************  BY NEWTON-RAPHSON METHOD  ****************'
		WRITE(*,*)'         ***********************************************************'

		CALL X0VCT(X)
		CALL PRINTX0(X)

		CALL MNEWT(NTRIAL,X,TOLX,TOLF,USRFUN)

		CALL USRFUN(X,FVEC,FJAC)
		CALL SAVEXf(X,FVEC)
		CALL PRINTXF(X,FVEC)
	!	CALL WRTTM()  !输出拖动力、拖动力矩
	ENDSUBROUTINE SLV06_NEWTON_RAPHSON

	SUBROUTINE usrfun(x,fvec,fjac)
		USE nrtype; USE nrutil
		USE nr
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
		CALL SLV6EQS(X,FVEC,.false.)
		CALL FDJAC(X,FVEC,FJAC)
	END SUBROUTINE usrfun
