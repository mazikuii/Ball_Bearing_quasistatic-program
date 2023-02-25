SUBROUTINE SLV15_ALCON1()
	USE FSBR
	USE nrtype
	USE INTFACE,ONLY:FUNCV
	IMPLICIT NONE
	EXTERNAL FT
	INTEGER,PARAMETER	::N=6
	INTEGER						:: UPR,UDIAG
	INTEGER						:: LRW,LIW
	INTEGER						:: INFO(9)
	REAL(DP)					:: X(N),XW(N),RWORK(3*N**2+25*N+19),IWORK(N+1)
	REAL(DP)					:: TAU,TAUMIN,TAUMAX
	REAL(DP)					:: EPS
	COMMON /UNIT/ UPR,UDIAG

	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************     BY ALCON2 SUBROUTINE   ****************'
	WRITE(*,*)'         ***********************************************************'

	UDIAG=2
	LRW=3*N**2+25*N+19
	LIW=N+1

	CALL X0VCT(X)
	XW(:)=1.0_DP

	TAU=BR%MU
	TAUMIN=1.0E-2_DP
	TAUMAX=BR%MU

	EPS=1.0E-3_DP
	INFO(1)=0
	INFO(2)=99
	INFO(3)=-1
	INFO(4)=1
	CALL ALCON1(FT,N,X,XW,TAU,TAUMIN,TAUMAX,EPS,INFO,RWORK,LRW,IWORK,LIW)

	CALL SLV6EQS(X,XW,.false.)
	CALL PRINTXF(X,XW)
ENDSUBROUTINE SLV15_ALCON1

