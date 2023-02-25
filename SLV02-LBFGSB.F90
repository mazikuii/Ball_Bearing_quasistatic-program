!测试程序，对于本程序而言L-BFGS-B算法并不比BRODYN算法或NEWTON算法优秀
SUBROUTINE SLV02_LBFGSB()
	USE FSBR
	USE INTFACE,ONLY:FV
	USE nr, ONLY : FDJAC_D
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER		:: NMAX=1024,MMAX=17
	INTEGER(I4B), PARAMETER		:: N=6 !方程个数
  CHARACTER*60     TASK, CSAVE
  LOGICAL      LSAVE(4)
  INTEGER	::  M, IPRINT,NBD(NMAX), IWA(3*NMAX), ISAVE(44)
  REAL(DP)::	F, FACTR, PGTOL,X(NMAX), L(NMAX), U(NMAX), G(NMAX), DSAVE(29) 
  REAL(DP)::   WA(2*MMAX*NMAX+4*NMAX+12*MMAX*MMAX+12*MMAX)
	!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(DP), DIMENSION(N)			:: XV,FVC
	REAL(DP), DIMENSION(N,N)		:: DF
	INTEGER(I4B)								:: I

	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            *************** BY L-BFGS-B METHOD  ***************'
	WRITE(*,*)'            ***************************************************'

	task = 'START'      !On first entry, it must be set to 'START'
	iprint = -1		!output control
  factr=1.0E-12_DP	!tolerances in the stopping criteria.
  pgtol=1.0E-12_DP	!tolerances in the stopping criteria.
  m=50

	CALL X0VCT(X)
	CALL PRINTX0(X)
	DO I=1,N
		nbd(I)=2		!2 if x(i) has both lower and upper bounds
	ENDDO
	l(1)=60.0_DP;			u(1)= 100.0_DP			
	l(2)=-5.0_DP;			u(2)= 5.0_DP			
	l(3)=10.0_DP;			u(3)= 300._DP			
	l(4)=200.0_DP;		u(4)= 540.0_DP			
	l(5)=0.60_DP;			u(5)= 1.50_DP			
	l(6)=0.60_DP;			u(6)= 1.50_DP	
!!!!!!!!!!!!!!!!!!!
	DO WHILE(.TRUE.)
		CALL SETULB(N,M,X,L,U,NBD,F,G,FACTR,PGTOL,WA,IWA,TASK,IPRINT,CSAVE,LSAVE,ISAVE,DSAVE)
		IF(TASK(1:2) .EQ. 'FG') THEN
			DO I=1,N
				XV(I)=X(I)
			ENDDO
			FVC=FV(XV)
			F=SUM(FVC)
			CALL FDJAC_D(XV,FVC,DF)
			DO I=1,N
				G(I)=SUM(DF(I,:))
			ENDDO
			CYCLE
		ENDIF
		IF(TASK(1:5) .EQ. 'NEW_X') THEN
         if (isave(34) .ge. 99) task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
			CYCLE
		ELSE
			EXIT
		ENDIF
	ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL SLV6EQS(X,FVC,.false.)
	CALL SAVEXf(X,FVC)
	CALL PRINTXF(X,FVC)
END SUBROUTINE SLV02_LBFGSB

FUNCTION FV(X)
	USE FSBR
	USE NRTYPE
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: X
	REAL(DP), DIMENSION(SIZE(X)) :: FV
	CHARACTER(LEN=1)			:: TAB=CHAR(9)
	REAL(DP)	::F(SIZE(X))
	INTEGER		::I

	CALL SLV6EQS(X,F,.false.)
	DO I=1,SIZE(X)
		!为方程加平方
		FV(I)=F(I)**2
	ENDDO
ENDFUNCTION FV

