SUBROUTINE SLV01_BROYDN()
	USE FSBR
	USE NR
	USE INTFACE,ONLY:FUNCV
	IMPLICIT NONE
	INTEGER(I4B) :: k
	REAL(DP), DIMENSION(BR%M) :: F
	REAL(DP), DIMENSION(BR%N) :: X
	LOGICAL(LGT) :: CHECK
	CHARACTER(LEN=1)			:: TAB=CHAR(9)
	REAL(DP)	::MU0
	REAL(DP)	::XW(6)

	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************  BY BRODYN METHOD  ****************'
	WRITE(*,*)'            ***************************************************'

	CALL X0VCT(X)
	CALL PRINTX0(X)					!迭代初值

	WRITE(*,*) '*********   MU   ******'
	WRITE(*,*) '        ',BR%MU
	WRITE(*,*) '***********************'
	CHECK=.TRUE.
	DO WHILE(CHECK)							!第一次迭代,必须收敛!
					!	BR%BALNO(:)%INN%TY=.0_DP  !key
					!	BR%BALNO(:)%OUT%TY=.0_DP  !key
		CALL BROYDN(X,CHECK)		
	ENDDO
	F=FUNCV(X)
	CALL SAVEXf(X,F)
	CALL PRINTXF(X,F)
	MU0=BR%MU										!保存MU初值
!	BR%MU=BR%MU*.99_DP
	DO WHILE(BR%MU > 1.0E-5_DP)				!控制MU最小值
		WRITE(*,*) '*********   MU   ******'
		WRITE(*,*) '        ',BR%MU
		WRITE(*,*) '***********************'
		CHECK=.TRUE.
		K=0
		DO WHILE(CHECK .AND. K<10 .AND. NOT(MAXIT))
			K=K+1
			CALL BROYDN(X,CHECK)		
		ENDDO
		IF(K>=10 .OR. MAXIT) THEN
			CALL X0VCT(X)				!迭代初值
			BR%MU=(MU0+BR%MU)*0.5
			WRITE(*,*) 'BISECTION!'
			IF(ABS(MU0-BR%MU)<SQRT(EPSILON(MU0))) THEN
				MU0=MU0+2*SQRT(EPSILON(MU0))
				BR%MU=MU0+SQRT(EPSILON(MU0))
				WRITE(*,*) 'PERTURBATION'
			ENDIF
			MAXIT=.FALSE.		!置否，否则循环无法继续
			CYCLE   !结束本次循环
		ELSE
	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************  BY BRODYN METHOD  ****************'
	WRITE(*,*)'            ***************************************************'
			F=FUNCV(X)
			CALL SAVEXf(X,F)				!保存收敛解
			CALL PRINTXF(X,F)
!			CALL SLV07_HYBRD()
			MU0=BR%MU
			BR%MU=BR%MU*.933_DP
		ENDIF
		WRITE(7,"(E23.15,A1,E23.15,A1,E23.15,A1,E23.15,A1,E23.15,A1,E23.15)")  &
			X(1),tab,X(2),tab,X(3),tab,X(4),tab,X(5),tab,X(6)
	ENDDO
ENDSUBROUTINE SLV01_BROYDN

FUNCTION FUNCV(X)
	USE FSBR
	USE NRTYPE
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: X
	REAL(DP), DIMENSION(SIZE(X)) :: FUNCV

	CALL SLV6EQS(X,FUNCV,.false.)
ENDFUNCTION FUNCV

