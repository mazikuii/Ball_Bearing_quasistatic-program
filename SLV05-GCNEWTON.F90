SUBROUTINE SLV05_GCNEWTON()
	USE FSBR
	USE NRTYPE
	USE NR
	USE INTFACE,ONLY:FUNCV
	IMPLICIT NONE
	INTEGER(I4B) :: k
	REAL(DP), DIMENSION(BR%M) :: F
	REAL(DP), DIMENSION(BR%N) :: X
	LOGICAL(LGT) :: CHECK
	CHARACTER(LEN=1)			:: TAB=CHAR(9)

	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************  BY NEWTON METHOD  ****************'
	WRITE(*,*)'            ***************************************************'

	CALL X0VCT(X)
	CALL PRINTX0(X)

	K=0
	CHECK=.TRUE.
	DO WHILE(CHECK .AND. K<13)
		CALL NEWT(X,CHECK)
		F=FUNCV(X)
		IF(CHECK) WRITE(*,*) 'CONVERGENCE PROBLEMS: NOT CONVERGED?!'
		CALL PRINTXF(X,F)
		K=K+1
	ENDDO

	IF(CHECK) THEN
		WRITE(*,*) 'a local minimum of the function fmin reached!'
	else
		WRITE(*,*) 'A ROOT OF FUNCV REACHED!'
	ENDIF
	CALL SAVEXf(X,f)
!	CALL WRTTM()  !输出拖动力、拖动力矩
ENDSUBROUTINE SLV05_GCNEWTON

!FUNCTION FUNCV(X)
!	SAME AS THAT USED IN BRODYN METHOD
!ENDFUNCTION FUNCV

