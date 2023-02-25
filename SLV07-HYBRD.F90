	SUBROUTINE SLV07_HYBRD()
	USE FSBR
	USE nrtype; USE nrutil
	USE nr
	USE INTFACE,ONLY:FUN
	IMPLICIT NONE
	REAL(DP), PARAMETER :: TOL=1.0e-16_DP
	INTEGER(I4B), PARAMETER :: N=6
	INTEGER(I4B) :: INFO,LWA
	REAL(DP), DIMENSION(N) :: fvec,x
	REAL(DP),ALLOCATABLE  :: WA(:)

	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************        BY HYBRD METHOD     ****************'
	WRITE(*,*)'         ***********************************************************'

	LWA=(N*(3*N+13))/2+1
	ALLOCATE(WA(LWA))
	CALL X0VCT(X)
	CALL PRINTX0(X)
	CALL HYBRD1(FUN,N,X,FVEC,TOL,INFO,WA,LWA)

	CALL SAVEXf(X,FVEC)
	CALL PRINTXF(X,FVEC)
!	CALL WRTTM()  !输出拖动力、拖动力矩
	ENDSUBROUTINE SLV07_HYBRD

	SUBROUTINE FUN(N,X,FVEC,IFLAG)
	USE NRTYPE; USE NRUTIL
	USE NR
	IMPLICIT NONE
	INTEGER			::N,IFLAG
	REAL(DP), DIMENSION(N) :: X,FVEC

	CALL SLV6EQS(X,FVEC,.false.)
	END SUBROUTINE FUN
