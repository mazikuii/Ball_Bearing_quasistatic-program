SUBROUTINE SLV17_BFGS()
!	driver for routine dfpmin
	USE nrtype
	USE nr
	USE INTFACE,ONLY:func,dfunc
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=6
	REAL(DP), PARAMETER :: GTOL=1.0e-12_dp
	INTEGER(I4B) :: iter
	REAL(DP) :: fret
	REAL(DP), DIMENSION(NDIM) :: p,F

	CALL X0VCT(P)
	CALL PRINTX0(P)					!µü´ú³õÖµ
	call dfpmin(p,GTOL,iter,fret,func,dfunc)

	write(*,'(1x,a,i3)') 'Iterations:',iter
	write(*,'(1x,a,e14.6)') 'Func. value at solution',fret
	CALL SLV6EQS(P,F,.false.)
	CALL PRINTXF(P,F)
	F=DFUNC(P)
	WRITE(*,*) 'CONVERGENCE GRADIENCE'
	CALL PRINTXF(P,F)

ENDSUBROUTINE SLV17_BFGS

FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(dP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP) :: func
	REAL(DP) :: F(SIZE(X))
	
	CALL SLV6EQS(X,F,.false.)
	FUNC=dot_product(F,F)
END FUNCTION func

FUNCTION dfunc(x)
!¼ÆËãD F/D LAMBDA
	USE FSBR
	USE nr
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: dfunc

	REAL(DP), DIMENSION(size(x)) :: H,XI
	REAL(DP) :: EPS,FI,F0
	INTEGER	::	I
	INTERFACE
		FUNCTION func(p)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: p
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE

	EPS=SQMPREC

	F0=FUNC(X)
	H=EPS*ABS(X)
	WHERE(H == 0.0_DP) H=EPS
	DO I=1,SIZE(X)
		XI=X
		XI(I)=X(I)+H(I)
		FI=FUNC(XI)
		dfunc(I)=(FI-F0)/H(I)
	ENDDO
END FUNCTION dfunc
