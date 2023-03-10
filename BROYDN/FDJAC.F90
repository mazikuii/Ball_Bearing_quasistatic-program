	SUBROUTINE FDJAC(X,FVEC,DF)
	USE FSBR
	USE NRTYPE; USE NRUTIL, ONLY : ASSERT_EQ
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: X
	REAL(DP), DIMENSION(:), INTENT(IN) :: FVEC
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: DF
	INTERFACE
		FUNCTION FUNCV(X)
		USE FSBR
		USE NRTYPE
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: X
		REAL(DP), DIMENSION(SIZE(X)) :: FUNCV
		END FUNCTION FUNCV
	END INTERFACE
!	REAL(DP), PARAMETER :: EPS=1.0E-4_DP
	REAL(DP) :: EPS
	INTEGER(I4B) :: J,N
	REAL(DP), DIMENSION(SIZE(X)) :: XSAV,XPH,H
	N=ASSERT_EQ(SIZE(X),SIZE(FVEC),SIZE(DF,1),SIZE(DF,2),'FDJAC')

!	EPS=1.0E-4_DP	!效果很差，很难收敛
	EPS=SQMPREC		!效果很好
!	EPS=1.0E-12_DP	!效果很好
!	EPS=EPSILON(X)

	XSAV=X
	H=EPS*ABS(XSAV)
	WHERE (H == 0.0_DP) H=EPS
	XPH=XSAV+H
	H=XPH-XSAV
	DO J=1,N
		X(J)=XPH(J)
		DF(:,J)=(FUNCV(X)-FVEC(:))/H(J)
		X(J)=XSAV(J)
	END DO
	END SUBROUTINE FDJAC

	SUBROUTINE FDJAC_D(X,FVEC,DF)
	USE FSBR
	USE NRTYPE; USE NRUTIL, ONLY : ASSERT_EQ
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: X
	REAL(DP), DIMENSION(:), INTENT(IN)		:: FVEC
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: DF
	INTERFACE
		FUNCTION FV(X)
		USE FSBR
		USE NRTYPE
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: X
		REAL(DP), DIMENSION(SIZE(X)) :: FV
		END FUNCTION FV
	END INTERFACE
	REAL(DP), PARAMETER :: EPS=1.0E-14_DP
	INTEGER(I4B) :: J,N
	REAL(DP), DIMENSION(SIZE(X)) :: XSAV,XPH,H
	N=ASSERT_EQ(SIZE(X),SIZE(FVEC),SIZE(DF,1),SIZE(DF,2),'FDJAC_D')
	XSAV=X
	H=EPS*ABS(XSAV)
	WHERE (H == 0.0_DP) H=EPS
	XPH=XSAV+H
	H=XPH-XSAV
	DO J=1,N
		X(J)=XPH(J)
		DF(:,J)=(FV(X)-FVEC(:))/H(J)
		X(J)=XSAV(J)
	END DO
	END SUBROUTINE FDJAC_D

