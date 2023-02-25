SUBROUTINE SLV16_HOMPACK()
	USE FSBR
  USE REAL_PRECISION, ONLY : R8
  USE HOMPACK90, ONLY : FIXPDF, FIXPNF, FIXPQF
  IMPLICIT NONE
  INTEGER, PARAMETER:: N=6, NDIMA=6
  REAL (KIND=R8):: A(N),ANSAE,ANSRE,ARCAE,ARCRE,& 
  ARCLEN,DTIME,SSPAR(8),Y(N+1)
  INTEGER:: IFLAG,II,J,NFE,NP1,TIMENEW(8),TIMEOLD(8),TRACE
  CHARACTER (LEN=6) NAME
! DEFINE ARGUMENTS FOR CALL TO HOMPACK PROCEDURE.
!
  NP1=N+1
  ARCRE=0.5D-4
  ARCAE=0.5D-4
  ANSRE=1.0D-4
  ANSAE=1.0D-4
  TRACE=1
  SSPAR=0.0
  IFLAG=-2
	CALL X0VCT(Y(2:NP1))
	CALL PRINTX0(Y(2:NP1))					!迭代初值
	A=Y(2:NP1)
	EMBED=.false.		!嵌入式

! CALL TO HOMPACK ROUTINE.
  CALL FIXPNF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A, &
   SSPAR,NFE,ARCLEN)
!	
! CALCULATE EXECUTION TIME.
!
        WRITE (6,45) NAME
45      FORMAT (//,7X,'TESTING',1X,A6)
        WRITE (6,50) Y(1),IFLAG,NFE,ARCLEN,(Y(J),J=2,NP1)
50      FORMAT(/' LAMBDA =',F16.12,'  FLAG =',I2,I8,' JACOBIAN ', &
        'EVALUATIONS',/,1X,4X, &
        'ARCLEN =',F16.12/(1X,4ES16.8))
!
! TEST REVERSE CALL SUBROUTINES  STEPNX  AND  ROOTNX  ON THE SAME
! PROBLEM.
!
	CALL SLV6EQS(Y(2:NP1),A,.false.)
	CALL PRINTXF(Y(2:NP1),A)
ENDSUBROUTINE SLV16_HOMPACK
!
! SAMPLE USER WRITTEN HOMOTOPY SUBROUTINES FOR TESTING FIXP*F.
!
SUBROUTINE F(X,V)
!********************************************************************
!
!      SUBROUTINE F(X,V) -- EVALUATES BROWN'S FUNCTION AT THE POINT
!         X, AND RETURNS THE VALUE IN V.
!
!********************************************************************
  USE REAL_PRECISION, ONLY : R8
  IMPLICIT NONE
  REAL (KIND=R8), INTENT(IN):: X(:)
  REAL (KIND=R8), INTENT(OUT):: V(:)

	CALL SLV6EQS(X,V,.false.)
END SUBROUTINE F

SUBROUTINE FJAC(X,V,K)
!********************************************************************
!
!      SUBROUTINE FJAC(X,V,K)  --  EVALUATES THE K-TH COLUMN OF
!         THE JACOBIAN MATRIX FOR BROWN'S FUNCTION EVALUATED AT
!         THE POINT X, RETURNING THE VALUE IN V.
!
!********************************************************************
	USE nr, ONLY : fdjac
  USE REAL_PRECISION, ONLY : R8
  REAL (KIND=R8):: X(:)
  REAL (KIND=R8), INTENT(OUT):: V(:)
  INTEGER, INTENT(IN):: K

  REAL (KIND=R8) :: FVEC(6),JAC(6,6)

	CALL SLV6EQS(X,FVEC,.false.)
	CALL FDJAC(X,FVEC,JAC)

	V=JAC(:,K)
END SUBROUTINE FJAC
! **********************************************************************
!
! THE REST OF THESE SUBROUTINES ARE NOT USED BY PROGRAM TESTF, AND ARE
! INCLUDED HERE SIMPLY FOR COMPLETENESS AND AS TEMPLATES FOR THEIR USE.
!
SUBROUTINE RHO(A,LAMBDA,X,V)
	USE FSBR
  USE REAL_PRECISION, ONLY : R8
  USE HOMPACK90_GLOBAL
  REAL (KIND=R8), INTENT(IN):: A(:),X(:)
  REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
  REAL (KIND=R8), INTENT(OUT):: V(:)

	REAL (KIND=R8)	::FVEC(6)
!
! EVALUATE  RHO(A,LAMBDA,X)  AND RETURN IN THE VECTOR  V .
!
	IF(EMBED) THEN
		CALL CONVT(BR%MU,LAMBDA)
		CALL SLV6EQS(X,V,.false.)
	ELSE
		CALL SLV6EQS(X,FVEC,.false.)
		V=LAMBDA*FVEC+(1.0-LAMBDA)*(X-A)
	ENDIF
END SUBROUTINE RHO

SUBROUTINE RHOJAC(A,LAMBDA,X,V,K)
	USE nr, ONLY : fdjac
	USE FSBR
  USE REAL_PRECISION, ONLY : R8
  USE HOMPACK90_GLOBAL
  REAL (KIND=R8), INTENT(IN):: A(:)
  REAL (KIND=R8):: X(:)
  REAL (KIND=R8), INTENT(IN OUT):: LAMBDA
  REAL (KIND=R8), INTENT(OUT):: V(:)
  INTEGER, INTENT(IN):: K
	CHARACTER(len=1)			:: tab=CHAR(9)
  REAL (KIND=R8) :: FVEC(6),JAC(6,6)
! RETURN IN THE VECTOR  V  THE KTH COLUMN OF THE JACOBIAN
! MATRIX [D RHO/D LAMBDA, D RHO/DX] EVALUATED AT THE POINT
! (A, LAMBDA, X).

	IF(EMBED) THEN
		CALL CONVT(BR%MU,LAMBDA)
		CALL SLV6EQS(X,FVEC,.false.)
		IF(K==1) THEN
			CALL DFDL(X,LAMBDA,FVEC,V)		!D RHO/D LAMBDA
		ELSE
			CALL FDJAC(X,FVEC,JAC)
			V=(LAMBDA+0.0)*JAC(:,K-1)      !错误所在处，奶奶的
			V(K-1)=V(K-1)+(1-LAMBDA)
		ENDIF
	ELSE
		CALL SLV6EQS(X,FVEC,.false.)
		IF(K==1) THEN
			V=FVEC-(X-A)			!D RHO/D LAMBDA
		ELSE
			CALL FDJAC(X,FVEC,JAC)
			V=(LAMBDA+0.0)*JAC(:,K-1)      !错误所在处，奶奶的
			V(K-1)=V(K-1)+(1-LAMBDA)
		ENDIF
	ENDIF
END SUBROUTINE RHOJAC      

SUBROUTINE DFDL(X,LAMBDA,FVEC,DF)
!计算D F/D LAMBDA
	USE FSBR
	USE NRTYPE; USE NRUTIL, ONLY : ASSERT_EQ
	USE INTFACE,ONLY:FUNCV
	IMPLICIT NONE
	REAL(DP), DIMENSION(6), INTENT(IN) :: FVEC
	REAL(DP), DIMENSION(6), INTENT(IN) :: X
	REAL(DP), INTENT(IN)				:: LAMBDA
	REAL(DP), DIMENSION(6), INTENT(OUT) :: DF
	REAL(DP) :: EPS
	REAL(DP) ::H,TEMP
	REAL(DP)	:: FVECN(SIZE(X))

	EPS=SQMPREC

	H=EPS*ABS(LAMBDA)
	IF(H == 0.0_DP) H=EPS
	TEMP=LAMBDA+H
	H=TEMP-LAMBDA
	CALL CONVT(BR%MU,TEMP)
	
	CALL SLV6EQS(X,FVECN,.false.)
	DF=(FVECN-FVEC)/H
END SUBROUTINE DFDL

SUBROUTINE CONVT(MU,LAMBDA)
	USE NRTYPE
	IMPLICIT NONE
	REAL(DP),INTENT(OUT)	:: MU
	REAL(DP),INTENT(IN)		:: LAMBDA

	MU=0.6_DP*LAMBDA+(1.0_DP-LAMBDA)*0.1_DP
ENDSUBROUTINE CONVT

	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RHOA(A,LAMBDA,X)
  USE REAL_PRECISION, ONLY : R8
  REAL (KIND=R8), INTENT(OUT):: A(:)
  REAL (KIND=R8), INTENT(IN):: LAMBDA,X(:)

! CALCULATE AND RETURN IN  A  THE VECTOR Z SUCH THAT
!  RHO(Z,LAMBDA,X) = 0 .
!
	A=0.0
  RETURN
END SUBROUTINE RHOA


SUBROUTINE FJACS(X)
  USE REAL_PRECISION, ONLY : R8
  USE HOMPACK90_GLOBAL
  REAL (KIND=R8), INTENT(IN):: X(:)
!
! If MODE = 1,
! evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
! the result in packed skyline storage format in QRSPARSE.  LENQR is the
! length of QRSPARSE, and ROWPOS contains the indices of the diagonal
! elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
! ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPOS
! is not used by this storage format.
!
! If MODE = 2,
! evaluate the N x N Jacobian matrix of F(X) at X, and return the result
! in sparse row storage format in QRSPARSE.  LENQR is the length of
! QRSPARSE, ROWPOS contains the indices of where each row begins within
! QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
! the corresponding elements in QRSPARSE.  Even if zero, the diagonal
! elements of the Jacobian matrix must be stored in QRSPARSE.
!
  RETURN
END SUBROUTINE FJACS

SUBROUTINE RHOJS(A,LAMBDA,X)
  USE REAL_PRECISION, ONLY : R8
  USE HOMPACK90_GLOBAL
  REAL (KIND=R8), INTENT(IN):: A(:),LAMBDA,X(:)
!
! If MODE = 1,
! evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
! the result in packed skyline storage format in QRSPARSE.  LENQR is the
! length of QRSPARSE, and ROWPOS contains the indices of the diagonal
! elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
! ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPOS
! is not used by this storage format.
!
! If MODE = 2,
! evaluate the N x N Jacobian matrix of F(X) at X, and return the result
! in sparse row storage format in QRSPARSE.  LENQR is the length of
! QRSPARSE, ROWPOS contains the indices of where each row begins within
! QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
! the corresponding elements in QRSPARSE.  Even if zero, the diagonal
! elements of the Jacobian matrix must be stored in QRSPARSE.
!
  RETURN
END SUBROUTINE RHOJS
! **********************************************************************
!
!  SUBROUTINE TO TEST THE REVERSE CALL SUBROUTINES  STEPNX  AND
!  ROOTNX.  THE TEST PROBLEM IS BROWN'S FUNCTION, ZERO FINDING.
!  THE OUTPUT IS SIMILAR TO THAT FROM THE TEST OF  FIXPNF, EXCEPT WITH
!  MORE JACOBIAN EVALUATIONS SINCE THE UNDEFINED FUNCTION OPTION OF
!  STEPNX  IS USED TO FORCE SMALLER STEPS.
!
