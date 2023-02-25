SUBROUTINE SLV09_SMSNO()
	USE FSBR
	IMPLICIT NONE
  EXTERNAL DEFLT, QDRTF, SMSNO
  INTEGER IMDCON
  INTEGER :: I, IV(60), LIV, LV, UIP(1)
  DOUBLE PRECISION D(BR%N), URP(BR%N,BR%N), X(BR%N)
	REAL(DP),ALLOCATABLE	::V(:)
	REAL(DP)	:: FV(BR%N)

	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************      BY SMSNO ALGORITHM    ****************'
	WRITE(*,*)'         ***********************************************************'

	LV=(77+BR%N*(BR%N+17)/2)
	LIV=60
	ALLOCATE(V(LV))

	CALL X0VCT(X)
	CALL PRINTX0(X)

	CALL DEFLT(2, IV, LIV, LV, V)
!	iv(1) = 0
	V(31)=1.0E-16_DP
	V(32)=1.0E-16_DP
	V(33)=1.0E-16_DP
	V(34)=1.0E-16_DP
	V(36)=1.0E-16_DP

  DO I = 1, BR%N
     D(I) = 1.0_DP/ABS(X(I))
     URP(I,1) = D(I)
  ENDDO

  CALL SMSNO(BR%N, D, X, QDRTF, IV, LIV, LV, V, UIP, URP, QDRTF)

	CALL SLV6EQS(X,FV,.false.)
	CALL SAVEXf(X,FV)
	CALL PRINTXF(X,FV)

ENDSUBROUTINE SLV09_SMSNO

SUBROUTINE QDRTF(N, X, NF, F, UIP, URP, UFP)
	USE FSBR
	IMPLICIT NONE
  INTEGER N, NF, UIP(1),I
  DOUBLE PRECISION X(N), F, URP(N,3)
  EXTERNAL UFP
	DOUBLE PRECISION R(BR%N)

	CALL SLV6EQS(X,R,.false.)
	DO I=1,N
		R(I)=R(I)**2
	ENDDO
	F=SUM(R)
ENDSUBROUTINE
