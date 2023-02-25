	SUBROUTINE SLV08_NL2SOL()
		USE FSBR
		USE NRTYPE; USE NRUTIL
		USE NR
		IMPLICIT NONE

		INTEGER IV(BR%N+60), UIPARM(1)
		REAL(DP) ::V(93+BR%N*BR%M+3*BR%N+BR%M*(3*BR%M+33)/2)
		REAL(DP) ::X(BR%N), FVEC(BR%N),URPARM(1)
		EXTERNAL MADR, MADJ
  
	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************     BY NL2SOL ALGORITHM    ****************'
	WRITE(*,*)'         ***********************************************************'

		CALL DFAULT(IV,V)
		IV(19) = 0
		IV(24) = 1

		V(31)=1.0E-38_DP
		V(32)=1.0E-15_DP
		V(33)=1.0E-15_DP
		V(34)=1.0E-16_DP
		V(36)=1.0E-11_DP
		CALL X0VCT(X)
		CALL PRINTX0(X)
		CALL NL2SNO(BR%M, BR%N, X, MADR, IV, V, UIPARM, URPARM, MADR)

		CALL SLV6EQS(X,FVEC,.false.)
		CALL SAVEXf(X,FVEC)
		CALL PRINTXF(X,FVEC)
	ENDSUBROUTINE SLV08_NL2SOL

  SUBROUTINE MADR(N, P, X, NF, R, UIPARM, URPARM, UFPARM)
		USE NRTYPE
		IMPLICIT NONE
    INTEGER N, P, NF, UIPARM(1)
    REAL(DP) :: X(P), R(N), URPARM(1)
    EXTERNAL UFPARM
		CALL SLV6EQS(X,R,.false.)
  END

  SUBROUTINE MADJ(N, P, X, NF, J, UIPARM, URPARM, UFPARM)
		USE NRTYPE
		IMPLICIT NONE
    INTEGER N, P, NF, UIPARM(1)
    REAL(DP) :: X(P), J(N,P), URPARM(1)
    EXTERNAL UFPARM
    J(1,1) = 2.0_DP
    J(1,2) = 1.0_DP
    J(2,1) = 2._DP
    J(2,2) = -1.0_DP
!    J(3,1) = 0.0
!    J(3,2) = -SIN(X(2))
  END
