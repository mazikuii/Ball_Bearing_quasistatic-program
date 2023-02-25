	SUBROUTINE SLV13_DNEQBJ()
	USE FSBR
	USE nrtype; USE nrutil
	USE nr
	USE INTFACE,ONLY:FUN12,DFUN12
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=6,ITMAX=100
	REAL(DP),DIMENSION(N)		:: XGUESS,X,XSCALE,FSCALE,FVEC
	REAL(DP)								::	IPARAM(6), RPARAM(5)
	integer	i

	WRITE(*,*)'         ***********************************************************'
	WRITE(*,*)'         ***************     BY DNEQBJ SUBROUTINE   ****************'
	WRITE(*,*)'         ***********************************************************'

	CALL X0VCT(XGUESS)
	CALL PRINTX0(XGUESS)

	do i=1,n
		XSCALE(i)=1.0_DP !		abs(xguess(i))
		FSCALE(i)=1.0_DP !		abs(fvec(i))
	enddo

	CALL N4QBJ (IPARAM, RPARAM)
!	IPARAM(1)=0
	IPARAM(2)=14
	IPARAM(3)=400
	IPARAM(4)=400
	IPARAM(6)=1

!	RPARAM(1) =1.0E-12_DP
!	RPARAM(2) =1.0E-12_DP
	RPARAM(3) =1.0E-12_DP
	RPARAM(4) =1.0E-0_DP
	RPARAM(5) =1.0E-0_DP

	CALL DNEQBJ(FUN12,DFUN12,N,XGUESS,XSCALE,FSCALE,IPARAM, RPARAM,X,FVEC)

	CALL SAVEXf(X,FVEC)
	CALL PRINTXF(X,FVEC)
!	CALL WRTTM()  !输出拖动力、拖动力矩
	ENDSUBROUTINE SLV13_DNEQBJ

	SUBROUTINE DFUN12(N, X, FJAC,LDFJAC)		
		USE NRTYPE
		USE NR
		IMPLICIT NONE
		INTEGER			::N,LDFJAC
		REAL(DP) :: X(N),FJAC(N,N),FVEC(N)

		CALL SLV6EQS(X,FVEC,.false.)
		CALL FDJAC(X,FVEC,FJAC)
	END SUBROUTINE DFUN12
