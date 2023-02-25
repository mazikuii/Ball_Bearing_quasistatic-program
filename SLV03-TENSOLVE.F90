SUBROUTINE SLV03_TENSOLVE()
	USE FSBR
	USE NRTYPE
	USE NR
	USE INTFACE,ONLY:frosen,jwood
	IMPLICIT NONE
	INTEGER(I4B) :: I
  INTEGER(I4B)	::maxm, maxn, maxp, m, n, itnlim, jacflg, method,global, ipr, lunc, lnem, lnen, lin, msg, termcd
  REAL(DP)			::gradtl, steptl, ftol, stepmx, dlt
  parameter       (maxm = 100, maxn = 30, maxp = 6)
  parameter       (lin = 3, lunc = 14, lnem = 51, lnen = 19)
  INTEGER(I4B)	::iwrkn(maxn,lin)
  REAL(DP)			::x0(maxn), xp(maxn), fp(maxm), gp(maxn),typx(maxn), typf(maxm),wrknen(maxn,lnen), wrkunc(maxp,lunc),wrknem(maxm,lnem)

	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************  BY TENSOLVE METHOD ***************'
	WRITE(*,*)'            ***************************************************'

! Set dimensions of the problem.
  m= 6
  n= 6
! Set default values for the TENSOLVE parameters.
  call tsdflt(m     , n     , itnlim, jacflg, gradtl, steptl,		&
              ftol  , method, global, stepmx, dlt   ,						&
              typx  , typf  , ipr   , msg    )
	if(m==6 .and. n==6) then
	! Set values for the initial point.      
		CALL X0VCT(X0)
		CALL PRINTX0(X0)
	! Alter some of the parameters.
		gradtl	= 1.0E-16_DP
		ftol		= 1.0E-16_DP
		steptl	= 1.0E-16_DP
		global	= 0   !TRUST REGION METHOD
		method	= 1   !1张量法更快速、准确。
		jacflg	=	0
!		MSG			= 2

	elseif(m==2 .and. n==2) then
!     Set values for the initial point.      
      x0(1)  = 01.0000_DP
      x0(2)  = 7.0_DP
	endif
! Call TENSOLVE.
  call tsneci(maxm  , maxn  , maxp  , x0    , m     , n     ,&
             typx  , typf  , itnlim, jacflg, gradtl, steptl, &
             ftol  , method, global, stepmx, dlt   , ipr   , &
             wrkunc, lunc  , wrknem, lnem  , wrknen, lnen  ,&
             iwrkn , lin   , frosen , jwood , msg   , &
             xp    , fp    , gp    , termcd )
	CALL SAVEXf(XP,FP)
	CALL PRINTXF(XP,FP)
ENDSUBROUTINE SLV03_TENSOLVE

SUBROUTINE FROSEN ( X, F, M, N )
	USE FSBR
	USE NRTYPE
	IMPLICIT NONE
	INTEGER   ::M, N
	REAL(DP)	::X(N), F(M)

	IF(SIZE(X)==6) THEN
		CALL SLV6EQS(X,F,.false.)
	ELSEIF(SIZE(X)==2) THEN
		F(1) = 1.0D4 * X(2)* X(1)-1.0D0
		F(2) = EXP(-X(1))+EXP(-X(2))-1.001D0
	ENDIF
END

SUBROUTINE JWOOD ( X, JAC, MAXM, M, N )
	USE nr, ONLY : fdjac
	implicit none
  INTEGER            MAXM, M, N
  DOUBLE PRECISION   X(N), JAC(MAXM,N)
	double precision   fvec(n),fjac(m,n)

	CALL SLV6EQS(X,FVEC,.false.)
	CALL FDJAC(X,FVEC,FJAC)
	
	jac(1:m,:)=fjac
END

