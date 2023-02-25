SUBROUTINE SLV20_NITSOL()
	USE FSBR
	USE NRTYPE
	USE INTFACE, ONLY:F6EQ,J6EQ
	IMPLICIT NONE
  
	external ddot, dnrm2

  INTEGER     MAXKD
  PARAMETER ( MAXKD=50 )

  INTEGER     MAXNX,     MAXNY,     MAXN
  PARAMETER ( MAXNX=128, MAXNY=128, MAXN=MAXNX*MAXNY ) 

  INTEGER     LRPAR
  PARAMETER ( LRPAR=(MAXNX+2)*(MAXNY+2)+MAXNX*MAXNY+5 )

  INTEGER     LRWORK 
  PARAMETER ( LRWORK=MAXN*(MAXKD+10)+MAXKD*(MAXKD+3))

  INTEGER N
  DOUBLE PRECISION X(MAXN)
  DOUBLE PRECISION FTOL
  DOUBLE PRECISION STPTOL
  INTEGER INPUT(10)
  INTEGER INFO(6)
  DOUBLE PRECISION RWORK(LRWORK) 
  DOUBLE PRECISION RPAR(LRPAR)
  INTEGER IPAR(2)
  INTEGER ITERM


	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************  BY NITSOL METHOD  ****************'
	WRITE(*,*)'            ***************************************************'

	N=6
	FTOL		=1.0E-15_DP
  STPTOL	=1.0E-15_DP

	input(3)=0

	CALL X0VCT(X)
	CALL PRINTX0(X)					!µü´ú³õÖµ
  CALL NITSOL(N, X, F6EQ, J6EQ, FTOL, STPTOL, &
      INPUT, INFO, RWORK, RPAR, IPAR, ITERM, DDOT, DNRM2)

	WRITE(*,*) 'iterm:',iterm
	IF(ITERM <0) WRITE(*,*) 'illegal value in input(k)'
	IF(ITERM==0) WRITE(*,*) 'normal termination: ||F||.le.ftol or ||step||.le.stptol.'
	IF(ITERM==1) WRITE(*,*) 'nnimax nonlinear iterations reached without success.'
	IF(ITERM==2) WRITE(*,*) 'failure to evaluate F.'
	IF(ITERM==3) WRITE(*,*) 'in nitjv, J*v failure.'
	IF(ITERM==4) WRITE(*,*) 'in nitjv, P(inverse)*v failure.'
	IF(ITERM==5) WRITE(*,*) 'in nitdrv, insufficient initial model norm reduction for adequate progress.'
	IF(ITERM==6) WRITE(*,*) 'in nitbt, failure to reach an acceptable step through backtracking.'
	CALL PRINTXF(X,RWORK)
      write(*,920) info(1)
      write(*,930) info(2)
      write(*,940) info(3)
      write(*,950) info(4)
      write(*,960) info(5)
      write(*,970) info(6)

 920  format(' No. function evaluations:     ', i9)
 930  format(' No. J*v evaluations:          ', i9) 
 940  format(' No. P(inverse)*v evaluations: ', i9)
 950  format(' No. linear iterations:        ', i9)
 960  format(' No. nonlinear iterations:     ', i9)
 970  format(' No. backtracks:               ', i9)

ENDSUBROUTINE SLV20_NITSOL

SUBROUTINE F6EQ(N, XCUR, FCUR, RPAR, IPAR, ITRMF)
	USE NRTYPE
	IMPLICIT NONE
	INTEGER I, ITRMF, J, J1, J2, N, IPAR(*), NX, NY 
	REAL(DP)	:: CL, CR, H2L, XCUR(N), FCUR(N), RPAR(*)

	CALL SLV6EQS(XCUR,FCUR,.false.)
ENDSUBROUTINE F6EQ

SUBROUTINE J6EQ(N,XCUR,FCUR,IJOB,V,Z,RPAR,IPAR,ITRMJV)
	IMPLICIT NONE
  integer n, ijob, ipar(*), itrmjv
  double precision xcur(n), fcur(n), v(n), z(n), rpar(*)

ENDSUBROUTINE J6EQ
