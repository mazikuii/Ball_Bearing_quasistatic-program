      SUBROUTINE SLV22_NLEQ2
      IMPLICIT DOUBLEPRECISION(S)
C
C     ____________________________________________________________
C
C     Testexample for NLEQ2: Computation of the Tchebychev
C     Polynom's evaluation points for quadrature-formulas.
C
C*  Written by        L. Weimann 
C*  Purpose           Testexample for code NLEQ2
C*  Version           2.3
C*  Revision          January 1992
C*  Latest Change     June 1992
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C
C     ____________________________________________________________
C
      INTEGER IRW
      PARAMETER (IRW=400)
      INTEGER IIW
      PARAMETER (IIW=61)
      INTEGER N
      PARAMETER (N=6)
      INTEGER I,N1
      DOUBLE PRECISION EPS
      INTEGER IOPT(50)
      INTEGER IERR,IFAIL
      DOUBLE PRECISION X(N),XSCAL(N),RW(IRW)
      INTEGER IW(IIW)
      REAL STIME,ETIME,CPTIME
      EXTERNAL F21
      EXTERNAL DF
C:    Begin
        EPS = 1.0D-16
        IOPT(:)=0
		IW(:)=0
		RW(:)=0.0D0
C       Execution mode: 0=Standard Mode, 1=Stepwise mode
        IOPT(2)=1
C       Jacobian: 0=(same as value 3)
C                 1=supplied by user routine JAC
C                 2=computed by numerical differentation (no feedback) 
C                 3=computed by numerical differentation (with feedback)
        IOPT(3)=1
C       Broyden updates: 0 = inhibit, 1=allow
C       IOPT(32)=1
C       Problem classification:
C       1 = linear , 2 = mildly nonlinear  3 = highly nonlinear
C       4 = extremely nonlinear
C       IOPT(31)=3
C       Set MPRERR, MPRMON, MPRSOL, MPRTIM
        IOPT(11)=3
        IOPT(13)=3
        IOPT(15)=2
        IOPT(19)=1
C       Set print units LUERR, LUMON, LUSOL, LUTIM
        IOPT(12)=6
        IOPT(14)=6
        IOPT(16)=6
        IOPT(20)=6
C       Solution output format:
C       0=standard format, 1= GRAZIL readable output
        IOPT(46)=0
C       Override maximum allowed number of iterations:
        IW(31)=200
C       Override initial pseudo-rank:
C       IW(32)=N
C       Override starting damping factor:
C       RW(21)=1.0D0
C       Override minimal allowed damping factor:
       RW(22)=1.0D-2
C       Override rank1-decision parameter SIGMA:
C       RW(23)=2.0D0
C       Override maximum permitted subcondition for DECCON:
C       RW(25)= 1.0D+16
		CALL X0VCT(X(1:N))
		CALL PRINTX0(X)					!µü´ú³õÖµ

        XSCAL(:) = 0.0D0
        IERR=-1
        I=0
31      IF (IERR.EQ.-1) THEN
          CALL NLEQ2(N,F21,DF,X,XSCAL,EPS,IOPT,IERR,IIW,IW,IRW,RW)
C         Clear workspace declared not to be used
          NIFREE=IW(16)
          IW(NIFREE:IIW)=0
          NRFREE=IW(17)
          RW(NRFREE:IRW)=0.0D0
          I=I+1
32         FORMAT('  Returned from call ',I4,' of NLEQ2')
          WRITE(6,32)I
C         IOPT(2)=0
          GOTO 31
        ENDIF

	CALL SLV6EQS(X,XSCAL,.false.)
	CALL SAVEXf(X,XSCAL)
	CALL PRINTXF(X,XSCAL)
      ENDSUBROUTINE SLV22_NLEQ2

