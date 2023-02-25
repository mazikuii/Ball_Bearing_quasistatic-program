      SUBROUTINE TT2
C
C     ------------------------------------------------------------
C
C     Testexample for ALCON2
C
C*  Purpose           Testexample for code ALCON2
C*  Version           0.9
C*  Revision          September 1985
C*  Latest Change     January 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C
C     ------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(10),XW(10),RWORK(5000),IWORK(1000),INFO(9)
      EXTERNAL FTEST,JTEST
      INTEGER UPR,UDIAG,UBIF
      COMMON /UNIT/ UPR,UDIAG
      UPR=6
      UDIAG=2
      OPEN(2,FILE='alcon2.diag')
      LRW=5000
      LIW=1000
      N=2
      X(1)=0.D0
      X(2)=0.D0
      XW(1)=1.D0
      XW(2)=1.D0
      TAU=0.D0
      TAUMIN=0.D0
      TAUMAX=1.D0
      UMAX=0.D0
      EPS=1.D-5
      INFO(1)=2
      INFO(2)=99
      INFO(3)=0
      INFO(4)=1
      CALL ALCON2(FTT2,JTT2,N,X,XW,TAU,TAUMIN,TAUMAX,UMAX,EPS,INFO,
     &   RWORK,LRW,IWORK,LIW)
      STOP
      END

      SUBROUTINE FTT2(X,TAU,F)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION F(2),X(2)
      DATA P1/1.D0/,P2/0.D0/,R/2.D0/
      F(1)=(X(1)**2-3.D0*R*X(1)+2.D0*R**2)*X(1)+
     &   (X(1)*X(2)-R*X(2))*X(2)-TAU*P1
      F(2)=(X(1)*X(2)-R*X(2))*X(1)+(X(2)**2-R*X(1)+2.D0)*X(2)-TAU*P2
      RETURN
      END

      SUBROUTINE JTT2(X,TAU,DFDX,DFDTAU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(2),DFDX(2,2),DFDTAU(2)
      DATA P1/1.D0/,P2/0.D0/,R/2.D0/
      DFDX(1,1)=3.D0*X(1)**2-6.D0*R*X(1)+2.D0*R**2+X(2)**2
      DFDX(1,2)=2.D0*X(1)*X(2)-2.D0*R*X(2)
      DFDX(2,1)=2.D0*X(1)*X(2)-2.D0*R*X(2)
      DFDX(2,2)=X(1)**2-2.D0*R*X(1)+3.D0*X(2)**2+2.D0
      DFDTAU(1)=-P1
      DFDTAU(2)=-P2
      RETURN
      END
