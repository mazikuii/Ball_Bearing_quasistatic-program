!FASTSIM???Գ???
SUBROUTINE FTEST(NO,FS,BR,VX,VY,PHI)
	USE GLOBAL
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: NO
	TYPE(FASTDATA),TARGET	:: FS
	TYPE(BEARSYSDATA),INTENT(INOUT),TARGET	:: BR
	REAL(DP) ::VX,VY,PHI,G1,LOAD,MU,A,B,C11,C22,C23
	REAL(DP) ::L1,L2,L3,NP,Z0,UX,UY,PHX,PHY
	REAL(DP)	::TX,TY
	LOGICAL	:: INOUT

G1=8.4E11_DP		!???Գ???
LOAD=1.0E3_DP
MU=0.3_DP
IF(NO==1) THEN
	A=6.0E-3_DP		!???Գ???
	B=6.0E-3_DP		!???Գ???
	C11=4.12_DP		!???Գ???
	C22=3.67_DP		!???Գ???
	C23=1.47_DP		!???Գ???
ELSEIF(NO==2) THEN
	A=7.5E-3_DP		!???Գ???
	B=1.5E-3_DP		!???Գ???
	C11=7.78_DP				!???Գ???
	C22=8.14_DP				!???Գ???
	C23=6.63_DP				!???Գ???
ELSE
	A=1.5E-3_DP		!???Գ???
	B=7.5E-3_DP		!???Գ???
	C11=3.37_DP				!???Գ???
	C22=2.63_DP				!???Գ???
	C23=0.603_DP				!???Գ???
ENDIF

L1=8.0_DP*A/3.0_DP/C11/G1
L2=8.0_DP*A/3.0_DP/C22/G1
L3=PI_D*A**2/4.0_DP/G1/SQRT(A*B)/C23
NP=2._DP*PI_D/3._DP     !????Hertz?Ӵ????м???
Z0=LOAD/A/B/NP

UX =A*VX/MU/Z0/L1
UY =A*VY/MU/Z0/L2
PHX=A*B*PHI/MU/Z0/L3 
PHY=A**2*PHI/MU/Z0/L3
!!!!!!!!!!!!!!!!!
!UX =1.0_dp
!UY =-2.0_dp
!PHX=2.0_dp
!PHY=4.0_dp
!!!!!!!!!!!!!!!!!
!FS% TRACTION	=.TRUE.			!???????϶???   ????ȷ!
FS% HEATHCOTE	=.FALSE.			!HEATHCOTE SLIP   FALSE
FS% HERTZ			=.TRUE.				!Hertz Contact Model

FS%UX =UX
FS%UH =0.0_dp
FS%UY =UY
FS%PHX=PHX
FS%PHY=PHY

CALL FASTSIM(FS)

INOUT=.FALSE.
CALL POSTFS(FS,BR,1,INOUT)

IF(INOUT) THEN
	TX=BR%BALNO(1)%OUT%TX/BR%BALNO(1)%OUT%P/BR%MU
	TY=BR%BALNO(1)%OUT%TY/BR%BALNO(1)%OUT%P/BR%MU
ELSE
	TX=BR%BALNO(1)%INN%TX/BR%BALNO(1)%INN%P/BR%MU
	TY=BR%BALNO(1)%INN%TY/BR%BALNO(1)%INN%P/BR%MU
ENDIF

WRITE(*,*) TX,TY
TX=TX*LOAD*MU
TY=TY*LOAD*MU
WRITE(*,*) TX,TY

ENDSUBROUTINE FTEST

!FASTSIM???Գ???
SUBROUTINE TEST(FS,BR,VX,VY,PHI)
	USE GLOBAL
	IMPLICIT NONE
	TYPE(FASTDATA),TARGET	:: FS
	TYPE(BEARSYSDATA),INTENT(INOUT),TARGET	:: BR
	REAL(DP) ::VX,VY,PHI,G1,LOAD,MU,A,B,C11,C22,C23
	REAL(DP) ::L1,L2,L3,NP,Z0,UX,UY,PHX,PHY
	REAL(DP)	::TX,TY

G1=8.4E10_DP		!???Գ???
LOAD=1.0E5_DP
MU=0.3_DP

A=6.0E-3_DP				!???Գ???
B=6.0E-3_DP				!???Գ???
C11=4.12_DP				!???Գ???
C22=3.67_DP				!???Գ???
C23=1.47_DP				!???Գ???

L1=8.0_DP*A/3.0_DP/C11/G1
L2=8.0_DP*A/3.0_DP/C22/G1
L3=PI_D*A**2/4.0_DP/G1/SQRT(A*B)/C23
NP=2._DP*PI_D/3._DP     !????Hertz?Ӵ????м???
Z0=LOAD/A/B/NP

UX =A*VX/MU/Z0/L1
UY =A*VY/MU/Z0/L2
PHX=A*B*PHI/MU/Z0/L3 
PHY=A**2*PHI/MU/Z0/L3

FS% HEATHCOTE	=.FALSE.			!HEATHCOTE SLIP   FALSE
FS% HERTZ			=.TRUE.				!Hertz Contact Model

FS%UX =UX
FS%UH =0.0_dp
FS%UY =UY
FS%PHX=PHX
FS%PHY=PHY

BR%MU=MU
BR%BALNO(:)%INN%P=LOAD
BR%BALNO(:)%OUT%P=LOAD
ENDSUBROUTINE TEST

SUBROUTINE FTEST2(FS,BR,VX,VY,PHI)
	USE GLOBAL
	IMPLICIT NONE
	TYPE(FASTDATA),TARGET	:: FS
	TYPE(BEARSYSDATA),INTENT(INOUT),TARGET	:: BR
	REAL(DP) ::VX,VY,PHI

	FS% HEATHCOTE	=.FALSE.		!HEATHCOTE SLIP   FALSE
	FS% HERTZ		=.TRUE.			!Hertz Contact Model
	BR% HERTZ		=FS% HERTZ		!Hertz Contact Model
	BR% SIMPMOMENT	=.TRUE.			!???ü򻯷????????????????ĵ?????
	BR% SYM			=.TRUE.			!?Ƿ??Գ??غ?

BR%G1=8.4E10_DP		
BR%G2=8.4E10_DP		
BR%BALNO(:)%OUT%P=1.0E5_DP
BR%BALNO(:)%INN%P=1.0E5_DP
BR%MU=0.3_DP
BR%BALNO(:)%OUT%AX=6.0E-3_DP
BR%BALNO(:)%OUT%BY=6.0E-3_DP		
BR%BALNO(:)%OUT%AB=BR%BALNO(1)%OUT%AX*BR%BALNO(1)%OUT%BY		
BR%BALNO(:)%INN%AX=6.0E-3_DP
BR%BALNO(:)%INN%BY=6.0E-3_DP
BR%BALNO(:)%INN%AB=BR%BALNO(1)%INN%AX*BR%BALNO(1)%INN%BY
BR%OCOM%KE=1.0_DP		
BR%ICOM%KE=1.0_DP		
BR%OCOM%C11=4.12_DP				
BR%OCOM%C22=3.67_DP
BR%OCOM%C23=1.47_DP
BR%ICOM%C11=4.12_DP
BR%ICOM%C22=3.67_DP
BR%ICOM%C23=1.47_DP

BR%BALNO(:)%OUT%XIX=VX		
BR%BALNO(:)%OUT%XIY=VY		
BR%BALNO(:)%OUT%PSI=PHI		
BR%BALNO(:)%INN%XIX=VX		
BR%BALNO(:)%INN%XIY=VY		
BR%BALNO(:)%INN%PSI=PHI		

CALL GETCREEP(BR)

BR%BALNO(:)%OUT%UH=0.0_DP
BR%BALNO(:)%INN%UH=0.0_DP

FS%UX =BR%BALNO(1)%OUT%UX
FS%UH =BR%BALNO(1)%OUT%UH
FS%UY =BR%BALNO(1)%OUT%UY
FS%PHX=BR%BALNO(1)%OUT%PHX
FS%PHY=BR%BALNO(1)%OUT%PHY
CALL FASTSIM(FS)
CALL POSTFS(FS,BR,1,.TRUE.)

FS%UX =BR%BALNO(1)%INN%UX
FS%UH =BR%BALNO(1)%INN%UH
FS%UY =BR%BALNO(1)%INN%UY
FS%PHX=BR%BALNO(1)%INN%PHX
FS%PHY=BR%BALNO(1)%INN%PHY
CALL FASTSIM(FS)
CALL POSTFS(FS,BR,1,.FALSE.)
WRITE(*,*) BR%BALNO(1)%OUT%TX,BR%BALNO(1)%OUT%TY
WRITE(*,*) BR%BALNO(1)%INN%TX,BR%BALNO(1)%INN%TY
ENDSUBROUTINE FTEST2

!????Ħ??ϵ????
SUBROUTINE FCOE(FS,BR)
	USE GLOBAL
	USE INTFACE,ONLY:POW
	IMPLICIT NONE
	TYPE(FASTDATA),TARGET	:: FS
	TYPE(BEARSYSDATA),INTENT(INOUT),TARGET	:: BR
	CHARACTER(len=1)			:: tab=CHAR(9)

	INTEGER	:: I,N

	FS% HEATHCOTE	=.FALSE.		!HEATHCOTE SLIP   FALSE
	FS% HERTZ		=.TRUE.			!Hertz Contact Model
	BR% HERTZ		=FS% HERTZ		!Hertz Contact Model
	BR% SIMPMOMENT	=.TRUE.			!???ü򻯷????????????????ĵ?????
	BR% SYM			=.TRUE.			!?Ƿ??Գ??غ?

	OPEN(5, file="t.XLS")  !?????????ļ? EXCEL?ļ?
	
	N=100
	DO I=1,N
		BR%MU=1.0E-10_dp*pow(10.0_dp,1.0_dp*I/10)
!		BR%MU=0.3_dp
		CALL GETCREEP(BR)
		CALL PREFS(FS,BR,1,NOT(.TRUE.))
		CALL FASTSIM(FS)
		CALL POSTFS(FS,BR,1,NOT(.TRUE.))

!		CALL VECTOR_COMBINING(FS)  !ʸ???ϳ?x??y??????????????ǰһ????
!		CALL PRINT_FIELD(FS)
	!	CALL PRINT_VFIELD(FS)
!		WRITE(*,*) BR%BALNO(1)%INN%TX,TAB,BR%BALNO(1)%INN%TY
		WRITE(5,"(e23.15,A1,e23.15,A1,e23.15)") BR%MU,TAB,BR%BALNO(1)%INN%TX,TAB,BR%BALNO(1)%INN%TY
	ENDDO
	CLOSE(5)
ENDSUBROUTINE FCOE

!???????????ٶȶԺ????϶?????Ӱ??
SUBROUTINE SPINY(FS,BR)
	USE GLOBAL
	USE INTFACE,ONLY:POW
	IMPLICIT NONE
	TYPE(FASTDATA),TARGET	:: FS
	TYPE(BEARSYSDATA),INTENT(INOUT),TARGET	:: BR
	CHARACTER(len=1)			:: tab=CHAR(9)

	INTEGER	:: I,N

	FS% HEATHCOTE	=.FALSE.		!HEATHCOTE SLIP   FALSE
	FS% HERTZ		=.TRUE.			!Hertz Contact Model
	BR% HERTZ		=FS% HERTZ		!Hertz Contact Model
	BR% SIMPMOMENT	=.TRUE.			!???ü򻯷????????????????ĵ?????
	BR% SYM			=.TRUE.			!?Ƿ??Գ??غ?

!	OPEN(5, file="spiny.XLS")  !?????????ļ? EXCEL?ļ?

	BR%BALNO(:)%OUT%P =1.0E1_DP 
	BR%MU =0.00000001 
	CALL PCHAXIS(BR)
	
	N=1
	DO I=1,N
		BR%BALNO(:)%OUT%XIX=0.0_DP
		BR%BALNO(:)%OUT%XIY=0.0_DP
		BR%BALNO(:)%OUT%PSI=-.0025_DP
		CALL GETCREEP(BR)		!????FASTSIM?????ĸ??ֲ???
		CALL PREFS(FS,BR,1,.TRUE.)  
		CALL FASTSIM(FS)
		CALL POSTFS(FS,BR,1,.TRUE.)

		CALL PRINT_SFIELD(FS)
		CALL PRINT_TFIELD(FS)
!		WRITE(5,"(e23.15,A1,e23.15)") FS%TFX/br%np,TAB,-FS%TFY	!/br%np
		WRITE(*,*) FS%TFX,TAB,-FS%TFY	
	ENDDO
!	CLOSE(5)
ENDSUBROUTINE SPINY

!??????Բ?ʡ???Բ????
SUBROUTINE EITEST(BR)
	USE GLOBAL
	USE INTFACE,ONLY:POW
	USE nr, ONLY : FDJAC_D
	IMPLICIT NONE
	TYPE(BEARSYSDATA),INTENT(INOUT),TARGET	:: BR
	INTERFACE
		FUNCTION FV(X)
		USE NRTYPE
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: X
		REAL(DP), DIMENSION(SIZE(X)) :: FV
		END FUNCTION FV
	END INTERFACE

	TYPE(PATCHABK)	 ,TARGET	:: ABK
	CHARACTER(LEN=1)			:: TAB=CHAR(9)
	INTEGER	:: I
	!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER(I4B), PARAMETER :: NMAX=1024,MMAX=17
  CHARACTER*60     TASK, CSAVE
  LOGICAL      LSAVE(4)
  INTEGER	::  M, IPRINT,NBD(NMAX), IWA(3*NMAX), ISAVE(44)
  REAL(DP)::	F, FACTR, PGTOL,X(NMAX), L(NMAX), U(NMAX), G(NMAX), DSAVE(29) 
  REAL(DP)::   WA(2*MMAX*NMAX+4*NMAX+12*MMAX*MMAX+12*MMAX)
	!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER :: N    !???̸???
  REAL(DP)::   t1, t2
	REAL(DP), DIMENSION(2) :: XV,FVC
	REAL(DP), DIMENSION(2,2) :: DF

	OPEN(7, FILE="EITEST.XLS")  !?????????ļ? EXCEL?ļ?
	task = 'START'      !On first entry, it must be set to 'START'
	iprint = 1		!output control
  factr=1.0E-1_DP	!tolerances in the stopping criteria.
  pgtol=0.0d-16	!tolerances in the stopping criteria.

	N=2
  m=20

	x(1)=4.5_DP
	x(2)=4.5_DP
	DO I=1,N
		nbd(I)=2.0_DP		!2 if x(i) has both lower and upper bounds
		l(I)=3.0_DP			!lower bounds
		u(I)=6.0_DP			!upper bounds
	ENDDO

	DO WHILE(.TRUE.)
		CALL SETULB(N,M,X,L,U,NBD,F,G,FACTR,PGTOL,WA,IWA,TASK,IPRINT,CSAVE,LSAVE,ISAVE,DSAVE)
		IF(TASK(1:2) .EQ. 'FG') THEN
			DO I=1,N
				XV(I)=X(I)
			ENDDO
			FVC=FV(XV)
			F=SUM(FVC)
			CALL FDJAC_D(XV,FVC,DF)
			DO I=1,N
				G(I)=SUM(DF(:,I))
			ENDDO
			CYCLE
		ENDIF

		IF(TASK(1:5) .EQ. 'NEW_X') THEN
		!the minimization routine has returned with a new iterate.
		!At this point have the opportunity of stopping the iteration 
		!or observing the values of certain parameters
			CYCLE
		ELSE
			EXIT
		ENDIF
	ENDDO

	WRITE(*,*) X(1),X(2)
	CLOSE(7)
END SUBROUTINE EITEST

!??ƽ??????ʽ?????з????г?
FUNCTION FV(X)
	USE GLOBAL
	USE NRTYPE
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: X
	REAL(DP), DIMENSION(SIZE(X)) :: FV
	INTEGER	::I

	FV(1)=-x(1)**3+5*x(1)**2-x(1)+2*x(2)-3.0_dp
	FV(2)=x(2)**3+x(2)**2-14*x(2)-x(1)-19.0_dp

	!Ϊ???̼?ƽ??
	DO I=1,SIZE(X)
		FV(I)=FV(I)**2
	ENDDO
ENDFUNCTION FV