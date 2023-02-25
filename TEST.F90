!FASTSIM测试程序
SUBROUTINE FTEST(NO,VX,VY,PHI)
	USE FSBR
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: NO
	REAL(DP) ::VX,VY,PHI,G1,LOAD,MU,A,B,C11,C22,C23
	REAL(DP) ::L1,L2,L3,NP,Z0,UX,UY,PHX,PHY
	REAL(DP)	::TX,TY
	LOGICAL	:: INOUT

	G1=8.4E10_DP		!测试程序
	LOAD=1.0E5_DP
	MU=0.3_DP
	IF(NO==1) THEN
		A=6.0E-3_DP		!测试程序
		B=6.0E-3_DP		!测试程序
		C11=4.12_DP		!测试程序
		C22=3.67_DP		!测试程序
		C23=1.47_DP		!测试程序
	ELSEIF(NO==2) THEN
		A=7.5E-3_DP		!测试程序
		B=1.5E-3_DP		!测试程序
		C11=7.78_DP				!测试程序
		C22=8.14_DP				!测试程序
		C23=6.63_DP				!测试程序
	ELSE
		A=1.5E-3_DP		!测试程序
		B=7.5E-3_DP		!测试程序
		C11=3.37_DP				!测试程序
		C22=2.63_DP				!测试程序
		C23=0.603_DP				!测试程序
	ENDIF

	L1=8.0_DP*A/3.0_DP/C11/G1
	L2=8.0_DP*A/3.0_DP/C22/G1
	L3=PI_D*A**2/4.0_DP/G1/SQRT(A*B)/C23
	NP=2._DP*PI_D/3._DP     !采用Hertz接触进行计算
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
	!FS% TRACTION	=.TRUE.			!无量钢拖动力   更精确!
	FS% HEATHCOTE	=.FALSE.			!HEATHCOTE SLIP   FALSE
	FS% HERTZ			=.TRUE.				!Hertz Contact Model

	FS%UX =UX
	FS%UH =0.0_dp
	FS%UY =UY
	FS%PHX=PHX
	FS%PHY=PHY

	CALL FASTSIM()

	INOUT=.FALSE.
	CALL POSTFS(1,INOUT)

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

SUBROUTINE FTEST2(ax,by,VX,VY,PHI)
	USE FSBR
	IMPLICIT NONE
	REAL(DP) ::ax,by,VX,VY,PHI

	FS% HEATHCOTE	=.true.		!HEATHCOTE SLIP   FALSE
	FS% HERTZ		=.TRUE.			!Hertz Contact Model
	BR% HERTZ		=FS% HERTZ		!Hertz Contact Model
	BR% SIMPMOMENT	=.TRUE.			!采用简化方法计算作用于球心的力矩
	BR% SYM			=.TRUE.			!是否对称载荷

	BR%G1=8.4E10_DP		
	BR%G2=8.4E10_DP		
	br%nu1=0.25_dp
	BR%BALNO(:)%OUT%P=1.0E5_DP
	BR%BALNO(:)%INN%P=1.0E5_DP
	BR%MU=0.3_DP
	BR%BALNO(:)%OUT%AX=ax
	BR%BALNO(:)%OUT%BY=by		
	BR%BALNO(:)%OUT%AB=AX*BY		
	BR%BALNO(:)%INN%AX=ax
	BR%BALNO(:)%INN%BY=by
	BR%BALNO(:)%INN%AB=AX*BY

	if(ax>=by) then 
		br%icom%xlarger=.true.
		BR%OCOM%KE=by/ax		
		BR%ICOM%KE=by/ax
	else	
		br%icom%xlarger=.false.
		BR%OCOM%KE=ax/by		
		BR%ICOM%KE=ax/by
	endif	
	call getcij(0.0_dp,0.0_dp)  

	BR%BALNO(:)%OUT%XIX=VX		
	BR%BALNO(:)%OUT%XIY=VY		
	BR%BALNO(:)%OUT%PSI=PHI		
	BR%BALNO(:)%INN%XIX=VX		
	BR%BALNO(:)%INN%XIY=VY		
	BR%BALNO(:)%INN%PSI=PHI		

	call getcreep_idx(1.0_dp,1.0_dp,1)	!计算fastsim所需的各种参数

	BR%BALNO(:)%OUT%UH=0.0000000_DP
	BR%BALNO(:)%INN%UH=0.000_DP

	FS%UX =BR%BALNO(1)%OUT%UX
	FS%UH =BR%BALNO(1)%OUT%UH
	FS%UY =BR%BALNO(1)%OUT%UY
	FS%PHX=BR%BALNO(1)%OUT%PHX
	FS%PHY=BR%BALNO(1)%OUT%PHY
	CALL FASTSIM()
	CALL POSTFS(1,.TRUE.)

	br%tmp1=ax*by*br%g1/(3.0_dp*BR%BALNO(1)%OUT%P*br%mu)*br%ocom%c11*vx
	br%tmp2=ax*by*br%g1/(3.0_dp*BR%BALNO(1)%OUT%P*br%mu)*br%ocom%c22*vy
	br%aj  =(ax*by)**1.5_dp*br%g1*br%ocom%c23/(BR%BALNO(1)%OUT%P*br%mu)*phi
	br%ai=BR%BALNO(1)%OUT%TX/BR%BALNO(1)%OUT%P/br%mu
	br%ao=BR%BALNO(1)%OUT%TY/BR%BALNO(1)%OUT%P/br%mu
!	write(*,*) br%tmp1,br%tmp2,br%aj
!	WRITE(*,*) BR%BALNO(1)%OUT%TX,BR%BALNO(1)%OUT%TY
!	WRITE(*,*) BR%BALNO(1)%OUT%TX/BR%BALNO(1)%OUT%P/br%mu,BR%BALNO(1)%OUT%TY/BR%BALNO(1)%OUT%P/br%mu
ENDSUBROUTINE FTEST2

SUBROUTINE FTEST3(VX,VY,PHI)
	USE FSBR
	IMPLICIT NONE
	REAL(DP) ::ax,by,VX,VY,PHI

	ax=BR%BALNO(1)%OUT%AX
	by=BR%BALNO(1)%OUT%BY		

	BR%BALNO(:)%OUT%XIX=VX		
	BR%BALNO(:)%OUT%XIY=VY		
	BR%BALNO(:)%OUT%PSI=PHI		
	BR%BALNO(:)%INN%XIX=VX		
	BR%BALNO(:)%INN%XIY=VY		
	BR%BALNO(:)%INN%PSI=PHI		

	call getcij(0.0_dp,0.0_dp)  
	call getcreep_idx(1.0_dp,1.0_dp,1)	!计算fastsim所需的各种参数

	FS%UX =BR%BALNO(1)%OUT%UX
	FS%UH =BR%BALNO(1)%OUT%UH
	FS%UY =BR%BALNO(1)%OUT%UY
	FS%PHX=BR%BALNO(1)%OUT%PHX
	FS%PHY=BR%BALNO(1)%OUT%PHY
	CALL FASTSIM()
	CALL POSTFS(1,.TRUE.)

	br%tmp1=ax*by*br%g1/(3.0_dp*BR%BALNO(1)%OUT%P*br%mu)*br%ocom%c11*vx
	br%tmp2=ax*by*br%g1/(3.0_dp*BR%BALNO(1)%OUT%P*br%mu)*br%ocom%c22*vy
	br%aj  =(ax*by)**1.5_dp*br%g1*br%ocom%c23/(BR%BALNO(1)%OUT%P*br%mu)*phi
	br%ai	=BR%BALNO(1)%OUT%TX/BR%BALNO(1)%OUT%P/br%mu
	br%ao	=BR%BALNO(1)%OUT%TY/BR%BALNO(1)%OUT%P/br%mu
!	write(*,*) br%tmp1,br%tmp2,br%aj
!	WRITE(*,*) BR%BALNO(1)%OUT%TX,BR%BALNO(1)%OUT%TY
!	WRITE(*,*) BR%BALNO(1)%OUT%TX/BR%BALNO(1)%OUT%P/br%mu,BR%BALNO(1)%OUT%TY/BR%BALNO(1)%OUT%P/br%mu
ENDSUBROUTINE FTEST3

!分析摩擦系数对
SUBROUTINE FCOE()
	USE FSBR
	USE INTFACE,ONLY:POW
	IMPLICIT NONE
	CHARACTER(len=1)			:: tab=CHAR(9)

	INTEGER	:: I,N

	FS% HEATHCOTE	=.FALSE.		!HEATHCOTE SLIP   FALSE
	FS% HERTZ		=.TRUE.			!Hertz Contact Model
	BR% HERTZ		=FS% HERTZ		!Hertz Contact Model
	BR% SIMPMOMENT	=.TRUE.			!采用简化方法计算作用于球心的力矩
	BR% SYM			=.TRUE.			!是否对称载荷

	OPEN(5, file="t.XLS")  !创建数据文件 EXCEL文件
	
	N=100
	DO I=1,N
		BR%MU=1.0E-10_dp*pow(10.0_dp,1.0_dp*I/10)
!		BR%MU=0.3_dp
		CALL GETCREEP2()
		CALL PREFS(1,NOT(.TRUE.))
		CALL FASTSIM()
		CALL POSTFS(1,NOT(.TRUE.))

!		CALL VECTOR_COMBINING(FS)  !矢量合成x、y，并将结果存入前一变量
!		CALL PRINT_FIELD(FS)
	!	CALL PRINT_VFIELD(FS)
!		WRITE(*,*) BR%BALNO(1)%INN%TX,TAB,BR%BALNO(1)%INN%TY
		WRITE(5,"(e23.15,A1,e23.15,A1,e23.15)") BR%MU,TAB,BR%BALNO(1)%INN%TX,TAB,BR%BALNO(1)%INN%TY
	ENDDO
	CLOSE(5)
ENDSUBROUTINE FCOE

!分析自旋角速度对横向拖动力的影响
SUBROUTINE SPINY()
	USE FSBR
	USE INTFACE,ONLY:POW
	IMPLICIT NONE
	CHARACTER(len=1)			:: tab=CHAR(9)

	INTEGER	:: I,N
	real(dp)	:: phi

	OPEN(5, file="spiny.dat")  !创建数据文件 EXCEL文件

!	BR%BALNO(:)%OUT%P =4.50E5_DP
!	BR%MU =0.1_dp 
	CALL GETPCHAXIS(pi_d*.2,pi_d*.2)
!	BR%BALNO(:)%OUT%AX =1.8E-3_DP
!	BR%BALNO(:)%OUT%BY =6.5E-3_DP
!	BR%OCOM%KE =BR%BALNO(1)%OUT%AX/BR%BALNO(1)%OUT%BY
!	BR%BALNO(:)%OUT%AB =BR%BALNO(1)%OUT%AX*BR%BALNO(1)%OUT%BY

		BR%BALNO(:)%OUT%XIX=0.00_DP
		BR%BALNO(:)%OUT%XIY=0.00_DP
	N=1000
	DO I=1,N
		phi=1.E-2_DP*i
		BR%BALNO(:)%OUT%PSI=phi/(BR%BALNO(1)%OUT%AB**1.5_dp*br%g1*br%ocom%c23/(br%mu*BR%BALNO(1)%OUT%P))
!		BR%BALNO(:)%OUT%PSI=-.1_DP*i
!		CALL GETCREEP(BR)		!计算FASTSIM所需的各种参数
		CALL GETCREEP2()		!计算FASTSIM所需的各种参数
		CALL PREFS(1,.TRUE.)  
		CALL FASTSIM()
		CALL POSTFS(1,.TRUE.)
!	CALL PRINT_SFIELD(FS)

		WRITE(5,"(e23.15,A1,e23.15)") phi,TAB,-FS%TFY	
	ENDDO
	CLOSE(5)
ENDSUBROUTINE SPINY

SUBROUTINE ELLIP()
	USE FSBR
	USE INTFACE,ONLY:POW
	IMPLICIT NONE
	TYPE(PATCHABK)	 ,TARGET	:: ABK
	CHARACTER(len=1)			:: tab=CHAR(9)
	REAL(DP)	:: P,P13,P23,A,B

	P=4.45_DP
	ABK%RX1		=6.35E-3_DP
	ABK%RY1		=6.35E-3_DP
	ABK%RX2		=-6.60E-3_DP
	ABK%RY2		=-38.9E-3_DP

	ABK%KP		=BR%KP
	ABK%E1		=2.28E11_DP
	ABK%E2		=2.28E11_DP
	ABK%NU1		=BR%NU1
	ABK%NU2		=BR%NU2
	!计算内圈接触面椭圆率、半轴
	CALL GETABK(ABK)

	p13=POW(P,1.0_dp/3.0_DP)
	p23=POW(P,2.0_dp/3.0_DP)
	A=ABK%AND	*p13
	B=ABK%BND	*p13
write(*,*) 'a',a,'b',b,'K',ABK%KE
!	ABK%KE
!	ABK%AND
!	ABK%BND
!	ABK%APPND
!	ABK%EI1K
!	ABK%EI2E
!	ABK%RX
!	ABK%RY
!	ABK%RE
!	ABK%XLARGER

ENDSUBROUTINE ELLIP