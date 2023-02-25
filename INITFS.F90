SUBROUTINE INITFS()
	USE FSBR
	IMPLICIT NONE
	REAL(DP),POINTER					:: HX(:),AREA(:)
	REAL(DP),POINTER					:: XC(:,:),YC(:,:)
	REAL(DP),POINTER					:: HY
	LOGICAL		:: SAVESLIP
	INTEGER		:: NX,NY
	INTEGER		:: I,J
	REAL(DP)	:: TEMP

!!!!!!!!请勿改动一下代码 完全正确
	NX= FS%NX
	NY= FS%NY
	HY=>FS%HY

	ALLOCATE(		FS%HX	 (-NY:NY)); 
	ALLOCATE(		FS%AREA(-NY:NY)); 
	ALLOCATE(		FS%XC(-NY-1:NY+1,-NX-1:NX+1));	
	ALLOCATE(		FS%YC(-NY-1:NY+1,-NX-1:NX+1));	

	ALLOCATE(	FS%TX(-NY-1:NY+1,-NX-1:NX+1));	
	ALLOCATE(	FS%TY(-NY-1:NY+1,-NX-1:NX+1));	
	ALLOCATE(	FS%SX(-NY-1:NY+1,-NX-1:NX+1));	
	ALLOCATE(	FS%SY(-NY-1:NY+1,-NX-1:NX+1));	
	HX	=>FS%HX
	AREA=>FS%AREA
	XC	=>FS%XC
	YC	=>FS%YC
	!此处不宜初始化任何输入参数
	!计算所需的压力作用点的坐标，以及网格的尺寸
	HY=1.0_DP/NY
	XC( NY+1,:)=	0.0_DP		!(0, 1)坐标
	YC( NY+1,:)=	1.0_DP		!(0, 1)坐标
	XC(-NY-1,:)=	0.0_DP		!(0,-1)坐标
	YC(-NY-1,:)=   -1.0_DP	!(0,-1)坐标
	DO I=1,NY
		TEMP=(I-0.5_DP)*HY
		HX(I) =SQRT(1.0_DP-TEMP**2)/NX
		HX(-I)=HX(I)
		YC(I,:)=TEMP    !计算压力点坐标
		YC(-I,:)=-TEMP
		DO J=1,NX
			TEMP=(J-0.5_DP)*HX(I)
			XC(I,J)=TEMP
			XC(I,-J)=-TEMP
			XC(-I,J)=TEMP
			XC(-I,-J)=-TEMP
		ENDDO
		XC(I, NX+1)	=XC( I, NX)+HX(I)*0.5_DP
		XC(I,-NX-1)	=XC( I,-NX)-HX(I)*0.5_DP
		XC(-I,NX+1)	=XC(I, NX+1)
		XC(-I,-NX-1)=XC(I,-NX-1)
		AREA(I)=HY*HX(I)
		AREA(-I)=AREA(I)
	ENDDO

	AREA(0)=HUGE(1.0_DP)
	HX(0)  =HUGE(1.0_DP)
	XC(0,:)=HUGE(1.0_DP)
	XC(:,0)=HUGE(1.0_DP)
	YC(0,:)=HUGE(1.0_DP)
	YC(:,0)=HUGE(1.0_DP)
END SUBROUTINE INITFS

