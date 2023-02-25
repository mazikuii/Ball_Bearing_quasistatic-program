!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!            FASTSIM ALGORITH       !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!Here in this surboutine namely FASTSIM hereafter, the FASTSIM ALGORITH 
!based on the simplied theory developed by J.J. Kalker, perhaps the most 
!famous resercher in the field wheel/rail contact and dynamic analysis, is 
!encapsulated in a single subroutine, as compared with the original one 
!probably programed by Prof. Kalker himself, to facilite its implementation 
!by other programs or routines.
!
!In addition to the original copy of FASTSIM, some modifications are made to 
!adapt it for more general applications including:
!*
!*
!input arguments:
!NX_F		: numbers of divided segments of semi-axis a, parallel with the rolling direction. 
!NY_F		: numbers of divided segments of semi-axis b, perpendicular to the rolling direction. 
!HX_F(:)	: array of segment length of in x direction
!HY_F		: segment length of y direction
!UX_F,UY_F	: dimensionless creepages in x and y direction
!PH_F		: simensionless spin, perpendicular to the contact patch
!output arguments:
!TX_F,TY_F	: traction components in x and y direction caused by creepages and spin within the contact patch
!MZX_F,MZX_F: moments perpendicular to the contact patch about the center of contact patch caused by the creep 
!             force components in x and y direction respectively.
!ORIFIELD_F(:,:)	: orientation field of slips within the contact patch, anti-clockwisely starting from the position direction of x
SUBROUTINE FASTSIM()			   !outputs
	USE FSBR
	IMPLICIT NONE
	!subroutine intput and output arguments
	!local arguments
	INTEGER	,POINTER :: NX,NY
	LOGICAL ,POINTER :: HERTZ,HEATHCOTE
	REAL(DP),POINTER :: HY,HX(:),SPOLE(:,:),AREA(:)
	REAL(DP),POINTER :: XC(:,:),YC(:,:),TX(:,:),TY(:,:),SX(:,:),SY(:,:)
	REAL(DP),POINTER :: UX,UH,UY,PHX,PHY
	REAL(DP)	:: P,G,CX,CY,DX,TEMP,X1,Y1,X0,Y0
	INTEGER		:: I,J

	NX					=>FS%NX
	NY					=>FS%NY
	HY					=>FS%HY
	HERTZ				=>FS%HERTZ
	HEATHCOTE			=>FS%HEATHCOTE

	UX					=>FS%UX
	UH					=>FS%UH
	UY					=>FS%UY
	PHX					=>FS%PHX
	PHY					=>FS%PHY

	HX					=>FS%HX
	AREA				=>FS%AREA
	XC					=>FS%XC
	YC					=>FS%YC
	TX					=>FS%TX
	TY					=>FS%TY
	SX					=>FS%SX
	SY					=>FS%SY
	SPOLE				=>FS%SPOLE
	!!!以下程序用来计算拖动力分布
	!CALCULATE SPIN POLE OF CONTACT PATCH.
!	IF(PHY==0.0_dp) THEN;SPOLE(1,1)=SIGN(HUGE(1.0_DP),-UY);ELSE;SPOLE(1,1)=-UY/PHY;ENDIF
!	IF(PHX==0.0_dp) THEN;SPOLE(1,2)=SIGN(HUGE(1.0_DP), UX);ELSE;SPOLE(1,2)= UX/PHX;ENDIF
	DX=0.0_DP
	TX(:,:)=0.0_DP  
	TY(:,:)=0.0_DP  
	SX(:,:)=0.0_DP  
	SY(:,:)=0.0_DP  
	DO I=NY,-NY,-1					!求解第J层网格，Y方向
	  IF(I/=0) THEN				!跳过0索引点
			X0=0.0_DP; Y0=0.0_DP
			DX=HX(I)
			DO J=NX,-NX,-1		!递推、积分    !此处边界很重要
				IF(J/=0) THEN									!跳过0索引点
					IF(HEATHCOTE) THEN
						CX=UX-YC(I,J)*PHX+YC(I,J)**2*UH  !UH必须为正
					ELSE
						CX=UX-YC(I,J)*PHX   !WITHOUT HEATHCOTE SLIP
					ENDIF
					CY=UY+XC(I,J)*PHY															 !RIGID SLIP
					IF(J==NX) THEN	   
						X1=-CX*DX*0.5_DP   
						Y1=-CY*DX*0.5_DP   
					ELSE
						X1=X0-CX*DX   
						Y1=Y0-CY*DX   
					ENDIF
					G=1.0_DP-XC(I,J)**2-YC(I,J)**2				 !Kalker接触 1-x*x-y*y
					IF(HERTZ) THEN
						G=SQRT(G)											!Hertz contact  (1-x*x-y*y)**1/2
					ENDIF
					TEMP=SQRT(X1**2+Y1**2)  
					P=TEMP/G    !在接触后沿，P_F=0/0，极端情况
					IF(P>1.0_DP) THEN  !滑动区
						!滑移速度分布
						SX(I,J)=-X1/DX*(1.0_DP-G/TEMP) !无量纲滑移速度，否则为零!必须适用变换前的压力值
						SY(I,J)=-Y1/DX*(1.0_DP-G/TEMP) !无量纲滑移速度，否则为零
						X1=X1/P 
						Y1=Y1/P 
					ENDIF
					TX(I,J)=X1
					TY(I,J)=Y1
					X0=X1
					Y0=Y1
				ENDIF
			ENDDO
			CX=UX-YC(I,J)*PHX
			CY=UY+XC(I,J)*PHY
			SX(I,J)=CX-X0/DX
			SY(I,J)=CY-Y0/DX
	  ENDIF
	ENDDO
sx(:,-NX-1)=0.0_dp
sy(:,-NX-1)=0.0_dp
ENDSUBROUTINE FASTSIM