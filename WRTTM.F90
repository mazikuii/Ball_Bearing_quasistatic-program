SUBROUTINE WRTTM(idx)
	USE FSBR
	USE INTFACE,ONLY:RAD2DEG,POW
	IMPLICIT NONE
	integer	::idx
	
	WRITE(*,*) "自由接触角:",BR%ALf0 ,RAD2DEG(BR%ALf0)
	WRITE(*,*) "实际接触角:",RAD2DEG(BR%balno(idx)%inn%alf),RAD2DEG(BR%balno(idx)%out%alf)
	WRITE(*,*) "静平衡接触角:",BR%BALNO(idx)%INN%ALfST,RAD2DEG(BR%BALNO(idx)%INN%ALfST)
	!WRITE(*,*) "无滑动外圈角速度:",OMGO,BR%OMGI*(1.0_DP-GAMA)/(1.0_DP+GAMA)
	!WRITE(*,*) "无滑动滚动体角速度:",OMGB,BR%OMGI*BR%DM/BR%D*(1.0_DP-GAMA)

	WRITE(*,*) "TX=",BR%BALNO(idx)%OUT%TX
	WRITE(*,*) "TY=",BR%BALNO(idx)%OUT%TY
	WRITE(*,*) "PO=",BR%BALNO(idx)%OUT%P
	WRITE(*,*) "MX=",BR%BALNO(idx)%OUT%MX
	WRITE(*,*) "MY=",BR%BALNO(idx)%OUT%MY
	WRITE(*,*) "MZ=",BR%BALNO(idx)%OUT%MZ
	WRITE(*,*) "拖动力合力:",SQRT(BR%BALNO(idx)%OUT%TX**2+BR%BALNO(idx)%OUT%TY**2)
	WRITE(*,*) "外接触滑动COULOMB力:",BR%MU*BR%BALNO(idx)%OUT%P
!	write(*,*) '椭圆率:',	BR%OCOM%KE,	SQRT(BR%OCOM%KE),1.0_DP/BR%OCOM%KE
	WRITE(*,*) 
	WRITE(*,*) "TX=",BR%BALNO(idx)%INN%TX
	WRITE(*,*) "TY=",BR%BALNO(idx)%INN%TY
	WRITE(*,*) "PI=",BR%BALNO(idx)%INN%P
	WRITE(*,*) "MX=",BR%BALNO(idx)%INN%MX
	WRITE(*,*) "MY=",BR%BALNO(idx)%INN%MY
	WRITE(*,*) "MZ=",BR%BALNO(idx)%INN%MZ
	write(*,*) "拖动力合力:",sqrt(BR%BALNO(idx)%INN%TX**2+BR%BALNO(idx)%INN%TY**2)
	WRITE(*,*) "内接触滑动COULOMB力:", BR%MU*BR%BALNO(idx)%INN%P
!	write(*,*) '椭圆率:',	BR%iCOM%KE,	SQRT(BR%iCOM%KE),1.0_DP/BR%iCOM%KE
	WRITE(*,*) 
	write(*,*) "滚动体离心力:",br%balno(idx)%fc
	write(*,*) "滚动体x轴陀螺力矩分量:",br%balno(idx)%mgx
	write(*,*) "滚动体z轴陀螺力矩分量:",br%balno(idx)%mgz
ENDSUBROUTINE WRTTM

SUBROUTINE OPTSLU()
	USE FSBR
	USE INTFACE,ONLY:RAD2DEG,rad2rev
	IMPLICIT NONE
	REAL(DP)::PCHI,PCHO,omgb,omgo
	real(dp)::t1,t2
	character(len=1)			:: tab=char(9)
	integer	:: i
	real(dp) DNRM2
	
	!输出到文件，便于提取数据
	write(12,*) '输出由JONES的计算方法得到的初解'
	write(12,*) BR%x0(1)*60_dp/(2.0_dp*PI_D),tab,"OMGO RPM"
	write(12,*) BR%x0(2)*60_dp/(2.0_dp*PI_D),tab,"OMGX RPM"
	write(12,*) BR%x0(3)*60_dp/(2.0_dp*PI_D),tab,'OMGY RPM'
	write(12,*) BR%x0(4)*60_dp/(2.0_dp*PI_D),tab,'OMGZ RPM'
	write(12,*) RAD2DEG(BR%x0(5)),tab,"ALPHAI"
	write(12,*) RAD2DEG(BR%x0(6)),tab,'ALPHAO'
	write(12,*) sqrt(BR%x0(3)**2+BR%x0(4)**2)*60_dp/(2.0_dp*PI_D),tab,'OMGB  RPM滚动体自转角速度'	!
	write(12,*) RAD2DEG(BR%BALNO(1)%INN%ALfST),tab, 'ALPHA DEGREE'
	write(12,*) 0.0_dp,tab, 'beta DEGREE'
	write(12,*)	'俯仰角JONES，近似计算'
	write(12,*)	RAD2DEG(br%pchang_octl),tab,'外滚道控制俯仰角HARRIS,进一步近似计算'
	write(12,*)	RAD2DEG(br%pchang_ictl),tab,'内滚道控制俯仰角HARRIS,进一步近似计算'
	write(12,*) br%omgb_octl,tab,'外滚道控制滚动体转速RPM'
	write(12,*) br%omgo_octl,tab,'外滚道控制外圈转速RPM'
	write(12,*) br%s2rinn,tab,'外滚道控制，内滚道旋滚比'
	write(12,*) br%omgb_ictl,tab,'内滚道控制滚动体转速RPM'
	write(12,*) br%omgo_ictl,tab,'内滚道控制外圈转速RPM'
	write(12,*) br%s2rout,tab,'内滚道控制，外滚道旋滚比'

	write(12,*) '输出方程的最终解'
	write(12,*) rad2rev(BR%x(2)),tab,'OMGX  RPM X轴分量(相对于保持架)'
	write(12,*) rad2rev(BR%x(3)),tab,'OMGY  RPM Y轴分量(相对于保持架)'
	write(12,*) rad2rev(BR%x(4)),tab,'OMGZ  RPM Z轴分量(相对于保持架)'
write(12,*) br%mu,tab,'br%mu'
write(12,*) br%aload,tab,'br%aload'
	if(br%absolute) then
		write(12,*) rad2rev(BR%x(1)),tab,'OMGC  RPM 保持架角速度'
	else
		write(12,*) rad2rev(BR%x(1)),tab,'OMGO  RPM 外圈角速度'
	endif
	write(12,*) RAD2DEG(BR%x(5)),tab,'ALPHAI DEGREE 内接触角'
	write(12,*) RAD2DEG(BR%x(6)),tab,'ALPHAO DEGREE 外接触角'
	write(12,*) RAD2DEG(atan(abs(BR%x(4))/SQRT(BR%x(2)**2+BR%x(3)**2))),tab,'俯仰角ALPHA DEGREE'
	write(12,*) RAD2DEG(atan(BR%x(2)/abs(BR%x(3)))),tab,'偏航角BETA DEGREE'
	write(12,*) sqrt(rad2rev(BR%x(2))**2+rad2rev(BR%x(3))**2+rad2rev(BR%x(4))**2),tab, '滚动体自转角速度'
	write(12,*) BR%BALNO(1)%out%spn2rol,tab,'外接触旋滚比'
	write(12,*) BR%BALNO(1)%inn%spn2rol,tab,'内接触旋滚比'
	WRITE(12,*) BR%BALNO(1)%OUT%TX,	tab,"TX"
	WRITE(12,*) BR%BALNO(1)%OUT%TY,	tab,"TY"
	WRITE(12,*) BR%BALNO(1)%OUT%P,	tab,"PO"
	WRITE(12,*) BR%BALNO(1)%OUT%MX,	tab,"MX"
	WRITE(12,*) BR%BALNO(1)%OUT%MY,	tab,"MY"
	WRITE(12,*) BR%BALNO(1)%OUT%MZ,	tab,"MZ"
	WRITE(12,*) SQRT(BR%BALNO(1)%OUT%TX**2+BR%BALNO(1)%OUT%TY**2),tab,"拖动力合力"
	WRITE(12,*) BR%BALNO(1)%INN%TX,	tab,	"TX"
	WRITE(12,*) BR%BALNO(1)%INN%TY,	tab,	"TY"
	WRITE(12,*) BR%BALNO(1)%INN%P,	tab,	"PI"
	WRITE(12,*) BR%BALNO(1)%INN%MX,	tab,	"MX"
	WRITE(12,*) BR%BALNO(1)%INN%MY,	tab,	"MY"
	WRITE(12,*) BR%BALNO(1)%INN%MZ,	tab,	"MZ"
	write(12,*) sqrt(BR%BALNO(1)%INN%TX**2+BR%BALNO(1)%INN%TY**2),tab,"拖动力合力"
	write(12,*) '******************************************'
	WRITE(12,*) RAD2DEG(BR%ALf0),tab,"自由接触角"
	WRITE(12,*) RAD2DEG(BR%X(5)),tab,"实际内接触角"
	WRITE(12,*) RAD2DEG(BR%X(6)),tab,"实际外接触角"
	WRITE(12,*) RAD2DEG(BR%BALNO(1)%INN%ALfST),tab,"静平衡接触角"
	!WRITE(*,*) "无滑动外圈角速度:",OMGO,BR%OMGI*(1.0_DP-GAMA)/(1.0_DP+GAMA)
	!WRITE(*,*) "无滑动滚动体角速度:",OMGB,BR%OMGI*BR%DM/BR%D*(1.0_DP-GAMA)
	WRITE(12,*) BR%MU*BR%BALNO(1)%OUT%P,tab,"外接触滑动COULOMB力"
!	write(12,*) '椭圆率:',	BR%OCOM%KE,	SQRT(BR%OCOM%KE),1.0_DP/BR%OCOM%KE
	WRITE(12,*) 
	WRITE(12,*) BR%MU*BR%BALNO(1)%INN%P,tab, "内接触滑动COULOMB力"
!	write(12,*) '椭圆率:',	BR%iCOM%KE,	SQRT(BR%iCOM%KE),1.0_DP/BR%iCOM%KE
	WRITE(12,*) 
	write(12,*) br%balno(1)%fc,	tab,"滚动体离心力"
	write(12,*) br%balno(1)%mgx,tab,"滚动体x轴陀螺力矩分量"
	write(12,*) br%balno(1)%mgz,tab,"滚动体z轴陀螺力矩分量"
	write(12,*) 
	do i=1,br%n
		write(12,"(e23.15,a1,a5,i1,a1)") br%balno(1)%f(i),tab,'br%f(',i,')'
	enddo
	write(12,*) DNRM2(size(br%balno(1)%f),br%balno(1)%f,1),tab,'norm of F(x)'

ENDSUBROUTINE OPTSLU

SUBROUTINE OPTLOG()	!输出输入新信息
	USE FSBR
	IMPLICIT NONE
	character(len=1)			:: tab=char(9)
	
	!输出程序求解信息，便于再现计算过程
	write(12,*) '**************程序输入数据*********************'
	write(12,*)  pntxo,					tab,'pntxo'
	write(12,*)  maxit,					tab,'maxit'
	write(12,*)  embed,					tab,'embed'
	write(12,*)  mprec,					tab,'mprec'
	write(12,*)	 sqmprec,				tab,'sqmprec'
	write(12,*)  FS% HEATHCOTE,	tab,'FS%HEATHCOTE'
	write(12,*)  FS% HERTZ	,		tab,'FS%HERTZ'
	write(12,*)  FS%	NY,				tab,'FS%NY'
	write(12,*)  FS%	NX,				tab,'FS%NX'
	write(12,*)  br%simpmoment,	tab,'br%simpmoment'
 	write(12,*)  br%sym,				tab,'br%sym'
	write(12,*)  br%hertz,			tab,'br%hertz'						!hertz压力分布控制字，true:hertz；false:non-hertz
	write(12,*)  br%OPTIN,			tab,'br%OPTIN'
	write(12,*)  br%unscl,			tab,'br%unscl'	!unscale coordinates? TRUE:UNSCALE;FALSE:NO
	write(12,*)  br%rcm,				tab,'br%rcm'		!raceway control method
	write(12,*)  br%absolute,		tab,'br%absolute'		!绝对运动,所有角速度均相对于惯性坐标系而言
	write(12,*)  br%cmdcreep,		tab,'br%cmdcreep'	!联合载荷作用下的蠕滑模型
	write(12,*)  br%hcmthd ,		tab,'br%hcmthd' !heathcote计算方法
	write(12,*)  br%z		,				tab,'br%z'		!滚动体个数
	write(12,*)  br%cyc,				tab,'br%cyc'
	write(12,*)  br%inert,			tab,'br%inert'	!0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
	write(12,*)  br%d,					tab,'br%d'
	write(12,*)  br%r		,				tab,'br%r'	!diameter of balls直径、半径
	write(12,*)  br%ri,					tab,'br%ri'
	write(12,*)  br%ro	,				tab,'br%ro'	!inner and outer raceway groove curvature radius
	write(12,*)  br%d_i,				tab,'br%d_i'
	write(12,*)  br%d_o	,				tab,'br%d_o'!inner and outer raceway diameter
	write(12,*)  br%kp,					tab,'br%kp'
	write(12,*)  br%kpst,				tab,'br%kpst'	!椭圆率、静态接触角计算精度
	write(12,*)  br%e1,					tab,'br%e1'
	write(12,*)  br%g1	,				tab,'br%g1'	!滚动体弹性、剪切模量
	write(12,*)  br%e2,					tab,'br%e2'
	write(12,*)  br%g2	,				tab,'br%g2'	!套圈弹性、剪切模量
	write(12,*)  br%mu	,				tab,'br%mu'			!摩擦系数
	write(12,*)  br%nu1	,				tab,'br%nu1'		!滚动体弹性模量
	write(12,*)  br%nu2	,				tab,'br%nu2'		!滚动体弹性模量
	write(12,*)  br%rho1,				tab,'br%rho1'			!密度
	write(12,*)  br%rho2	,			tab,'br%rho2'		!密度
	write(12,*)  br%aload	,			tab,'br%aload'	!轴承轴向载荷
	write(12,*)  br%rload	,			tab,'br%rload'	!轴承径向载荷
	write(12,*)  br%mload	,			tab,'br%mload'	!轴承径向载荷
	write(12,*)  br%omgi,				tab,'br%omgi'		!内圈角速度
	write(12,*)  br%omgo,				tab,'br%omgo'		!外圈角速度
	write(12,*)  br%startang,		tab,'br%startang' !起始滚动体的方位角
	write(12,*)  br%m,					tab,'br%m'
	write(12,*)  br%n,					tab,'br%n'
	write(12,*)  br%nsec,				tab,'br%nsec'
	write(12,*)  br%intmthd	,		tab,'br%intmthd'!1:接触角->变形量;2:变形量->接触角 忽略摩擦力
	write(12,*)  br%crpmthd	,		tab,'br%crpmthd'!1:接触角->变形量;2:变形量->接触角　考虑摩擦力
	write(12,*)  br%init,				tab,'br%init'
	write(12,*) '***********************'
ENDSUBROUTINE OPTLOG