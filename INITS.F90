subroutine initparas()
	use fsbr
	use intface,only:rev2rad
	implicit none
	integer :: bn

	br%n=8
	br%m=br%n
	br%nsec=20		!插值椭圆率等参数
	fs%heathcote		=.true.		!heathcote slip   false
	br%hcmthd				=3				!1-johnson最简单;2-johnson+mzk;3-mzk最复杂
	br%s2rmthd			=1		!1-不考虑套圈角速度；2-考虑套圈角速度
	fs%hertz				=.true.		!hertz contact model
	br%hertz				=fs%hertz	!hertz contact model
	br%simpmoment		=.true.		!采用简化方法计算作用于球心的力矩
	br%optin				=.false.		!true:output inner contact;false:outer
	br%unscl				=.true.		!true:unscale coordinates;false:unchanged,unit circle

	pntxo=.true.	!输出初值
	maxit=.false.

	!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
	br%inert				=4
	br%sym					=.true.		!是否对称载荷
	br%cmdcreep			=.true.	!联合载荷：T:蠕滑计算模型;F:理想计算模型
	br%rcm					=2	!1:mzk siplified;2:inner raceway control;3:outer raceway control

	br%absolute		=.false.	!true:绝对运动;false:相对运动
	br%intmthd			=2 !1:接触角->变形量;2:变形量->接触角
	br%crpmthd			=2 !1:接触角->变形量;2:变形量->接触角
	!接触角途径精度低、速度慢，但计算过程稳定，适用面广
	!变形量途径精度高、速度快，但不稳定，适用面较小，有待改进

	bn=7
	br%xtran=2.2_dp
	br%kp		=1.0e-16_dp
	br%kpst	=1.0e-16_dp
	mprec=epsilon(1.0_dp)
	sqmprec=sqrt(mprec)
	if(bn==1) then
		fs%	ny=20
		fs%	nx=20
		br%z		=16
		br%aload	=100._dp 
		br%omgi		=10000._dp
		br%mu			=.006_dp				!摩擦系数
		br%e1			=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !密度
		!no clearanceffr
		br%d		=6.35e-3_dp*2  
		br%d_i	=38.9e-3_dp*2-6.35e-3_dp*4		-0.e-3_dp	 
		br%d_o	=38.9e-3_dp*2 
		br%ri		=6.60e-3_dp 
		br%ro		=6.60e-3_dp 
	elseif(bn==2) then
		fs%	ny		=20
		fs%	nx		=20
		br%z			=16
		br%aload	=9000._dp
		if(br%absolute) then
			br%omgi		=4000.0_dp   !速度边界条件
			br%omgo		=0.00_dp   !速度边界条件	omgc未知
		else
			br%omgi		=200._dp		!omgc=0,omgo未知
		endif 
		br%mu			=.0001_dp				!摩擦系数
		br%e1			=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !密度
		!case from hamrock and dowson
		br%d		=6.35e-3_dp*2 
		br%d_i	=38.9e-3_dp*2-6.35e-3_dp*4-0.35e-3_dp	
		br%d_o	=38.9e-3_dp*2 
		br%ri		=6.60e-3_dp 
		br%ro		=6.60e-3_dp 
	elseif(bn==3) then
		fs%	ny		=50
		fs%	nx		=50
		br%z			=16
		br%aload	=1000._dp		!unit newton
		br%omgi		=560._dp			!rev/sec  !当转速较大时，滚动体有可能完全脱离内圈，从而形成负接触角
		br%mu			=.01_dp				!摩擦系数
		br%e1			=2.1e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.865e3_dp  !密度   kg/m_sup3
		!case from harris: angular contact ball bearing 218
		br%d		=22.23e-3_dp
		br%d_i	=102.79e-3_dp
		br%d_o	=147.73e-3_dp
		br%ri		=11.63e-3_dp
		br%ro		=11.63e-3_dp
	elseif(bn==4) then
		br%z			=16
		br%aload	=100._dp 
		br%omgi		=100._dp
		br%mu			=.001_dp				!摩擦系数
		br%e1			=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !密度
		!!!!!!!!!!!!!!
		br%d		=12.7e-3_dp
		br%d_i		=52.291e-3_dp
		br%d_o		=77.706e-3_dp
		br%ri		=6.6e-3_dp
		br%ro		=6.6e-3_dp
	elseif(bn==5) then
		br%z		=16
		br%aload	=100._dp 
		br%omgi		=100._dp
		br%mu			=.001_dp				!摩擦系数
		br%e1		=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !密度
		!case from 赵联春
		br%d		=15.081e-3_dp
		br%d_i	=49.912e-3_dp 
		br%d_o	=80.088e-3_dp
		br%ri		=7.665e-3_dp 
		br%ro		=8.01e-3_dp
	elseif(bn==6) then    !cjme 
		fs%	ny=200
		fs%	nx=200
		br%z		=3
		br%aload	=1000._dp		!unit newton
		br%omgi		=5000._dp		!rev/sec
		br%mu			=.03_dp				!摩擦系数
		br%e1		=2.1e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.865e3_dp  !密度   kg/m_sup3
		!!!!!!!!!!!!!!!!
		br%d		=022.2250e-3_dp
		br%d_i	=102.7938e-3_dp 
		br%d_o	=147.7264e-3_dp
		br%ri		=011.6281e-3_dp 
		br%ro		=011.6281e-3_dp
	elseif(bn==7) then
		fs%ny		=25
		fs%nx		=25
		br%z		=16

		br%aload		=40000._dp
		br%rload		=1000.0_dp
		br%mload		=100.0_dp
		br%startang	=0.0_dp
		if(br%absolute) then
			br%omgi		=2000._dp   !速度边界条件
			br%omgo		=0._dp   !速度边界条件	omgc未知
		else
			br%omgi		=7000._dp		!omgc=0,omgo未知
		endif
		br%mu				=0.1_dp				!摩擦系数
		!case from hamrock and dowson
		br%d		=022.2250e-3_dp
		br%d_i	=102.7938e-3_dp
		br%d_o	=147.7264e-3_dp
		br%ri		=011.6281e-3_dp
		br%ro		=011.6281e-3_dp
!		br%d		=6.35e-3_dp*2 
!		br%d_i	=38.9e-3_dp*2-6.35e-3_dp*4-0.45e-3_dp	
!		br%d_o	=38.9e-3_dp*2 
!		br%ri		=6.60e-3_dp 
!		br%ro		=6.60e-3_dp 
		br%e1				=2.06e11_dp
		br%nu1			=0.3_dp
		br%rho1			=7.8e3_dp  !密度
	endif

	br%omgi		=rev2rad(br%omgi)		!rad/sec
	br%omgo		=rev2rad(br%omgo)		!rad/sec

	br%e2			=br%e1
	br%nu2		=br%nu1
	br%rho2		=br%rho1  !密度
	if(br%e1/=br%e2 .or. br%nu1/=br%nu2 .or. br%rho1/=br%rho2) pause 'error in initparas!'
	if(br%m/=br%n) stop
endsubroutine initparas