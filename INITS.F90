subroutine initparas()
	use fsbr
	use intface,only:rev2rad
	implicit none
	integer :: bn

	br%n=8
	br%m=br%n
	br%nsec=20		!��ֵ��Բ�ʵȲ���
	fs%heathcote		=.true.		!heathcote slip   false
	br%hcmthd				=3				!1-johnson���;2-johnson+mzk;3-mzk���
	br%s2rmthd			=1		!1-��������Ȧ���ٶȣ�2-������Ȧ���ٶ�
	fs%hertz				=.true.		!hertz contact model
	br%hertz				=fs%hertz	!hertz contact model
	br%simpmoment		=.true.		!���ü򻯷����������������ĵ�����
	br%optin				=.false.		!true:output inner contact;false:outer
	br%unscl				=.true.		!true:unscale coordinates;false:unchanged,unit circle

	pntxo=.true.	!�����ֵ
	maxit=.false.

	!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
	br%inert				=4
	br%sym					=.true.		!�Ƿ�Գ��غ�
	br%cmdcreep			=.true.	!�����غɣ�T:�们����ģ��;F:�������ģ��
	br%rcm					=2	!1:mzk siplified;2:inner raceway control;3:outer raceway control

	br%absolute		=.false.	!true:�����˶�;false:����˶�
	br%intmthd			=2 !1:�Ӵ���->������;2:������->�Ӵ���
	br%crpmthd			=2 !1:�Ӵ���->������;2:������->�Ӵ���
	!�Ӵ���;�����ȵ͡��ٶ���������������ȶ����������
	!������;�����ȸߡ��ٶȿ죬�����ȶ����������С���д��Ľ�

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
		br%mu			=.006_dp				!Ħ��ϵ��
		br%e1			=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !�ܶ�
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
			br%omgi		=4000.0_dp   !�ٶȱ߽�����
			br%omgo		=0.00_dp   !�ٶȱ߽�����	omgcδ֪
		else
			br%omgi		=200._dp		!omgc=0,omgoδ֪
		endif 
		br%mu			=.0001_dp				!Ħ��ϵ��
		br%e1			=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !�ܶ�
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
		br%omgi		=560._dp			!rev/sec  !��ת�ٽϴ�ʱ���������п�����ȫ������Ȧ���Ӷ��γɸ��Ӵ���
		br%mu			=.01_dp				!Ħ��ϵ��
		br%e1			=2.1e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.865e3_dp  !�ܶ�   kg/m_sup3
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
		br%mu			=.001_dp				!Ħ��ϵ��
		br%e1			=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !�ܶ�
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
		br%mu			=.001_dp				!Ħ��ϵ��
		br%e1		=2.06e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.8e3_dp  !�ܶ�
		!case from ������
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
		br%mu			=.03_dp				!Ħ��ϵ��
		br%e1		=2.1e11_dp
		br%nu1		=0.3_dp
		br%rho1		=7.865e3_dp  !�ܶ�   kg/m_sup3
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
			br%omgi		=2000._dp   !�ٶȱ߽�����
			br%omgo		=0._dp   !�ٶȱ߽�����	omgcδ֪
		else
			br%omgi		=7000._dp		!omgc=0,omgoδ֪
		endif
		br%mu				=0.1_dp				!Ħ��ϵ��
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
		br%rho1			=7.8e3_dp  !�ܶ�
	endif

	br%omgi		=rev2rad(br%omgi)		!rad/sec
	br%omgo		=rev2rad(br%omgo)		!rad/sec

	br%e2			=br%e1
	br%nu2		=br%nu1
	br%rho2		=br%rho1  !�ܶ�
	if(br%e1/=br%e2 .or. br%nu1/=br%nu2 .or. br%rho1/=br%rho2) pause 'error in initparas!'
	if(br%m/=br%n) stop
endsubroutine initparas