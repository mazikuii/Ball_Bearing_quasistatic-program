subroutine slvcmbd()
	!������Ħ����̬�Ӵ���
	use fsbr
	use nr
	use neintface,only:slv3eqs_st,slv3eqs_inert,slv3eqs_creep
	USE INTFACE,ONLY:RAD2DEG
	use shell32
	implicit none
	real(dp)							:: contang,load
	real(dp)							:: dlti,dlto
	real(dp)							:: alfi,alfo
	real(dp)		:: appnd,dy
	real(dp)	:: x(3),f(3)
	integer	i
	character(len=1)			:: tab=char(9)
	
	contang	=br%balno(1)%inn%alfst
	load		=br%balno(1)%inn%p
	!���ܵ����������غ�,���һ���������м���
	!��������λ�Ƶĳ�ֵ
	call ratint(br%alf,br%icom%appnd,contang,appnd,dy)
	dlti=appnd*load**(2.0_dp/3.0_dp)
	call ratint(br%alf,br%ocom%appnd,contang,appnd,dy)
	dlto=appnd*load**(2.0_dp/3.0_dp)
	br%dlta	=(br%a+dlti+dlto)*sin(contang)-br%a*sin(br%alf0)
	br%dltr	=1.0e-8_dp		!��ֵС�˺�
	br%theta=1.0e-8_dp

	!�����ڸ��������غ�����������λ��,����λ�ƺͽ�λ��
	!��������������������ܼ򵥣�ͬһ�����������Ӵ�����ȡ�
	call comstat_st(br%dlta,br%dltr,br%theta)
	x(1)=br%dlta
	x(2)=br%dltr
	x(3)=br%theta
	call slv3eqs_st(x,f,.true.)	!�����м����
	write(*,*) '*************  static analysis succeeded! ************'
	write(*,*) x
	write(*,*)
	write(*,*) f
	write(*,*) 'max moment:',br%aload*br%dm*.5_dp
	write(*,*) '******************************************************'
	if(br%inert>=1) then
		!������������,����Ӵ��ǲ���ȣ�����ʱ��Ҫ���ǹ�����Ĺ�ת���ٶ�
		!���ò����ǽ��ٶȵ�λ������Ϊ��ֵ
		br%balno(:)%inn%app=br%balno(:)%inn%app*1.000004_dp		!ϵ���ܹؼ���Ӱ��������
		br%balno(:)%out%app=br%balno(:)%out%app*1.000004_dp
		call comstat_inert(br%dlta,br%dltr,br%theta)
		x(1)=br%dlta
		x(2)=br%dltr
		x(3)=br%theta
		call slv3eqs_inert(x,f,.true.)	!�����м����
		write(*,*) '*******  quasi-static analysis succeeded!   **********'
		write(*,*) x
		write(*,*)
		write(*,*) f
		write(*,*) '******************************************************'
	endif
	if(br%cmdcreep) then
		call comstat_creep(br%dlta,br%dltr,br%theta)
		x(1)=br%dlta
		x(2)=br%dltr
		x(3)=br%theta
		call slv3eqs_creep(x,f,.true.)
	endif

	br%dlta=x(1)
	br%dltr=x(2)
	br%theta=x(3)
	br%f3(:)=f(:)
endsubroutine slvcmbd