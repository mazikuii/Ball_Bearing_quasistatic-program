subroutine slvcmbd()
	!计算无摩擦静态接触角
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
	!若受到联合作用载荷,则进一步进行下列计算
	!计算轴向位移的初值
	call ratint(br%alf,br%icom%appnd,contang,appnd,dy)
	dlti=appnd*load**(2.0_dp/3.0_dp)
	call ratint(br%alf,br%ocom%appnd,contang,appnd,dy)
	dlto=appnd*load**(2.0_dp/3.0_dp)
	br%dlta	=(br%a+dlti+dlto)*sin(contang)-br%a*sin(br%alf0)
	br%dltr	=1.0e-8_dp		!初值小了好
	br%theta=1.0e-8_dp

	!计算在给定联合载荷作用下轴向位移,径向位移和角位移
	!不考虑离心力，因而求解很简单，同一滚动体的内外接触角相等。
	call comstat_st(br%dlta,br%dltr,br%theta)
	x(1)=br%dlta
	x(2)=br%dltr
	x(3)=br%theta
	call slv3eqs_st(x,f,.true.)	!保存中间变量
	write(*,*) '*************  static analysis succeeded! ************'
	write(*,*) x
	write(*,*)
	write(*,*) f
	write(*,*) 'max moment:',br%aload*br%dm*.5_dp
	write(*,*) '******************************************************'
	if(br%inert>=1) then
		!若考虑离心力,内外接触角不相等，计算时需要考虑滚动体的公转角速度
		!采用不考虑角速度的位移量作为初值
		br%balno(:)%inn%app=br%balno(:)%inn%app*1.000004_dp		!系数很关键，影响收敛性
		br%balno(:)%out%app=br%balno(:)%out%app*1.000004_dp
		call comstat_inert(br%dlta,br%dltr,br%theta)
		x(1)=br%dlta
		x(2)=br%dltr
		x(3)=br%theta
		call slv3eqs_inert(x,f,.true.)	!保存中间变量
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