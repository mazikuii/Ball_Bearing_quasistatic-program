subroutine initbr()
	use fsbr
	use intface,only:pow
	implicit none
	type(patchabk)	 ,target	:: abk
	integer	:: i
	real(dp)	:: dt

	allocate(br%balno(br%z))
	allocate(br%x(br%n))     !用来保存方程组的解
	allocate(br%x0(br%n))     !用来保存方程组的初值，用于比较计算
	do i=1,br%z
		allocate(br%balno(i)%f(br%n))     !保存方程余量
	enddo
	allocate(br%ocom%ke		(br%nsec+1))
	allocate(br%ocom%and	(br%nsec+1))
	allocate(br%ocom%bnd	(br%nsec+1))
	allocate(br%ocom%appnd(br%nsec+1))
	allocate(br%ocom%ei1	(br%nsec+1))
	allocate(br%ocom%ei2	(br%nsec+1))
	allocate(br%ocom%rx		(br%nsec+1))
	allocate(br%ocom%ry		(br%nsec+1))
	allocate(br%ocom%k		(br%nsec+1))
	allocate(br%icom%ke		(br%nsec+1))
	allocate(br%icom%and	(br%nsec+1))
	allocate(br%icom%bnd	(br%nsec+1))
	allocate(br%icom%appnd(br%nsec+1))
	allocate(br%icom%ei1	(br%nsec+1))
	allocate(br%icom%ei2	(br%nsec+1))
	allocate(br%icom%rx		(br%nsec+1))
	allocate(br%icom%ry		(br%nsec+1))
	allocate(br%icom%k		(br%nsec+1))
	allocate(br%alf				(br%nsec+1))
	allocate(br%kn				(br%nsec+1))

	br%rx1		= br%d/2.0_dp
	br%ry1		= br%d/2.0_dp
	br%r			= br%d/2.0_dp
	br%icom%rx2	= br%d_i/2.0_dp
	br%icom%ry2	=-br%ri
	br%ocom%rx2	=-br%d_o/2.0_dp
	br%ocom%ry2	=-br%ro
	br%g1=br%e1/2.0_dp/(1.0_dp+br%nu1)
	br%g2=br%e2/2.0_dp/(1.0_dp+br%nu2)
	if(br%g1/=br%g2) pause "error in inits!"

	br%m1=br%rho1*pi_d*br%d**3.0_dp/6.0_dp
	br%ip=br%rho1*pi_d*br%d**5.0_dp/60.0_dp
	br%dm =0.5_dp*(abs(br%d_i)+abs(br%d_o))
	br%pd =abs(br%d_o)-abs(br%d_i)-2.0_dp*br%d
	br%fo =abs(br%ro/br%d)
	br%fi =abs(br%ri/br%d)
	if(br%hcmthd==3) then !最精确的heathcote近似算方法 mzk
		br%ocom%rc=2.0_dp*br%fo*br%d/(2.0_dp*br%fo+1.0_dp)	!common contact radius of curvature
		br%icom%rc=2.0_dp*br%fi*br%d/(2.0_dp*br%fi+1.0_dp)	!common contact radius of curvature
	endif
	if(br%sym) then			!载荷对称性判断
		br%cyc=1
	else
		br%cyc=br%z
	endif 
	if(br%hertz) then       
		br%np=2.0_dp*pi_d/3.0_dp
	else
		br%np=pi_d/2.0_dp
	endif

	br%a  =abs(br%ro)+abs(br%ri)-abs(br%d)
	br%ai =abs(br%ri)-abs(br%d)/2.0_dp
	br%ao =abs(br%ro)-abs(br%d)/2.0_dp
	br%alf0=acos(1.0_dp-br%pd/(2.0_dp*br%a))

	abk%rx1		=br%rx1
	abk%ry1		=br%ry1
	abk%kp		=br%kp
	abk%e1		=br%e1	!quasiidentical
	abk%e2		=br%e2	!quasiidentical
	abk%nu1		=br%nu1
	abk%nu2		=br%nu2
	
	dt=0.5_dp*pi_d*.95_dp/br%nsec
	abk%ry2		=br%icom%ry2		!不变
	do i=1,br%nsec+1
		!计算内圈接触面椭圆率、半轴
		br%alf(i)=dt*(i-1)
		abk%rx2		=(br%d_i+2.0_dp*br%ri*(1.0_dp-cos(br%alf(i))))/(2.0_dp*cos(br%alf(i)))
		call getabk(abk)
		br%icom%ke(i)			=abk%ke
		br%icom%and(i)		=abk%and
		br%icom%bnd(i)		=abk%bnd
		br%icom%appnd(i)	=abk%appnd
		br%icom%ei1(i)		=abk%ei1
		br%icom%ei2(i)		=abk%ei2
		br%icom%rx(i)			=abk%rx
		br%icom%ry(i)			=abk%ry
		br%icom%k(i)			=pow(1.0_dp/abk%appnd,1.5_dp)
	enddo
	br%icom%xlarger		=abk%xlarger

	abk%ry2					=br%ocom%ry2
	do i=1,br%nsec+1
		!计算外圈接触面椭圆率、半轴
		abk%rx2					=-(br%d_o-2.0_dp*br%ro*(1.0_dp-cos(br%alf(i))))/(2.0_dp*cos(br%alf(i)))
		call getabk(abk)
		br%ocom%ke(i)			=abk%ke
		br%ocom%and(i)		=abk%and
		br%ocom%bnd(i)		=abk%bnd
		br%ocom%appnd(i)	=abk%appnd	
		br%ocom%ei1(i)		=abk%ei1
		br%ocom%ei2(i)		=abk%ei2
		br%ocom%rx(i)			=abk%rx
		br%ocom%ry(i)			=abk%ry
		br%ocom%k(i)			=pow(1.0_dp/abk%appnd,1.5_dp)

		br%kn(i)=1.0_dp/pow(br%icom%appnd(i)+br%ocom%appnd(i),1.5_dp)
	enddo
	br%ocom%xlarger	=abk%xlarger

	br%startang=mod(abs(br%startang),TWOPI_D/br%z)
	br%balno(1)%posang=br%startang
	do i=1,br%z-1
		br%balno(i+1)%posang=br%startang+i*TWOPI_D/br%z
		if(br%balno(i+1)%posang>pi_d) br%balno(i+1)%posang=br%balno(i+1)%posang-TWOPI_D
	enddo
end subroutine initbr

subroutine zeros()
	use fsbr
	implicit none
	integer i
	!global
	pntxo=.false.
	maxit=.false.
	embed=.false.
	mprec=0.0_dp
	sqmprec=0.0_dp

	!fs
	fs%ux=0.0_dp
	fs%uy=0.0_dp
	fs%phx=0.0_dp
	fs%phy=0.0_dp
	fs%uh=0.0_dp
	fs%spole(:,:)=0.0_dp  
	fs%tx(:,:)=0.0_dp
	fs%ty(:,:)=0.0_dp
	fs%sx(:,:)=0.0_dp
	fs%sy(:,:)=0.0_dp
	fs%tfx=0.0_dp
	fs%tfy=0.0_dp


	!br
	br%dlta=0.0_dp
	br%dltr=0.0_dp
	br%theta=0.0_dp
	br%dltast=0.0_dp
	br%dlta0=0.0_dp
	br%s2rinn=0.0_dp
	br%s2rout=0.0_dp
	br%pchang=0.0_dp
	br%pchang_ictl=0.0_dp
	br%pchang_octl=0.0_dp
	br%omgo_ictl=0.0_dp
	br%omgo_octl=0.0_dp
	br%aj=0.0_dp
	br%gama=0.0_dp
	br%tmp1=0.0_dp
	br%tmp2=0.0_dp
	br%iconang_ictl=0.0_dp
	br%oconang_ictl=0.0_dp
	br%iconang_octl=0.0_dp
	br%oconang_octl=0.0_dp
	br%index=0
	br%init=.false.
  br%x(:)=0.0_dp
	br%x0(:)=0.0_dp
	do i=1,br%z
		br%balno(i)%f(:)=0.0_dp
	enddo
	br%icom%c11=0.0_dp
	br%icom%c22=0.0_dp
	br%icom%c23=0.0_dp
	br%icom%c32=0.0_dp
	br%icom%c33=0.0_dp
	br%ocom%c11=0.0_dp
	br%ocom%c22=0.0_dp
	br%ocom%c23=0.0_dp
	br%ocom%c32=0.0_dp
	br%ocom%c33=0.0_dp
	br%balno(:)%e=0.0_dp
	br%balno(:)%e0_ictl=0.0_dp
	br%balno(:)%e0_octl=0.0_dp
	br%balno(:)%est=0.0_dp
	br%balno(:)%mgx=0.0_dp
	br%balno(:)%mgz=0.0_dp
	br%balno(:)%fc=0.0_dp
	br%balno(:)%omgc=0.0_dp
	br%balno(:)%omgb=0.0_dp
	br%balno(:)%omgx=0.0_dp
	br%balno(:)%omgy=0.0_dp
	br%balno(:)%omgz=0.0_dp
	br%balno(:)%inn%ax=0.0_dp
	br%balno(:)%inn%by=0.0_dp
	br%balno(:)%inn%app=0.0_dp
	br%balno(:)%inn%ab=0.0_dp
	br%balno(:)%inn%rdmd=0.0_dp
	br%balno(:)%inn%p=0.0_dp
	br%balno(:)%inn%alfst=0.0_dp !静态接触角contact angle 包括受离心力的情况
	br%balno(:)%inn%alf=0.0_dp !contact angle
	br%balno(:)%inn%xix=0.0_dp
	br%balno(:)%inn%xiy=0.0_dp
	br%balno(:)%inn%psi=0.0_dp	!蠕滑率
	br%balno(:)%inn%ux=0.0_dp
	br%balno(:)%inn%uh=0.0_dp
	br%balno(:)%inn%uy=0.0_dp
	br%balno(:)%inn%phx=0.0_dp
	br%balno(:)%inn%phy=0.0_dp   !无量纲蠕滑量，采用kalker量纲分析方法
	br%balno(:)%inn%tx=0.0_dp
	br%balno(:)%inn%ty=0.0_dp !滚动体受到的集中拖动力
	br%balno(:)%inn%mx=0.0_dp
	br%balno(:)%inn%my=0.0_dp
	br%balno(:)%inn%mz=0.0_dp !滚动体受到的力矩
	br%balno(:)%inn%rd=0.0_dp			!????????
	br%balno(:)%inn%rbd=0.0_dp		!滚动体变形后无量纲半径
	br%balno(:)%inn%spn2rol=0.0_dp  !spin to roll ratio
	br%balno(:)%out%ax=0.0_dp
	br%balno(:)%out%by=0.0_dp
	br%balno(:)%out%app=0.0_dp
	br%balno(:)%out%ab=0.0_dp
	br%balno(:)%out%rdmd=0.0_dp
	br%balno(:)%out%p=0.0_dp
	br%balno(:)%out%alfst=0.0_dp !静态接触角contact angle 包括受离心力的情况
	br%balno(:)%out%alf=0.0_dp !contact angle
	br%balno(:)%out%xix=0.0_dp
	br%balno(:)%out%xiy=0.0_dp
	br%balno(:)%out%psi=0.0_dp	!蠕滑率
	br%balno(:)%out%ux=0.0_dp
	br%balno(:)%out%uh=0.0_dp
	br%balno(:)%out%uy=0.0_dp
	br%balno(:)%out%phx=0.0_dp
	br%balno(:)%out%phy=0.0_dp   !无量纲蠕滑量，采用kalker量纲分析方法
	br%balno(:)%out%tx=0.0_dp
	br%balno(:)%out%ty=0.0_dp !滚动体受到的集中拖动力
	br%balno(:)%out%mx=0.0_dp
	br%balno(:)%out%my=0.0_dp
	br%balno(:)%out%mz=0.0_dp !滚动体受到的力矩
	br%balno(:)%out%rd=0.0_dp			!????????
	br%balno(:)%out%rbd=0.0_dp		!滚动体变形后无量纲半径
	br%balno(:)%out%spn2rol=0.0_dp  !spin to roll ratio
end subroutine zeros