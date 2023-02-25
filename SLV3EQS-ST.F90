!方程组求子过程
!整体方程组，组合载荷作用下
subroutine comstat_st(dlta,dltr,theta)
	use	nrtype
	use fsbr
	use neintface
	use intface,only:rad2deg,rad2rev,deg2rad
	implicit none
	real(dp),intent(inout):: dlta,dltr,theta
	integer,parameter:: n=3
	real(dp),allocatable::x(:),fp(:),gp(:),typx(:),typy(:)
	logical::analjac,restart,maxtaken
	integer::lshmsd,global,method,termcode
	integer::itncount,itnlimit,itconsec,consecmax,nmt,nfnt,njnt
	real(dp)::ftol,xtol,mtol,btol,eta0,stepmx,dlt

	allocate(x(n),typx(n),typy(n),fp(n),gp(n))
	call nedflt(analjac,restart,maxtaken,&
							lshmsd,global,method,termcode,&
							itncount,itnlimit,itconsec,consecmax,nmt,&
							ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
							x,typx,typy)
	ftol			=1.e-16_dp
	xtol			=1.e-16_dp
	mtol			=1.e-16_dp
	btol			=1.e-16_dp
	stepmx		=001.00_dp
	method		=2	!1-newton;2-broyden
	global		=4	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=3	!1-mono;2-averged;3-max;4-convex
	eta0			=.85_dp
	nmt				=6
	dlt				=1.0_dp		!信赖域法的半径
	consecmax	=420
	itnlimit	=20000

	typy	=1.0_dp
	typx	=1.0_dp
	typy			=1.0e-14_dp	!对计算结果和收敛性影响不是很大
	!非常关键,直接关系到精度和收敛性
	typx			=1.0e-8_dp 
!	if(br%rload==0.0_dp .or. br%mload==0.0_dp)	then
!		typx(2)		=1.0_dp	!系数变换,否则会出现异常
!		typx(3)		=1.0_dp
!	endif
	x(1)=dlta
	x(2)=dltr
	x(3)=theta

	call neslv4(x,fp,gp,qs3eqs_st,jac24,typx,typy,&
						 analjac,restart,maxtaken,&
						 lshmsd,global,method,termcode,&
						 itncount,itnlimit,itconsec,consecmax,nmt,&
						 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
						 nfnt,njnt)
	dlta	=x(1)		!当离心力较大时竟然会出现赋值
	dltr	=x(2)
	theta	=x(3)
	deallocate(x,typx,typy,fp,gp)
endsubroutine comstat_st

!该函数为slv3eqs的封装函数,便于对slv3eqs的调用形式进行随意修改.
function qs3eqs_st(x)
	use fsbr
	use nrtype
	use neintface,only:slv3eqs_st
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	qs3eqs_st(size(x(:)))
	real(dp)::f(size(x))

	call slv3eqs_st(x,f,.false.)
	qs3eqs_st(:)=f(:)
endfunction qs3eqs_st

!实施子过程,不考虑惯性效应,静平衡
subroutine slv3eqs_st(x,f,sv)
	use fsbr
	use nrtype
	use nr
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	f(:)
	logical	::	sv

	real(dp)	:: qn,salf,calf
	integer		:: i
	real(dp)	:: dltn,posang,a,ri,ro
	real(dp)	:: cospa,salf0,calf0
	real(dp)  :: dlta,dltr,theta
	real(dp)	:: temp,sn,cs,e,gama
	real(dp)	:: kn,dy,k

	dlta	=x(1)
	dltr	=x(2)
	theta	=x(3)

	a=br%a
	salf0=sin(br%alf0)
	calf0=cos(br%alf0)
	ri=br%d_i/2.0_dp+br%ri	!内滚道中心轨迹圆半径
	ro=br%d_o/2.0_dp-br%ro	!外滚道中心轨迹圆半径

	f(:)=0.0_dp
	do i=1,br%z
		posang=br%balno(i)%posang
		cospa=cos(posang)
		sn=a*salf0+dlta+ri*sin(theta)*cospa
		cs=a*calf0+dltr*cospa-ri*(1.0_dp-cos(theta))*cospa
		temp=sqrt(sn**2+cs**2)
		dltn=temp-a
		if(dltn<=0.0_dp) then	!滚动体与滚道脱离接触,按MZK方法计算.
			dltn=0.0_dp	
			qn=0.0_dp
		else !没有脱离,按JONES方法计算
			salf=sn/temp
			calf=cs/temp
			call ratint(br%alf,br%kn,asin(salf),kn,dy)
			qn=kn*(dltn**1.5_dp)
		endif
		if(sv) then
			if(dltn==0.0_dp) then	!滚动体与滚道脱离接触,按MZK方法计算.
				gama=acos(sn/temp)  !gama  借用变量
				sn=acos((temp**2+br%ai**2-br%ao**2)/(2.0_dp*temp*br%ai))  !theta1
				cs=acos((temp**2+br%ao**2-br%ai**2)/(2.0_dp*temp*br%ao))	!theta2
				br%balno(i)%inn%alf=PIO2_D-gama+sn
				br%balno(i)%out%alf=PIO2_D-gama-cs
				br%balno(i)%e=ro+br%ao*sin(gama+cs)
			else
				br%balno(i)%inn%alf=asin(salf)
				br%balno(i)%out%alf=asin(salf)
				call ratint(br%alf,br%ocom%k,asin(salf),k,dy)
				br%balno(i)%e=ro+(br%ao+(qn/k)**(2.0_dp/3.0_dp))*calf
			endif
			br%balno(i)%inn%p=qn
			br%balno(i)%out%p=qn
			call ratint(br%alf,br%icom%k,asin(salf),k,dy)
			br%balno(i)%inn%app=(qn/k)**(2.0_dp/3.0_dp)
			call ratint(br%alf,br%ocom%k,asin(salf),k,dy)
			br%balno(i)%out%app=(qn/k)**(2.0_dp/3.0_dp)
		endif
		f(1)=f(1)+qn*salf
		f(2)=f(2)+qn*calf*cospa
		f(3)=f(3)+qn*ri*(salf*cos(theta)-calf*sin(theta))*cospa		!qs3eqs(3)=qs3eqs(3)+qn*ri*sin(asin(salf)-theta)*cospa
	enddo
	f(1)=br%aload-f(1)
	f(2)=br%rload-f(2)
	f(3)=br%mload-f(3)
endsubroutine slv3eqs_st