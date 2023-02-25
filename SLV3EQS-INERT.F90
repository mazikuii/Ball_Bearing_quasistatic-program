!方程组求子过程
!整体方程组，组合载荷作用下
subroutine comstat_inert(dlta,dltr,theta)
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
	stepmx		=1.00_dp
	method		=2	!1-newton;2-broyden
	global		=4	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=2	!1-mono;2-averged;3-max;4-convex
	eta0			=.85_dp
	nmt				=6
	consecmax	=420
	itnlimit	=20000
	dlt				=1.000_dp

!	typy			=1.0e-12_dp	
	!非常关键,直接关系到精度和收敛性
	typx			=1.0e-12_dp 
	x(1)=dlta
	x(2)=dltr
	x(3)=theta

	call neslv(x,fp,gp,qs3eqs_inert,jac24,typx,typy,&
						 analjac,restart,maxtaken,&
						 lshmsd,global,method,termcode,&
						 itncount,itnlimit,itconsec,consecmax,nmt,&
						 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
						 nfnt,njnt)
	dlta	=x(1)		!当离心力较大时竟然会出现赋值
	dltr	=x(2)
	theta	=x(3)
	deallocate(x,typx,typy,fp,gp)
endsubroutine comstat_inert

!该函数为slv3eqs的封装函数,便于对slv3eqs的调用形式进行随意修改.
function qs3eqs_inert(x)
	use fsbr
	use nrtype
	use neintface,only:slv3eqs_inert
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	qs3eqs_inert(size(x(:)))
	real(dp)::f(size(x))

	call slv3eqs_inert(x,f,.false.)
	qs3eqs_inert(:)=f(:)
endfunction qs3eqs_inert

!实施子过程,考虑离心力作用,忽略摩擦力.
subroutine slv3eqs_inert(x,f,sv)
	use fsbr
	use nrtype
	use nr,only:ratint
	use neintface,only:conangs,slv2eqs_asym
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	f(:)
	logical	::	sv

	real(dp)	:: qi,salf,calf
	integer		:: i
	real(dp)	:: posang,a,ri,ro
	real(dp)	:: cospa
	real(dp)  :: dlta,dltr,theta
	real(dp)	:: sn,cs,e
	real(dp)	:: alf(2),fvec(2)
	integer :: j
	real(dp)	:: k,dy
	logical	:: harris=.false.

	dlta				=x(1)
	dltr				=x(2)
	theta				=x(3)

	a=br%a
	ri=br%d_i/2.0_dp+br%ri	!内滚道中心轨迹圆半径
	ro=br%d_o/2.0_dp-br%ro	!外滚道中心轨迹圆半径
	f(:)=0.0_dp		!词句很关键,千万不能漏掉
	do i=1,br%z
		posang=br%balno(i)%posang
		cospa=cos(posang)
		sn=a*sin(br%alf0)+dlta+ri*sin(theta)*cospa
		cs=a*cos(br%alf0)+dltr*cospa-ri*(1.0_dp-cos(theta))*cospa
		br%aj=sqrt(sn**2+cs**2)
		br%gama=atan(cs/sn)
		if(br%aj<a .and. harris) then
			alf(:)=0.0_dp
			qi=0.0_dp
			br%aj=0.0_dp
			br%gama=0.0_dp
		else
			if(br%intmthd==1) then
				alf(1)=br%balno(i)%inn%alf
				alf(2)=br%balno(i)%out%alf
				if(alf(1)==alf(2)) then
					alf(1)=alf(1)*1.0001_dp
					alf(2)=alf(2)/1.0001_dp
				endif
				call ratint(br%alf,br%icom%k,alf(1),k,dy)
				br%tmp1=k
				call ratint(br%alf,br%ocom%k,alf(2),k,dy)
				br%tmp2=k
!				call conangs_asym(alf(1),alf(2))  !NESLV方法,效果不好
				call HYBRD_asym(alf)	!比NESLV效果好
				call slv2eqs_asym(alf,fvec,.true.)
				br%balno(i)%inn%alf=alf(1)
				br%balno(i)%out%alf=alf(2)
				qi=br%aj	!临时借用,传递内接触力
			elseif(br%intmthd==2) then
				alf(1)=br%balno(i)%inn%app
				alf(2)=br%balno(i)%out%app
				if(alf(1)<=0.0_dp) alf(1)=1.0e-5_dp	!系数很关键，影响收敛性
				if(alf(2)<=0.0_dp) alf(2)=1.0e-5_dp
				call ratint(br%alf,br%icom%k,alf(1),k,dy)
				br%tmp1=k
				call ratint(br%alf,br%ocom%k,alf(2),k,dy)
				br%tmp2=k
!				call conangs_asym(alf(1),alf(2))  !NESLV方法,效果不好
				call HYBRD_asym(alf)	!比NESLV效果好
				call slv2eqs_asym(alf,fvec,.true.)
				br%balno(i)%inn%app=alf(1)	!保存delta,备下次循环用
				br%balno(i)%out%app=alf(2)	!保存delta,备下次循环用
				qi=br%aj	!临时借用,传递内接触力
				alf(1)=br%tmp1	!临时借用,传递alfi
				alf(2)=br%tmp2	!临时借用,传递alfo
			endif
		endif
		salf=sin(alf(1))
		calf=cos(alf(1))
		if(sv) then
			br%balno(i)%inn%alf=alf(1)
			br%balno(i)%out%alf=alf(2)
			br%balno(i)%inn%p=br%aj
			br%balno(i)%out%p=br%gama
			br%balno(i)%e=ro+(br%ao+br%balno(i)%out%app)*alf(2)
		endif
		f(1)=f(1)+qi*salf
		f(2)=f(2)+qi*calf*cospa
		f(3)=f(3)+qi*ri*(salf*cos(theta)-calf*sin(theta))*cospa		
	enddo
	f(1)=br%aload-f(1)
	f(2)=br%rload-f(2)
	f(3)=br%mload-f(3)
endsubroutine slv3eqs_inert