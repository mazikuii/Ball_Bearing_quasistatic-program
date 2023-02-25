!方程组求子过程
!整体方程组，组合载荷作用下
subroutine comstat_creep(dlta,dltr,theta)
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
	INTEGER I
	allocate(x(n),typx(n),typy(n),fp(n),gp(n))
	call nedflt(analjac,restart,maxtaken,&
							lshmsd,global,method,termcode,&
							itncount,itnlimit,itconsec,consecmax,nmt,&
							ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
							x,typx,typy)
	ftol			=1.e-12_dp
	xtol			=1.e-12_dp
	mtol			=1.e-12_dp
	btol			=1.e-12_dp
	stepmx		=1._dp
	method		=2	!1-newton;2-broyden
	!1的效果反而最好，2的求解有问题，皆不出来
	!3和4都在当处于极限大力矩状态时出现发散的问题，但当轴承受载的不对称性较小时的求解
	!精度还是很高的．
	global		=1	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=3	!1-mono;2-averged;3-max;4-convex
	eta0			=1._dp
	nmt				=6
	dlt				=100.0_dp		!信赖域法的半径
	consecmax	=420
	itnlimit	=20000

	typy	=1.0_dp
	typx	=1.0_dp
!	typy			=1.0e-14_dp	!对计算结果和收敛性影响不是很大
	!非常关键,直接关系到精度和收敛性
	typx			=1.0e-10_dp 

	x(1)=dlta
	x(2)=dltr
	x(3)=theta
	br%init	=.true.

	do i=1,2
	if(i==2) then
		global		=4	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	endif
	call neslv4(x,fp,gp,qs3eqs_creep,jac24,typx,typy,&
						 analjac,restart,maxtaken,&
						 lshmsd,global,method,termcode,&
						 itncount,itnlimit,itconsec,consecmax,nmt,&
						 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
						 nfnt,njnt)
	enddo
	dlta	=x(1)		!当离心力较大时竟然会出现赋值
	dltr	=x(2)
	theta	=x(3)
	deallocate(x,typx,typy,fp,gp)
endsubroutine comstat_creep

!该函数为slv3eqs的封装函数,便于对slv3eqs的调用形式进行随意修改.
function qs3eqs_creep(x)
	use fsbr
	use nrtype
	use neintface,only:slv3eqs_creep
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	qs3eqs_creep(size(x(:)))
	real(dp)::f(size(x))

	call slv3eqs_creep(x,f,.false.)
	qs3eqs_creep(:)=f(:)
write(*,*) 
write(*,*) '...................'
write(*,*) x
write(*,*) f
write(*,*) '...................'
endfunction qs3eqs_creep

!实施子过程,不考虑惯性效应,静平衡
subroutine slv3eqs_creep(x,f,sv)
	use fsbr
	use nr,only:ratint
	use nrtype
	use neintface,only:conangs,slv2eqs_asym,slv6eqs_creep,slv6eqs_idx
	use intface,only:pow
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)	::	f(:)
	logical		::	sv

	real(dp)	::	qi,salf,calf
	integer		::	i
	real(dp)	::	posang,a,ri,ro
	real(dp)	::	cospa
	real(dp)  ::	dlta,dltr,theta
	real(dp)	::	sn,cs,e
	real(dp)	::	fvec(6)
	integer		::	j

	real(dp)	::	est,dlti,dlto
	real(dp)	::	pchi,pcho,omgb,omgo,omgc
	real(dp)	::	pchang
	real(dp)	::	alfi,alfo
	real(dp)	::	rdi,rdo,tmp,bt1,bt2
	real(dp)	::	tx,ty,mz
	real(dp)	::	kiko,ki,ko,dy

	dlta				=x(1)
	dltr				=x(2)
	theta				=x(3)

	br%dlta=dlta		!保存当前位移量，传递到滚动体求解程序段
	br%dltr=dltr
	br%theta=theta
	a=br%a
	ri=br%d_i/2.0_dp+br%ri	!内滚道中心轨迹圆半径
	ro=br%d_o/2.0_dp-br%ro	!外滚道中心轨迹圆半径
	f(:)=0.0_dp		!词句很关键,千万不能漏掉
	do i=1,br%z
		br%index=i
		posang=br%balno(i)%posang
		cospa=cos(posang)
		sn=a*sin(br%alf0)+dlta+ri*sin(theta)*cospa
		cs=a*cos(br%alf0)+dltr*cospa-ri*(1.0_dp-cos(theta))*cospa
		br%aj=sqrt(sn**2+cs**2)
		br%gama=atan(cs/sn)

		if(br%init) then
			est=br%balno(i)%e
			rdi=br%r-br%balno(i)%inn%app*0.5_dp
			rdo=br%r-br%balno(i)%out%app*0.5_dp
			!1:mzk 近似方法，最不好；2：内滚道控制理论,内圈主动：3：外滚道控制理论，内圈主动
			alfi=br%balno(i)%inn%alf
			alfo=br%balno(i)%out%alf
			call ratint(br%alf,br%icom%k	,alfi,ki		,dy)
			call ratint(br%alf,br%ocom%k	,alfo,ko		,dy)
			kiko=(ki/ko)**(2.0_dp/3.0_dp)
			br%tmp1=(br%aj-br%a)/(1.0_dp+kiko)
			br%tmp2=br%tmp1*kiko
			if(br%rcm==1) then
				omgc=((est-rdi*cos((alfi+alfo)*.5_dp))*br%omgi+(est+rdo*cos((alfi+alfo)*.5_dp))*br%omgo)/2.0_dp/est
				omgb=(omgc-br%omgi)*(est/rdi-cos((alfi+alfo)*.5_dp))
				pchang=(alfi+alfo)*.5_dp
			elseif(br%rcm==2) then
				!内滚道控制!pchi 内滚道控制时滚动体的姿态角
				pchi=atan(est*sin(alfi)/(est*cos(alfi)-rdi))
				omgc=(est*(br%omgo+br%omgi)+rdi*(br%omgo-br%omgi)*cos(alfi))/2.0_dp/est
				omgb=(omgc-br%omgi)*sin(alfi)/(sin(pchi)*cos(alfi)-cos(pchi)*sin(alfi))
				pchang=pchi
			elseif(br%rcm==3) then
				!外滚道控制
				!pcho 外滚道控制时滚动体的姿态角
				pcho=atan(est*sin(alfo)/(est*cos(alfo)+rdo))
				omgc=(est*(br%omgo+br%omgi)+rdi*(br%omgo-br%omgi)*cos(alfi))/2.0_dp/est
				omgb=(omgc-br%omgo)*sin(alfo)/(sin(pcho)*cos(alfo)-cos(pcho)*sin(alfo))
				pchang=pcho
			endif
			br%x(1)=omgc
			br%x(2)=0.001_dp
			br%x(3)=omgb*cos(pchang)
			br%x(4)=omgb*sin(pchang)
			if(br%crpmthd==1) then
				br%x(5)=br%balno(i)%inn%alf	*1.01_dp
				br%x(6)=br%balno(i)%out%alf	/1.01_dp
			elseif(br%crpmthd==2) then
				br%x(5)=br%balno(i)%inn%app
				br%x(6)=br%balno(i)%out%app
				br%x(5)=br%x(5)	!*1.02_dp		!MZK需要
				br%x(6)=br%x(6)	!*1.02_dp
			endif
		else
			br%x(1)=br%balno(i)%omgc
			br%x(2)=br%balno(i)%omgx
			br%x(3)=br%balno(i)%omgy	
			br%x(4)=br%balno(i)%omgz
			if(br%crpmthd==1) then
				br%x(5)=br%balno(i)%inn%alf
				br%x(6)=br%balno(i)%out%alf
			elseif(br%crpmthd==2) then
				br%x(5)=br%balno(i)%inn%app
				br%x(6)=br%balno(i)%out%app
			endif
		endif

		call slv6eqs_creep(br%x)	!计算第I个滚动体在给定几何条件下的平衡状态
!		call wrttm(i)  !输出拖动力、拖动力矩
		call slv6eqs_idx(br%x,fvec,.true.,i,br%crpmthd)
write(*,*) br%x(1),br%x(3),br%x(4)
!write(*,*) i,sqrt(dot_product(fvec,fvec))

		qi=br%balno(i)%inn%p
		alfi=br%balno(i)%inn%alf
		alfo=br%balno(i)%out%alf
		salf=sin(alfi)
		calf=cos(alfi)

		f(1)=f(1)+qi*salf
		f(2)=f(2)+qi*calf*cospa
		f(3)=f(3)+qi*ri*(salf*cos(theta)-calf*sin(theta))*cospa		
	enddo
	f(1)=br%aload-f(1)
	f(2)=br%rload-f(2)
	f(3)=br%mload-f(3)
	if(br%init) br%init=.false.
endsubroutine slv3eqs_creep