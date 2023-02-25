!计算对称载荷下滚动体接触角
subroutine conangs_sym(alfi,alfo)
	use	nrtype
	use fsbr
	use neintface
	use intface,only:rad2deg,rad2rev,deg2rad
	implicit none
	real(dp),intent(inout):: alfi,alfo
	integer,parameter:: n=2
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
	stepmx		=1.0000_dp
	method		=2	!1-newton;2-broyden
	global		=4	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=3	!1-mono;2-averged;3-max;4-convex
	eta0			=1_dp
	nmt				=6

	consecmax	=420
	itnlimit	=20000

	x(1)=alfi
	x(2)=alfo

	call neslv(x,fp,gp,qs2eqs_sym,jac24,typx,typy,&
						 analjac,restart,maxtaken,&
						 lshmsd,global,method,termcode,&
						 itncount,itnlimit,itconsec,consecmax,nmt,&
						 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
						 nfnt,njnt)
	alfi=abs(x(1))		!当离心力较大时竟然会出现赋值
	alfo=abs(x(2))
	call slv2eqs_sym(x,fp,.true.)

!	if(br%rcm==1) write(*,*) 'MZK 简化方法计算初值'
!	if(br%rcm==2) write(*,*) '内滚道控制方法计算初值'
!	if(br%rcm==3) write(*,*) '外滚道控制方法计算初值'

!	write(*,*) 'alphai',rad2deg(alfi),alfi
!	write(*,*) 'alphao',rad2deg(alfo),alfo
!	write(*,*) 'resudal',fp
!	write(*,*) 'omgb',br%balno(1)%omgb,rad2rev(br%balno(1)%omgb)
!	write(*,*) 'omgc',br%balno(1)%omgc,rad2rev(br%balno(1)%omgc)
!	write(*,*) 'pitch angle', br%pchang,rad2deg(br%pchang)
!	write(*,*) 'fc', br%balno(1)%fc
!	write(*,*)

	deallocate(x,typx,typy,fp,gp)
endsubroutine conangs_sym

!计算在离心力、无摩擦状态下内、外接触角，作为后续方程组求解的初值
!该子过程的功能:已知轴向力,计算离心力作用下的内外接触角
function qs2eqs_sym(x)
	use fsbr
	use nrtype
	use neintface,only:slv2eqs_sym
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	qs2eqs_sym(size(x(:)))
	real(dp)	:: f(size(x))
	
	call slv2eqs_sym(x,f,.false.)
	qs2eqs_sym=f
endfunction qs2eqs_sym

!计算轴向载荷和离心力作用下，滚动体的内外接触角。
!特点：对称，所有内外接触角均相等。
recursive subroutine slv2eqs_sym(x,f,sv)
	use	nrtype
	use nr
	use fsbr
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)						::f(:)
	logical		:: sv
	real(dp)	::alfi,alfo
	real(dp)	::fc,p_i,p_o
	real(dp)	::dlti,dlto
	real(dp)	::fo,fi
	real(dp)::pchi,pcho,omgb,omgo,e,omgc
	real(dp)	:: rdi,rdo
	real(dp)	::appnd,dy

	alfi	=abs(x(1))
	alfo	=abs(x(2))

	p_i=br%aload/br%z/sin(alfi)
	p_o=br%aload/br%z/sin(alfo)

	call ratint(br%alf,br%icom%appnd,alfi,appnd,dy)
	dlti=appnd*p_i**(2.0_dp/3.0_dp)
	call ratint(br%alf,br%ocom%appnd,alfo,appnd,dy)
	dlto=appnd*p_o**(2.0_dp/3.0_dp)
	rdi=br%r-dlti/2.0_dp
	rdo=br%r-dlto/2.0_dp
	e=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+dlti)*cos(alfi)

	if(br%absolute) then	!一般工况
		!rcm 1:mzk 近似方法，最不好；2：内滚道控制理论,内圈主动：3：外滚道控制理论，内圈主动
		if(br%rcm==1) then
			omgc=((e-rdi*cos(alfi))*br%omgi+(e+rdo*cos(alfo))*br%omgo)/2.0_dp/e
			!omgb wrt. cage
			omgb=(omgc-br%omgi)*(e/br%r-cos(alfi))
			br%balno(1)%omgc=omgc
			br%balno(1)%omgb=omgb
			br%pchang=(alfi+alfo)/2.0_dp		!average
		elseif(br%rcm==2) then
			!内滚道控制!pchi 内滚道控制时滚动体的姿态角
			pchi=atan(e*sin(alfi)/(e*cos(alfi)-rdi))
			omgc=(rdo*rdi*(cos(alfo)*br%omgo-cos(alfi)*br%omgi)+e*(rdi*br%omgo+rdo*br%omgi))/ &
					 (rdo*rdi*(cos(alfo)-cos(alfi))+e*(rdo+rdi))
			omgb=(omgc-br%omgi)*sin(alfi)/(sin(pchi)*cos(alfi)-cos(pchi)*sin(alfi))
			br%balno(1)%omgc=omgc
			br%balno(1)%omgb=omgb
			br%pchang=pchi
		elseif(br%rcm==3) then
			!外滚道控制
			!pcho 外滚道控制时滚动体的姿态角
			!计算公式参考harris著作，注意计算程序中的物理量的彼此代换
			pcho=atan(e*sin(alfo)/(e*cos(alfo)+rdo))
			omgc=(rdo*rdi*(cos(alfo)*br%omgo-cos(alfi)*br%omgi)+e*(rdi*br%omgo+rdo*br%omgi))/ &
					 (rdo*rdi*(cos(alfo)-cos(alfi))+e*(rdo+rdi))
			omgb=(omgc-br%omgo)*sin(alfo)/(sin(pcho)*cos(alfo)-cos(pcho)*sin(alfo))
			!外滚道控制外圈角速度'
			br%balno(1)%omgc=omgc
			br%balno(1)%omgb=omgb
			br%pchang=pcho
		endif
	else	!特殊工况
		!rcm 1:mzk 近似方法，最不好；2：内滚道控制理论,内圈主动：3：外滚道控制理论，内圈主动
		if(br%rcm==1) then
			omgb=br%omgi*(e-rdi*cos(alfi))/rdi
			omgo=omgb*rdo/(e+rdo*cos(alfo))
			br%pchang=(alfi+alfo)/2.0_dp
		elseif(br%rcm==2) then
			!内滚道控制
			!pchi 内滚道控制时滚动体的姿态角
			!计算公式参考harris著作，注意计算程序中的物理量的彼此代换
!			pchi=atan(e*sin(alfi)/(e*cos(alfi)-rdi))
			pchi=alfi-atan(sin(alfi)/(cos(alfi)-e/rdi))
!			omgb=br%omgi*(-e+rdi*cos(alfi))/	&
!									rdi/(cos(pchi)*cos(alfi)+sin(pchi)*sin(alfi))
!			omgb=br%omgi*(e/rdi-cos(alfi))/cos(pchi-alfi)
			omgb=br%omgi*sin(alfi)/sin(pchi-alfi)
!			omgo=omgb/((e+rdo*cos(alfo))/	&
!									rdo/(cos(alfo)*cos(pchi)+sin(pchi)*sin(alfo)))
			omgo=omgb*cos(pchi-alfo)/(e/rdo+cos(alfo))
			!内滚道控制外圈角速度
			br%pchang=pchi
		elseif(br%rcm==3) then
			!外滚道控制
			!pcho 外滚道控制时滚动体的姿态角
			!计算公式参考harris著作，注意计算程序中的物理量的彼此代换
!			pcho=atan(e*sin(alfo)/(e*cos(alfo)+rdo))
			pcho=alfo-atan(sin(alfo)/(cos(alfo)+e/rdo))
!			omgb=br%omgi*(-e+rdi*cos(alfi))/	&
!											rdi/(cos(pcho)*cos(alfi)+sin(pcho)*sin(alfi))
			omgb=br%omgi*(e/rdi-cos(alfi))/cos(pcho-alfi)
!			omgo=omgb/((e+rdo*cos(alfo))/	&
!											rdo/(cos(alfo)*cos(pcho)+sin(pcho)*sin(alfo)))
			omgo=-omgb*sin(pcho-alfo)/sin(alfo)
			!外滚道控制外圈角速度'
			br%pchang=pcho
		endif
		br%omgo=abs(omgo)
		br%balno(1)%omgb=abs(omgb)
		omgc=omgo
	endif
	
	fc=br%m1*e*(omgc**2.0_dp)
	if(sv) then
		br%balno(:)%fc=fc
		br%balno(:)%out%p =p_o 
		br%balno(:)%inn%p =p_i 
		br%balno(:)%out%app =dlto 
		br%balno(:)%inn%app =dlti 
		br%balno(:)%e=e	!滚动体节圆半径
		br%dlta0=(br%ao+dlto)*sin(alfo)+(br%ai+dlti)*sin(alfi)-(0.5_dp*br%d_i+br%ri-0.5_dp*br%d_o+br%ro)*tan(br%alf0)
	endif

	f(1)=p_i*cos(alfi)+fc-p_o*cos(alfo)
	f(2)=((br%ao+dlto)*cos(alfo)+(br%ai+dlti)*cos(alfi))/(br%a*cos(br%alf0))-1.0_dp
endsubroutine slv2eqs_sym
