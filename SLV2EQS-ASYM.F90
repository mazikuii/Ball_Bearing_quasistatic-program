!���㲻�Գ��غ��¹�����Ӵ���
subroutine conangs_asym(alfi,alfo)
	use	nrtype
	use fsbr
	use neintface
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
	ftol			=1.e-12_dp
	xtol			=1.e-12_dp
	mtol			=1.e-12_dp
	btol			=1.e-12_dp
	stepmx		=1.0000_dp
	method		=2	!1-newton;2-broyden
	global		=4	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=4	!1-mono;2-averged;3-max;4-convex
	eta0			=1_dp
	nmt				=6

	consecmax	=420
	itnlimit	=20000

!	typy			=1.0e-12_dp
	typx			=1.0e-10_dp
	x(1)=alfi
	x(2)=alfo

	call neslv4(x,fp,gp,qs2eqs_asym,jac24,typx,typy,&
						 analjac,restart,maxtaken,&
						 lshmsd,global,method,termcode,&
						 itncount,itnlimit,itconsec,consecmax,nmt,&
						 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
						 nfnt,njnt)
	alfi=abs(x(1))		!���������ϴ�ʱ��Ȼ����ָ�ֵ
	alfo=abs(x(2))

	deallocate(x,typx,typy,fp,gp)
endsubroutine conangs_asym

!����������������Ħ��״̬���ڡ���Ӵ��ǣ���Ϊ�������������ĳ�ֵ
!���ӹ��̵Ĺ���:��֪������,���������������µ�����Ӵ���
function qs2eqs_asym(x)
	use fsbr
	use nrtype
	use neintface,only:slv2eqs_asym
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	qs2eqs_asym(size(x(:)))
	real(dp)	:: f(size(x))
	
	call slv2eqs_asym(x,f,.false.)
	qs2eqs_asym=f
endfunction qs2eqs_asym

!���������غ�������ĳ�������������Ӵ��ǡ�
!���ڹ�����������λ����֪������·��ƽӴ�����������ٶȱ߽�����
!�����ض������µĽӴ��ǣ���ʱ����Ӵ��ǻ�����ȡ�
!�����õ�����Ӵ��غ�
subroutine slv2eqs_asym(x,f,sv)
	use	nrtype
	use fsbr
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)						::	f(:)
	logical		::sv
	real(dp)	::alfi,alfo
	real(dp)	::fc,p_i,p_o
	real(dp)	::dlti,dlto
	real(dp)	::fo,fi
	real(dp)::pchi,pcho,omgb,omgo,e,omgc
	real(dp)	:: rdi,rdo

	real(dp)	::aj,gama,bt1,bt2,a,ai,ao
	real(dp)	::ki,ko
	real(dp)	::tmp

	a=br%a
	ai=br%ai
	ao=br%ao
	ki=br%tmp1			!�ն�ϵ����Ϊ�Ӵ��ǵĺ���
	ko=br%tmp2			!�ն�ϵ����Ϊ�Ӵ��ǵĺ���
	aj=br%aj				!���ڲ������ݣ����õĴ��ݷ�ʽ
	gama=br%gama		!���ڲ������ݣ����õĴ��ݷ�ʽ
	if(br%intmthd==1) then
		alfi	=abs(x(1))
		alfo	=abs(x(2))
		dlti= aj*cos(alfo+gama)/sin(alfi-alfo)-ai
		dlto=-aj*cos(alfi+gama)/sin(alfi-alfo)-ao
		if(dlti<0.0_dp .or. dlto<0.0_dp) then !��ѧλ���ļ���
!			dlti=0.0_dp
!			dlto=0.0_dp
		else	!��������

		endif
		rdi=br%r-dlti/2.0_dp
		rdo=br%r-dlto/2.0_dp
	elseif(br%intmthd==2) then
		dlti=x(1)
		dlto=x(2)
		rdi=br%r-dlti/2.0_dp
		rdo=br%r-dlto/2.0_dp
		tmp=(aj*aj+(ai+dlti)**2-(ao+dlto)**2)/(2.0_dp*aj*(ai+dlti))
		if(tmp>1.0_dp) tmp=.9999999_dp
		bt1=acos(tmp)
		tmp=(aj*aj+(ao+dlto)**2-(ai+dlti)**2)/(2.0_dp*aj*(ao+dlto))
		if(tmp>1.0_dp) tmp=.9999999_dp
		bt2=acos(tmp)
		alfi=PIO2_D-gama+bt1
		alfo=PIO2_D-gama-bt2
	endif
	p_i=ki*dlti**1.5_dp
	p_o=ko*dlto**1.5_dp
	e=abs(br%ocom%rx2)-abs(br%ocom%ry2)+(ao+dlto)*cos(alfo)
	!rcm 1:mzk ���Ʒ�������ã�2���ڹ�����������,��Ȧ������3��������������ۣ���Ȧ����
	if(br%rcm==1) then
		omgc=((e-rdi*cos(alfi))*br%omgi+(e+rdo*cos(alfo))*br%omgo)/2.0_dp/e
		br%balno(1)%omgc=omgc
	elseif(br%rcm==2) then
		!�ڹ�������!pchi �ڹ�������ʱ���������̬��
		pchi=atan(e*sin(alfi)/(e*cos(alfi)-rdi))
		omgc=(rdo*rdi*(cos(alfo)*br%omgo-cos(alfi)*br%omgi)+e*(rdi*br%omgo+rdo*br%omgi))/ &
				 (rdo*rdi*(cos(alfo)-cos(alfi))+e*(rdo+rdi))
		br%balno(1)%omgc=omgc
	elseif(br%rcm==3) then
		!���������
		!pcho ���������ʱ���������̬��
		!���㹫ʽ�ο�harris������ע���������е��������ı˴˴���
		pcho=atan(e*sin(alfo)/(e*cos(alfo)+rdo))
		omgc=(rdo*rdi*(cos(alfo)*br%omgo-cos(alfi)*br%omgi)+e*(rdi*br%omgo+rdo*br%omgi))/ &
				 (rdo*rdi*(cos(alfo)-cos(alfi))+e*(rdo+rdi))
		br%balno(1)%omgc=omgc
	endif

	if(sv) then
		br%aj		=p_i    !��ʱ���ã�ע��
		br%gama	=p_o
		br%tmp1=alfi
		br%tmp2=alfo
	endif
	fc=br%m1*e*(omgc**2.0_dp)
	br%balno(1)%fc=fc

	f(1)=p_i*cos(alfi)+fc-p_o*cos(alfo)
	f(2)=p_i*sin(alfi)-p_o*sin(alfo)
endsubroutine slv2eqs_asym