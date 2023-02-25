subroutine slv6eqs_idx(x,f,sv,idx,mthd)
	use fsbr
	use nrtype
	use nr,only:ratint
	USE INTFACE,ONLY:RAD2DEG
	implicit none
	real(dp), dimension(:)	:: x
	real(dp), dimension(:)	:: f
	logical		::sv
	integer		::idx
	integer		::mthd	!1:接触角->变形量;2:变形量->接触角
	character(len=1)	:: tab=char(9)
	real(dp)	::alfi,alfo
	real(dp)	::omgi,omgo,omgc
	real(dp)	::omgx,omgy,omgz
	real(dp)	::xixo,xixi,xiyo,xiyi,psio,psii
	real(dp)	::num,demo !局 部变量，分子、分母
	real(dp)	::cai,cao,sai,sao
	real(dp)	::txi,tyi,mxi,myi,mzi
	real(dp)	::txo,tyo,mxo,myo,mzo
	real(dp)	::mgz,mgx,fc,p_i,p_o
	real(dp)	::dlti,dlto
	real(dp)	::ro,ri
	real(dp)	::e
	logical		::outer
	integer		::i
	real(dp)	::tmp,k,dy
	real(dp)	::inn21,out21
	real(dp)	::aj,gama,bt1,bt2,posang

	aj=br%aj
	gama=br%gama
	posang=br%balno(idx)%posang
	inn21=br%tmp1	!临时借用
	out21=br%tmp2	!临时借用
	if(br%absolute) then
		omgi=br%omgi		!known
		omgo=br%omgo		!known
		omgc=x(1)				!unknown
	else
		omgi=br%omgi		!known
		omgo=x(1)
		omgc=omgo
	endif
	omgx	=x(2)    
	omgy	=x(3)
	omgz	=x(4)

	if(mthd==1) then
		alfi=x(5)
		alfo=x(6)
		dlti= aj*cos(alfo+gama)/sin(alfi-alfo)-br%ai
		dlto=-aj*cos(alfi+gama)/sin(alfi-alfo)-br%ao
		dlti=abs(dlti)		!!!!!!!!!!!没有道理的道理
		dlto=abs(dlto)		!!!!!!!!!!!没有道理的道理
	elseif(mthd==2) then
		x(5)=abs(x(5))				!!!!!!!!!!!没有从根本上解决该问题
		x(6)=abs(x(6))				!!!!!!!!!!!没有从根本上解决该问题

		dlti	=x(5)
		dlto	=x(6)
!		if(dlti<inn21.or. dlto<out21) then
!			dlti=inn21
!			dlto=out21
!			write(*,*) dlti,dlto, 'wrong'
!		endif	

		tmp=(aj*aj+(br%ai+dlti)**2-(br%ao+dlto)**2)/(2.0_dp*aj*(br%ai+dlti))
!		if(tmp>1.0_dp) tmp=.9999999999_dp	!MZK方法需要该语句
		if(tmp>1.0_dp .or. tmp<.56_dp) tmp=.99999999_dp	!MZK方法需要该语句
		bt1=acos(tmp)
		tmp=(aj*aj+(br%ao+dlto)**2-(br%ai+dlti)**2)/(2.0_dp*aj*(br%ao+dlto))
!		if(tmp>1.0_dp) tmp=.9999999999_dp
		if(tmp>1.0_dp .or. tmp<.56_dp) tmp=.99999999_dp	!MZK方法需要该语句
		bt2=acos(tmp)
		alfi=pio2_d-gama+bt1
		alfo=pio2_d-gama-bt2
	endif

	cai=cos(alfi)
	sai=sin(alfi)
	cao=cos(alfo)
	sao=sin(alfo)

	call ratint(br%alf,br%icom%k	,alfi,k		,dy)
	p_i=k*dlti**1.5_dp
	call ratint(br%alf,br%ocom%k	,alfo,k		,dy)
	p_o=k*dlto**1.5_dp
	br%balno(idx)%inn%p =p_i
	br%balno(idx)%out%p =p_o
	call getpchaxis_idx(idx,alfi,alfo)

	e=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+dlti)*cai
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!creepage and spin
	ri=br%r-dlti/2.0_dp
	ro=br%r-dlto/2.0_dp
	if(br%absolute) then  !绝对速度
		tmp=omgy*cao+omgz*sao
		num	= e*(omgc-omgo)+ro*(tmp-omgo*cao+omgc*cao)
		demo=(e*(omgo-omgc)+ro*(tmp+omgo*cao-omgc*cao))/2.0_dp		
		xixo=-num/demo
		num	=ro*omgx   
		xiyo=-num/demo
		num	=(omgy*sao-omgz*cao-omgo*sao +omgc*sao)
		psio=-num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(idx)%out%spn2rol=num/tmp 	!详见学位论文
			else
		    br%balno(idx)%out%spn2rol=num/(tmp+(omgc-omgo)*cao)	!详见学位论文
			endif
		endif
!		if(sv) br%balno(:)%out%spn2rol=(omgb*sin(alphao-alpha)+(omgc-omgo)*sin(alphao))/(omgb*cos(alphao-alpha)+(omgc-omgo)*cos(alphao))	!化简公式
		if(br%hcmthd==1) then
			out21=0.0_dp
		else
			out21=omgi*cai/tmp
		endif

		tmp=omgy*cai+omgz*sai
		num = (omgi-omgc)*e+ri*(tmp-omgi*cai+omgc*cai)
		demo=((omgi-omgc)*e-ri*(tmp+omgi*cai-omgc*cai))/2.0_dp		
		xixi=num/demo
		num =ri*omgx  
		xiyi=-num/demo
		num	=(omgz*cai-omgy*sai+omgi*sai -omgc*sai)
		psii=num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(idx)%inn%spn2rol=num/tmp 	!旋滚比 滚动体相对于内圈
			else
		    br%balno(idx)%inn%spn2rol=num/(tmp+(omgc-omgi)*cai)	!旋滚比 滚动体相对于内圈
			endif
		endif
!		if(sv) br%balno(:)%inn%spn2rol=(omgb*sin(alpha-alphai)+(omgi-omgc)*sin(alphai))/(omgb*cos(alpha-alphai)-(omgi-omgc)*cos(alphai))	!化简公式	!旋滚比 滚动体相对于内圈
		if(br%hcmthd==1) then
			inn21=0.0_dp
		else
			inn21=omgo*cao/tmp
		endif
	else
		tmp=omgy*cao+omgz*sao
		num	= e*omgo-ro*(tmp-omgo*cao)
		demo=(e*omgo+ro*(tmp+omgo*cao))/2.0_dp		
		xixo=num/demo
		num	=ro*omgx   
		xiyo=num/demo
		num	=(omgz*cao-omgy*sao+omgo*sao)
		psio=num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(idx)%out%spn2rol=num/tmp 	!旋滚比
			else
		    br%balno(idx)%out%spn2rol=num/(tmp-omgo*cao)	!旋滚比
			endif
		endif
		if(br%hcmthd==1) then
			out21=0.0_dp
		else
			out21=omgi*cai/tmp
		endif

		tmp=omgy*cai+omgz*sai
		num = omgi*e-ri*(tmp+omgi*cai)
		demo=(omgi*e+ri*(tmp-omgi*cai))/2.0_dp		
		xixi=num/demo
		num =-ri*omgx  
		xiyi=num/demo
		num	=(omgy*sai-omgz*cai+omgi*sai)
		psii=num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(idx)%inn%spn2rol=num/tmp 	!旋滚比
			else
		    br%balno(idx)%inn%spn2rol=num/(tmp-omgi*cai)	!旋滚比
			endif
		endif
		if(br%hcmthd==1) then
			inn21=0.0_dp
		else
			inn21=omgo*cao/tmp
		endif
	endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	br%balno(idx)%out%xix=xixo
	br%balno(idx)%out%xiy=xiyo
	br%balno(idx)%out%psi=psio
	br%balno(idx)%inn%xix=xixi
	br%balno(idx)%inn%xiy=xiyi
	br%balno(idx)%inn%psi=psii

	call getcij(br%balno(idx)%inn%alf,br%balno(idx)%out%alf)  
	call getcreep_idx(inn21,out21,idx)	!计算fastsim所需的各种参数

	if(br%optin) then
		outer=.true.
		call prefs(idx,outer)
		call fastsim()
		call postfs(idx,outer)
		outer=.false.
		call prefs(idx,outer)  
		call fastsim()
		call postfs(idx,outer)
	else
		outer=.false.
		call prefs(idx,outer)  
		call fastsim()
		call postfs(idx,outer)
		outer=.true.
		call prefs(idx,outer)
		call fastsim()
		call postfs(idx,outer)
	endif

	fc =0.0_dp			!赋初值
	mgx=0.0_dp	!赋初值
	mgz=0.0_dp	!赋初值
	!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
	selectcase(br%inert)
	case(-3)
		mgx=br%ip*omgz*omgc
		mgz=br%ip*omgx*omgc
	case(-2)
		mgx=br%ip*omgz*omgc
	case(-1)
		mgz=br%ip*omgx*omgc
	case(0)
	 !default fc=mgx=mgz=0.0_dp
	case(1)
		fc =br%m1*e*(omgc**2.0_dp)
	case(2)
		fc =br%m1*e*(omgc**2.0_dp)
		mgx=br%ip*omgz*omgc
	case(3)
		fc =br%m1*e*(omgc**2.0_dp)
		mgz=br%ip*omgx*omgc
	case(4)
		fc =br%m1*e*(omgc**2.0_dp)
		mgx=br%ip*omgz*omgc
		mgz=br%ip*omgx*omgc
	endselect

	txo=br%balno(idx)%out%tx
	tyo=br%balno(idx)%out%ty
	mxo=br%balno(idx)%out%mx
	myo=br%balno(idx)%out%my
	mzo=br%balno(idx)%out%mz
	txi=br%balno(idx)%inn%tx
	tyi=br%balno(idx)%inn%ty
	mxi=br%balno(idx)%inn%mx
	myi=br%balno(idx)%inn%my
	mzi=br%balno(idx)%inn%mz

	f(1)=	txo-txi
	f(2)=-p_o*cao+p_i*cai-tyo*sao-tyi*sai+fc
	f(3)=	p_o*sao-p_i*sai-tyo*cao-tyi*cai
!	f(4)=	mxo-mxi-mgx  
	f(4)=	mxo-mxi+mgx		!注意方向
	f(5)=-mzo*cao+mzi*cai-myi*sai-myo*sao+mgz 
	f(6)=	mzo*sao-mzi*sai-myi*cai-myo*cao	

	if(sv) then
		br%balno(idx)%fc=fc
		br%balno(idx)%mgx=mgx
		br%balno(idx)%mgz=mgz
	  br%balno(idx)%e=e					!滚动体节圆半径
		br%balno(idx)%inn%p =p_i 
		br%balno(idx)%out%p =p_o 
		br%balno(idx)%omgc=omgc
		br%balno(idx)%omgx=omgx
		br%balno(idx)%omgy=omgy
		br%balno(idx)%omgz=omgz
		br%balno(idx)%inn%app=dlti
		br%balno(idx)%out%app=dlto
		br%balno(idx)%inn%alf=alfi
		br%balno(idx)%out%alf=alfo
		br%balno(idx)%f(:)=f(:)
	endif
endsubroutine slv6eqs_idx