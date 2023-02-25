subroutine slv6eqs(x,f,sv)
	use fsbr
	use nrtype
	implicit none
	real(dp), dimension(:) :: x
	real(dp), dimension(:) :: f
	logical				::sv
	character(len=1)			:: tab=char(9)
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
	real(dp)	::temp
	real(dp)	::inn21,out21

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
	alfi=x(5)
	alfo=x(6)
	cai=cos(alfi)
	sai=sin(alfi)
	cao=cos(alfo)
	sao=sin(alfo)
	tyo=br%balno(1)%out%ty
	tyi=br%balno(1)%inn%ty

	p_o=(br%aload/br%z+tyo*cao)/sao
	p_i=(br%aload/br%z-tyi*cai)/sai
	br%balno(:)%out%p =p_o 
	br%balno(:)%inn%p =p_i 
	call getpchaxis(alfi,alfo)

	dlto=br%balno(1)%out%app
	dlti=br%balno(1)%inn%app
	ro=br%r-dlto/2.0_dp !接触中心变形后与滚动体中心的距离
	ri=br%r-dlti/2.0_dp !同上

	e=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+br%balno(1)%inn%app)*cai
	if(sv) br%balno(:)%e=e	!滚动体节圆半径
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!creepage and spin
	if(br%absolute) then  !绝对速度
		temp=omgy*cao+omgz*sao
		num	= e*(omgc-omgo)+ro*(temp-omgo*cao+omgc*cao)
		demo=(e*(omgo-omgc)+ro*(temp+omgo*cao-omgc*cao))/2.0_dp		
		xixo=-num/demo
		num	=ro*omgx   
		xiyo=-num/demo
		num	=(omgy*sao-omgz*cao-omgo*sao +omgc*sao)
		psio=-num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(:)%out%spn2rol=num/temp	!详见学位论文
			else
		    br%balno(:)%out%spn2rol=num/(temp+(omgc-omgo)*cao)	!详见学位论文
			endif
		endif
!		if(sv) br%balno(:)%out%spn2rol=(omgb*sin(alphao-alpha)+(omgc-omgo)*sin(alphao))/(omgb*cos(alphao-alpha)+(omgc-omgo)*cos(alphao))	!化简公式
		if(br%hcmthd==1) then
			out21=0.0_dp
		else
			out21=omgi*cai/temp
		endif

		temp=omgy*cai+omgz*sai
		num = (omgi-omgc)*e+ri*(temp-omgi*cai+omgc*cai)
		demo=((omgi-omgc)*e-ri*(temp+omgi*cai-omgc*cai))/2.0_dp		
		xixi=num/demo
		num =ri*omgx  
		xiyi=-num/demo
		num	=(omgz*cai-omgy*sai+omgi*sai -omgc*sai)
		psii=num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(:)%inn%spn2rol=num/temp 	!旋滚比 滚动体相对于内圈
			else
		    br%balno(:)%inn%spn2rol=num/(temp+(omgc-omgi)*cai)	!旋滚比 滚动体相对于内圈
			endif
		endif
!		if(sv) br%balno(:)%inn%spn2rol=(omgb*sin(alpha-alphai)+(omgi-omgc)*sin(alphai))/(omgb*cos(alpha-alphai)-(omgi-omgc)*cos(alphai))	!化简公式	!旋滚比 滚动体相对于内圈
		if(br%hcmthd==1) then
			inn21=0.0_dp
		else
			inn21=omgo*cao/temp
		endif
	else
		temp=omgy*cao+omgz*sao
		num	= e*omgo-ro*(temp-omgo*cao)
		demo=(e*omgo+ro*(temp+omgo*cao))/2.0_dp		
		xixo=num/demo
		num	=ro*omgx   
		xiyo=num/demo
		num	=(omgz*cao-omgy*sao+omgo*sao)
		psio=num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(:)%out%spn2rol=num/temp 	!旋滚比
			else
		    br%balno(:)%out%spn2rol=num/(temp-omgo*cao)	!旋滚比
			endif
		endif
		if(br%hcmthd==1) then
			out21=0.0_dp
		else
			out21=omgi*cai/temp
		endif

		temp=omgy*cai+omgz*sai
		num = omgi*e-ri*(temp+omgi*cai)
		demo=(omgi*e+ri*(temp-omgi*cai))/2.0_dp		
		xixi=num/demo
		num =-ri*omgx  
		xiyi=num/demo
		num	=(omgy*sai-omgz*cai+omgi*sai)
		psii=num/demo
		if(sv) then
			if(br%s2rmthd==1) then  !no consideration omgi
		    br%balno(:)%inn%spn2rol=num/temp	!旋滚比
			else
		    br%balno(:)%inn%spn2rol=num/(temp-omgi*cai)	!旋滚比
			endif
		endif
		if(br%hcmthd==1) then
			inn21=0.0_dp
		else
			inn21=omgo*cao/temp
		endif
	endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	br%balno(:)%out%xix=xixo
	br%balno(:)%out%xiy=xiyo
	br%balno(:)%out%psi=psio
	br%balno(:)%inn%xix=xixi
	br%balno(:)%inn%xiy=xiyi
	br%balno(:)%inn%psi=psii

	call getcij(alfi,alfo)  
	call getcreep3(inn21,out21)	!计算fastsim所需的各种参数

	!循环计算接触拖动力
	do i=1,br%cyc
		if(br%optin) then
			outer=.true.
			call prefs(i,outer)
			call fastsim()
			call postfs(i,outer)
			outer=.false.
			call prefs(i,outer)  
			call fastsim()
			call postfs(i,outer)
		else
			outer=.false.
			call prefs(i,outer)  
			call fastsim()
			call postfs(i,outer)
			outer=.true.
			call prefs(i,outer)
			call fastsim()
			call postfs(i,outer)
		endif
	enddo 

		!拟静力学因素
		!离心力
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
	if(sv) then
		br%balno(:)%fc=fc
		br%balno(:)%mgx=mgx
		br%balno(:)%mgz=mgz
		br%balno(:)%inn%alf=alfi
		br%balno(:)%out%alf=alfo
	endif

	txo=br%balno(1)%out%tx
	tyo=br%balno(1)%out%ty
	mxo=br%balno(1)%out%mx
	myo=br%balno(1)%out%my
	mzo=br%balno(1)%out%mz
	txi=br%balno(1)%inn%tx
	tyi=br%balno(1)%inn%ty
	mxi=br%balno(1)%inn%mx
	myi=br%balno(1)%inn%my
	mzi=br%balno(1)%inn%mz


!	f(1)=(mxo-mxi-mgx)
	f(1)=(mxo-mxi+mgx)		!注意方向
	f(2)=(-mzo*cao+mzi*cai-myi*sai-myo*sao+mgz +mzo*sao-mzi*sai-myi*cai-myo*cao)
	f(3)=(txo-txi)
	f(4)=(-p_o*cao+p_i*cai-tyo*sao-tyi*sai+fc)
	f(5)=(txi*(e-ri*cai)+mzi*sai-txo*(e+ro*cao)-mzo*sao) 
	f(6)=((br%ao+dlto)*cao+(br%ai+dlti)*cai)/(br%a*cos(br%alf0))-1.0_dp
	!	funcv(7)=br%constt*(br%aload/br%z-p_o*sao+tyo*cao)
	!	funcv(8)=br%constt*(br%aload/br%z-p_i*sai-tyi*cai)
!	write(6,"(e23.15,a1,e23.15,a1,e23.15,a1,e23.15,a1,e23.15,a1,e23.15,a1,e23.15)")  &
	!	omgo,tab,omgb,tab,alpha,tab,beta,tab,alphai,tab,alphao,tab,p_i,tab,p_o
!		x(1),tab,x(2),tab,x(3),tab,x(4),tab,x(5),tab,x(6),tab,0.5_dp*(x(5)+x(6))
	!	txi,tab,tyi,tab,xiyo,tab,xiyi,tab,psio,tab,psii,tab,funcv(6),tab,funcv(8)
endsubroutine slv6eqs