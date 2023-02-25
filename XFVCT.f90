subroutine x0vct(x)
	use fsbr
	use nr
	use intface,only:pow,rad2deg,deg2rad
	implicit none
	real(dp), dimension(:) :: x
	real(dp) ::est,uist,uost
	real(dp)::pchi,pcho,omgb,omgo,omgc
	real(dp)::t1,t2
	real(dp)::pchang
	real(dp)::e,rdi,rdo
	real(dp)::alfst,appndi,appndo
	real(dp)::dy
	integer	 ::n,i

	n=size(x)

	!将br某些在计算过程中改变且会影响到下一轮计算的变量置零；
	!此点对于采用多重算法相互递进的调用非常重要；
!	call zerobr()  !该代码的取舍非常关键。
	alfst=br%balno(1)%inn%alfst
	if(all(x==0.0_dp)) then
		!mzk 近似计算
		if(br%absolute) then	  	
			!一般工况
			!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
			if(br%inert>=1) then	
				!一般工况、考虑惯性力
				x(5)=br%balno(1)%inn%alfst
				x(6)=br%balno(1)%out%alfst
				call conangs_sym(x(5),x(6))
				omgc=br%balno(1)%omgc
				omgb=br%balno(1)%omgb
				pchang=br%pchang
				x(1)=omgc
				x(2)=-0.001_dp
				x(3)=omgb*cos(pchang)
				x(4)=omgb*sin(pchang)  
			else	
				!一般工况、不考虑惯性力
				!1:mzk 近似方法，最不好；2：内滚道控制理论,内圈主动：3：外滚道控制理论，内圈主动
				call ratint(br%alf,br%icom%appnd,alfst,appndi,dy)
				call ratint(br%alf,br%ocom%appnd,alfst,appndo,dy)
				uist=appndi*pow(br%balno(1)%inn%p,2.0_dp/3.0_dp)/2.0_dp
				uost=appndo*pow(br%balno(1)%out%p,2.0_dp/3.0_dp)/2.0_dp
				est=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+uist*2.0_dp)*cos(alfst)
				if(br%rcm==1) then
					omgc=((est-(br%r-uist)*cos(alfst))*br%omgi+(est+(br%r-uost)*cos(alfst))*br%omgo)/2.0_dp/est
					omgb=(omgc-br%omgi)*(est/(br%r-uist)-cos(alfst))
					pchang=alfst
				elseif(br%rcm==2) then
					!内滚道控制!pchi 内滚道控制时滚动体的姿态角
					pchi=atan(est*sin(alfst)/(est*cos(alfst)-(br%r-uist)))
					omgc=(est*(br%omgo+br%omgi)+(br%r-uist)*(br%omgo-br%omgi)*cos(alfst))/2.0_dp/est
					omgb=(omgc-br%omgi)*sin(alfst)/(sin(pchi)*cos(alfst)-cos(pchi)*sin(alfst))
					pchang=pchi
				elseif(br%rcm==3) then
					!外滚道控制
					!pcho 外滚道控制时滚动体的姿态角
					pcho=atan(est*sin(alfst)/(est*cos(alfst)+(br%r-uost)))
					omgc=(est*(br%omgo+br%omgi)+(br%r-uist)*(br%omgo-br%omgi)*cos(alfst))/2.0_dp/est
					omgb=(omgc-br%omgo)*sin(alfst)/(sin(pcho)*cos(alfst)-cos(pcho)*sin(alfst))
					pchang=pcho
				endif
				x(1)=omgc
				x(2)=0.001_dp
				x(3)=omgb*cos(pchang)	
				x(4)=omgb*sin(pchang)
				x(5)=alfst
				x(6)=alfst
			endif
		else  
			!特殊工况
			!-3:mgx+mgz;-2:mgx;-1:mgz;0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
			if(br%inert>=1) then	
				!以下代码用于计算旋滚比
				i=br%rcm  !保存
				!外滚道控制，内接触旋滚比
				br%rcm=3  !设置为外控
				x(5)=br%balno(1)%inn%alfst
				x(6)=br%balno(1)%out%alfst
				call conangs_sym(x(5),x(6))
				pcho=br%pchang	!与上式非常接近
				omgb=br%balno(1)%omgb
!				t1=omgb*(cos(pcho)*sin(x(5))-sin(pcho)*cos(x(5)))+br%omgi*sin(x(5))
!				t2=omgb*(cos(pcho)*cos(x(5))+sin(pcho)*sin(x(5)))
!				br%s2rinn=t1/t2
				br%s2rinn=-tan(pcho-x(5))+sin(x(5))/(br%balno(1)%e/(br%r-br%balno(1)%inn%app*.5_dp)-cos(x(5)))
				br%balno(:)%e0_octl=br%balno(1)%e
				br%omgo_octl=br%omgo
				br%omgb_octl=omgb
				br%pchang_octl=pcho
				br%iconang_octl=x(5)
				br%oconang_octl=x(6)
				br%dlta0_octl=br%dlta0

				!内滚道控制，外接触旋滚比
				br%rcm=2  !设置为内控
				call conangs_sym(x(5),x(6))
				PCHI=br%pchang
				omgb=br%balno(1)%omgb
				omgo=br%omgo
!				t1=omgb*(sin(pchi)*cos(x(6))-cos(pchi)*sin(x(6)))+omgo*sin(x(6))
!				t2=omgb*(cos(pchi)*cos(x(6))+sin(pchi)*sin(x(6)))
!				br%s2rout=t1/t2
				br%s2rout=tan(pchi-x(6))+sin(x(6))/(br%balno(1)%e/(br%r-br%balno(1)%out%app*.5_dp)+cos(x(6)))
				br%balno(:)%e0_ictl=br%balno(1)%e
				br%omgo_ictl=br%omgo
				br%omgb_ictl=omgb
				br%pchang_ictl=pchi
				br%iconang_ictl=x(5)
				br%oconang_ictl=x(6)
				br%dlta0_ictl=br%dlta0

				!回复初始状态，继续正常计算
				br%rcm=i
				!特殊工况、考虑惯性力
				x(5)=br%balno(1)%inn%alfst
				x(6)=br%balno(1)%out%alfst
				call conangs_sym(x(5),x(6))
				omgo=br%omgo
				omgb=br%balno(1)%omgb
				pchang=br%pchang
				x(1)=omgo
				x(2)=-0.001_dp
				x(3)=omgb*cos(pchang)
				x(4)=omgb*sin(pchang)  
				if(n>6) then
					x(7)=br%balno(1)%inn%p
					x(8)=br%balno(1)%out%p
				endif
			else  	
				!特殊工况、不考虑惯性力
				!1:mzk 近似方法，最不好；2：内滚道控制理论,内圈主动：3：外滚道控制理论，内圈主动
				if(br%rcm==1) then
					uist=appndi*pow(br%balno(1)%inn%p,2.0_dp/3.0_dp)/2.0_dp
					uost=appndo*pow(br%balno(1)%out%p,2.0_dp/3.0_dp)/2.0_dp
					est=abs(br%icom%rx2)+abs(br%icom%ry2)-(br%ai+uist*2.0_dp)*cos(alfst)
					omgb=br%omgi*(est-(br%r-uist)*cos(alfst))/(br%r-uist)
					omgo=omgb*(br%r-uost)/(est+(br%r-uost)*cos(alfst))
					pchang=alfst
				elseif(br%rcm==2) then
					!内滚道控制
					!pchi 内滚道控制时滚动体的姿态角
					!计算公式参考harris著作，注意计算程序中的物理量的彼此代换
					pchi=atan(br%dm/2.0*sin(br%balno(1)%inn%alfst)/(br%dm/2.0*cos(br%balno(1)%inn%alfst)-br%r))
					omgb=br%omgi*(-br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(pchi)*cos(br%balno(1)%inn%alfst)+sin(pchi)*sin(br%balno(1)%inn%alfst))
					omgo=omgb/((br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(br%balno(1)%inn%alfst)*cos(pchi)+sin(pchi)*sin(br%balno(1)%inn%alfst)))
					!内滚道控制外圈角速度
					omgo=abs(omgo)
					omgb=abs(omgb)
					pchang=pchi
				elseif(br%rcm==3) then
					!外滚道控制
					!pcho 外滚道控制时滚动体的姿态角
					!计算公式参考harris著作，注意计算程序中的物理量的彼此代换
					pcho=atan(br%dm/2.0*sin(br%balno(1)%inn%alfst)/(br%dm/2.0*cos(br%balno(1)%inn%alfst)+br%r))
					omgb=br%omgi*(-br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(pcho)*cos(br%balno(1)%inn%alfst)+sin(pcho)*sin(br%balno(1)%inn%alfst))
					omgo=omgb/((br%dm/2.0+br%r*cos(br%balno(1)%inn%alfst))/	&
											br%r/(cos(br%balno(1)%inn%alfst)*cos(pcho)+sin(pcho)*sin(br%balno(1)%inn%alfst)))
					!外滚道控制外圈角速度'
					omgo=abs(omgo)
					omgb=abs(omgb)
					pchang=pcho
				endif
				x(1)=omgo
				x(2)=0.001_dp
				x(3)=omgb*cos(pchang)
				x(4)=omgb*sin(pchang)
				x(5)=alfst
				x(6)=alfst
				if(n>6) then
					x(7)=br%balno(1)%inn%p
					x(8)=br%balno(1)%out%p
				endif
			endif
		endif

		!保存理想拟静力学分析参数，以备后用．
		br%x0(:)=x(:)
		br%x0(2)=0.0_dp	!经典理论方法，偏航角不存在．
	else
		x(:)=br%x(:)
	endif
end subroutine x0vct

subroutine printx0(x)
	use fsbr
	use intface,only:rad2deg
	implicit none
	integer :: i
	real(dp), dimension(br%n) :: x

	if(pntxo) then
		write(*,*) '    *************************初     值***********************'
		write(*,*) '        index','          x'

		do i=1,size(x)
			if(i==5 .or. i==6) then
				write(*,*) i,x(i),'rad',rad2deg(x(i)),'deg'
			else
				write(*,*) i,x(i)
			endif
		enddo

		write(*,*) '    *********************************************************'
		pntxo=.false.  !只输出第一个求解算法的初值；
	endif
endsubroutine printx0

subroutine printxf(x,f)
	use fsbr
	implicit none
	integer :: i
	real(dp), dimension(br%n) :: x
	real(dp), dimension(br%m) :: f
	real(dp) DNRM2

	write(*,*) '        index','          x','                     f'
	do i=1,size(x)
		write(*,*) i,x(i),f(i)
	enddo
	write(*,*) 'norm of F(x):',DNRM2(size(x),f,1)
	write(*,*) 'omgb',sqrt(x(3)**2+x(4)**2)
endsubroutine printxf

subroutine savexf(x,f)
	use fsbr
	implicit none
	integer :: i,j
	real(dp), dimension(br%n) :: x,f

	do i=1,br%n
		br%x(i)		=x(i)
		do j=1,br%z
			br%balno(j)%f(i)		=f(i)
		enddo
	enddo
endsubroutine savexf
