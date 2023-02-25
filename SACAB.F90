!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!           quasi-static simulation tool								!!!!!
!!!!!         for angular contact ball bearings							!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program sacab
	use fsbr
	use dfport 	
	USE INTFACE,ONLY:RAD2DEG,rad2rev,rev2rad
	use shell32
	implicit none
	type(patchabk)	 ,target	:: abk
	character(len=1)			:: tab=char(9)
	real(4) i,j,k,ta(2)
	character*24 today
	logical	:: t=.false.  !测试或滚动接触计算程序
	real::f_open
	real(dp),allocatable::zeropos(:,:)
	real(dp)::vx,vy,vz,rb,rx,ry
	integer :: num=1
	real(dp)	DNRM2
	real(dp)	alfa,beta,omgx,omgy,omgz,omgb,alfi,alfo,posang,tori,toro,torii,toroi
	logical	:: rep !避免重复计及
	integer	:: typ
	call initparas()  
	call initfs()			
	call initbr()
	call zeros()
	call getstatic()	!计算无摩擦纯轴向载荷静态接触角

	if(t==.false.) then
		open(12, file="final.xls")		!用于输出最终数据,注意物理量的单位
!		call slv01_broydn()					!broydn
!		fs% heathcote		=.true.			!heathcote slip   false
!		call slv02_lbfgsb()					!l-bfgs-b  对于本问题的小维度而言，浪费！
!		call slv03_tensolve()				!tensolve
!		call slv04_sa()							!sa
!		call slv05_gcnewton()				!globally convergent newton
!		call slv06_newton_raphson()	!newton-raphson  !同hybrid,可作为后续精度提高程序
!		call slv08_nl2sol()					!nl2sol  
!		call slv09_smsno()					!子程序值得借鉴
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		call slv10_dneqnf()					!dneqnf    imsl  同hybrid
!		call slv11_dneqnj()					!dneqnj    imsl  同hybrid
!		call slv12_dneqbf()					!broyden's update and finite-difference jacobian   全局收敛性不理想 imsl
!		call slv13_dneqbj()					!broyden's update and analytic jacobian
!	以上二个算法来自imsl，程序提供给用户的控制字段太少；很遗憾看不到源代码。
!		call slv14_alcon2()
!		call slv15_alcon1()					!continuation algorithm by p. deuflhard
!		call slv16_hompack()				!homotopy
!		call slv17_bfgs()						!bfgs method from nr f90
!		call slv18_frprmn()					!from nr 效果很差
!		call slv19_uncmin()			!不理想
!		call slv20_nitsol()		!nitsol
!		call slv23_giant_gbit()  !不能用、不会用!

		do i=1,3
!			call slv21_nleq1()					!
!			call slv22_nleq2()					!
!			call slv07_hybrd()						!minpack
		enddo
		open(13, file="temp.xls")		!用于输出最终数据,注意物理量的单位
		if(br%sym) then !.and. not(br%cmdcreep)) then
			!纯轴向载荷蠕滑计算及数据输出
			do i=1,30
!			do i=1,1
!				br%omgi=50._dp*i
				br%aload=40000._dp
				br%omgi=10000.0_dp
				br%omgi=rev2rad(br%omgi)
				br%mu=0.01_dp*i
!				br%aload=199._dp*500._dp
				call getstatic()	!使每次计算互不影响，计算静平衡参数
				call slv24_neslv()
				br%dlta=(br%ao+br%balno(1)%out%app)*sin(BR%BALNO(1)%out%alf)+(br%ai+br%balno(1)%inn%app)*sin(BR%BALNO(1)%inn%alf)-(0.5_dp*br%d_i+br%ri-0.5_dp*br%d_o+br%ro)*tan(br%alf0)
				write(13,"(57(e24.16,a1))")			&
					rad2rev(br%omgi),tab,	&		!内圈相对转速
					br%balno(1)%est,tab,&			!静平衡节圆半径
					br%dltast,tab,&						!静平衡内圈轴向位移
					rad2rev(abs(br%omgi)+abs(BR%omgo_ictl)),tab,	&		!内控内圈绝对运动
					br%balno(1)%e0_ictl,tab,&			!内控节圆半径
					br%dlta0_ictl,tab,&						!内控轴向位移
					rad2rev(BR%omgo_ictl),tab,&		!内控外圈转速
					rad2rev(BR%omgb_ictl),tab,&		!内控滚动体转速
					RAD2DEG(BR%iconang_ictl),tab,&	!内控内接触角
					RAD2DEG(BR%oconang_ictl),tab,&	!内控外接触角
					RAD2DEG(br%pchang_ictl),tab,&	!内控俯仰角
					br%s2rout,tab,&	!内滚道控制，外接触旋滚比
					rad2rev(abs(br%omgi)+abs(BR%omgo_octl)),tab,	&		!内控外圈绝对运动
					br%balno(1)%e0_octl,tab,&			!外控节圆半径
					br%dlta0_octl,tab,&						!外控轴向位移
					rad2rev(BR%omgo_octl),tab,&		!外控外圈转速
					rad2rev(BR%omgb_octl),tab,&		!外控滚动体转速
					RAD2DEG(BR%iconang_octl),tab,&	!外控内接触角
					RAD2DEG(BR%oconang_octl),tab,&	!外控外接触角
					RAD2DEG(br%pchang_octl),tab,&	!外控俯仰角
					br%s2rinn,tab,&	!外滚道控制，内接触旋滚比
					0.0_dp,tab,&							!不考虑摩擦偏航角
					rad2rev(abs(br%omgi)+abs(BR%x(1))),tab,	&		!蠕滑内圈绝对运动
					br%balno(1)%e,tab, &		!滚动体节圆半径
					br%dlta,tab,&						!动平衡内圈轴向位移
					rad2rev(BR%x(1)),tab,&  !相对运动，外圈的相对运动
					RAD2DEG(BR%x(5)),tab,&  !蠕滑内接触角
					RAD2DEG(BR%x(6)),tab,&  !蠕滑外接触角
!					RAD2DEG(atan(abs(BR%x(4))/SQRT(BR%x(2)**2+BR%x(3)**2))),tab,&	!俯仰角 绝对值
!					RAD2DEG(atan(BR%x(2)/abs(BR%x(3)))),tab,&		!偏航角
					RAD2DEG(atan(BR%x(4)/SQRT(BR%x(2)**2+BR%x(3)**2))),tab,&	!俯仰角 绝对值
					RAD2DEG(atan(BR%x(2)/BR%x(3))),tab,&		!偏航角
					RAD2DEG(BR%x(5)-atan(sin(br%x(5))/(cos(br%x(5))-br%balno(1)%e/(br%r-.5_dp*br%balno(1)%inn%app)))),tab,&	  !理论内控俯仰角
					RAD2DEG(BR%x(6)-atan(sin(br%x(6))/(cos(br%x(6))+br%balno(1)%e/(br%r-.5_dp*br%balno(1)%out%app)))),tab,&	  !理论外控俯仰角
					rad2rev(sqrt(BR%x(2)**2+BR%x(3)**2+BR%x(4)**2)),tab, &
					rad2rev(BR%x(2)),tab, &
					rad2rev(BR%x(3)),tab, &
					rad2rev(BR%x(4)),tab, &
					BR%BALNO(1)%inn%spn2rol,tab,&
					BR%BALNO(1)%out%spn2rol,tab,&
					br%z*(br%balno(1)%inn%tx*(br%balno(1)%e-(br%r-br%balno(1)%inn%app/2.0_dp)*cos(br%balno(1)%inn%alf))+br%balno(1)%inn%mz*sin(br%balno(1)%inn%alf)),tab,&
					br%z*(br%balno(1)%out%tx*(br%balno(1)%e+(br%r-br%balno(1)%out%app/2.0_dp)*cos(br%balno(1)%out%alf))+br%balno(1)%out%mz*sin(br%balno(1)%out%alf)),tab,&
					BR%BALNO(1)%OUT%TX,	tab,&
					BR%BALNO(1)%OUT%TY,	tab,&
					BR%BALNO(1)%OUT%P,	tab,&
					BR%BALNO(1)%OUT%MX,	tab,&
					BR%BALNO(1)%OUT%MY,	tab,&
					BR%BALNO(1)%OUT%MZ,	tab,&
					SQRT(BR%BALNO(1)%OUT%TX**2+BR%BALNO(1)%OUT%TY**2),tab,&
					BR%BALNO(1)%INN%TX,	tab,&
					BR%BALNO(1)%INN%TY,	tab,&
					BR%BALNO(1)%INN%P,	tab,&
					BR%BALNO(1)%INN%MX,	tab,&
					BR%BALNO(1)%INN%MY,	tab,&
					BR%BALNO(1)%INN%MZ,	tab,&
					sqrt(BR%BALNO(1)%INN%TX**2+BR%BALNO(1)%INN%TY**2),tab,&
					br%balno(1)%mgx,tab,&
					br%balno(1)%mgz,tab,&
					DNRM2(size(br%balno(1)%f),br%balno(1)%f,1),tab
			enddo
		else
			!联合载荷蠕滑计算及数据输出
			typ=1
			do k=1,1
!				br%omgi=1000._dp+000._dp*k
				if(typ==1) then
					br%omgi		=10000._dp
					br%aload	=40000._dp
					br%rload	=4000.0_dp
					br%mload	=300.0_dp
					br%mu			=0.025_dp*k
				elseif(typ==2) then
					br%omgi		=12000._dp
					br%aload	=40000._dp
					br%rload	=1000.0_dp*k
					br%mload	=300.0_dp
					br%mu			=0.1_dp
				elseif(typ==3) then
					br%omgi		=1000._dp*k
					br%aload	=40000._dp
					br%rload	=5000.0_dp
					br%mload	=300.0_dp
					br%mu			=0.1_dp
				elseif(typ==4) then
					br%omgi		=12000._dp
					br%aload	=40000._dp
					br%rload	=1000.0_dp
					br%mload	=920.0_dp*k
					br%mu			=0.1_dp
				endif
				write(13,"(a5,a1,e23.15)") 'aload',tab,br%aload
				write(13,"(a5,a1,e23.15)") 'rload',tab,br%rload
				write(13,"(a5,a1,e23.15)") 'mload',tab,br%mload
				write(13,"(a5,a1,e23.15)") 'omgi'	,tab,br%omgi
				write(13,"(a5,a1,e23.15)") 'mu'		,tab,br%mu
				write(13,"(a2,a1,i5,a1,a1,i5)") 'nx',tab,fs%nx,'ny',tab,fs%ny
				br%omgi=rev2rad(br%omgi)
				call getstatic()	!使每次计算互不影响，计算静平衡参数
				call slvcmbd()
				write(*,*) '.......               FINISHED!              ..........'
				write(*,*) 'x'
				write(*,*) br%dlta,br%dltr,br%theta
				write(*,*) 'f'
				write(*,*) br%f3
				write(*,*) dnrm2(size(br%f3),br%f3,1)
				rep=.false.
				tori=0.0_dp
				toro=0.0_dp
				do j=1,br%z+1
					i=j
					if(i==br%z+1) then
						i=1
						rep=.true.
					endif
					omgx=rad2rev(br%balno(i)%omgx)
					omgy=rad2rev(br%balno(i)%omgy)
					omgz=rad2rev(br%balno(i)%omgz)
					omgb=sqrt(omgx**2+omgy**2+omgz**2)
					alfa=rad2deg(asin(omgz/omgb))
					beta=rad2deg(atan(omgx/omgy))
					alfi=br%balno(i)%inn%alf
					alfo=br%balno(i)%out%alf
					if(br%balno(i)%posang<0.0_dp .or. int(j)==br%z+1) then
						posang=br%balno(i)%posang+TWOPI_D
					else
						posang=br%balno(i)%posang
					endif
					torii=br%balno(i)%inn%tx*(br%balno(i)%e-(br%r-br%balno(i)%inn%app/2.0_dp)*cos(br%balno(i)%inn%alf))+br%balno(i)%inn%mz*sin(br%balno(i)%inn%alf)
					toroi=br%balno(i)%out%tx*(br%balno(i)%e+(br%r-br%balno(i)%out%app/2.0_dp)*cos(br%balno(i)%out%alf))+br%balno(i)%out%mz*sin(br%balno(i)%out%alf)
					if(rep==.false.) then
						tori=tori+torii
						toro=toro+toroi
					endif
					write(13,"(i2,a1,\)") int(j),tab
					write(13,"(39(e24.16,a1))")	&
													posang,tab,						&
													rad2deg(posang),tab,						&
													rad2rev(br%omgi),tab,&
													rad2rev(br%balno(i)%omgc),tab,						&
													omgb,tab,						&
													omgx,tab,						&
													omgy,tab,						&
													omgz,tab,						&
													alfa,tab,						&
													beta,tab,						&
													rad2deg(alfi),tab,		&
													rad2deg(alfo),tab,		&
													rad2deg(alfi-atan(sin(alfi)/(cos(alfi)-br%balno(i)%e/(br%r-.5_dp*br%balno(i)%inn%app)))),tab,&	  !理论内控俯仰角
													rad2deg(alfo-atan(sin(alfo)/(cos(alfo)+br%balno(i)%e/(br%r-.5_dp*br%balno(i)%out%app)))),tab,&	  !理论外控俯仰角
													br%balno(i)%e,tab,									&
													br%balno(i)%inn%app,tab,						&
													br%balno(i)%out%app,tab,						&							
													br%balno(i)%out%tx,	tab,&
													br%balno(i)%out%ty,	tab,&
													br%balno(i)%out%p,	tab,&
													br%balno(i)%out%mx,	tab,&
													br%balno(i)%out%my,	tab,&
													br%balno(i)%out%mz,	tab,&
													sqrt(br%balno(i)%out%tx**2+br%balno(i)%out%ty**2),tab,&
													br%balno(i)%inn%tx,	tab,&
													br%balno(i)%inn%ty,	tab,&
													br%balno(i)%inn%p,	tab,&
													br%balno(i)%inn%mx,	tab,&
													br%balno(i)%inn%my,	tab,&
													br%balno(i)%inn%mz,	tab,&
													sqrt(br%balno(i)%inn%tx**2+br%balno(i)%inn%ty**2),tab,&
													br%balno(i)%mgx,tab,								&
													br%balno(i)%mgz,tab,								&
													br%balno(i)%fc,tab,									&
													abs(br%balno(i)%inn%spn2rol),tab,				&
													abs(br%balno(i)%out%spn2rol),tab,				&
													torii,tab,&
													toroi,tab,&
													dnrm2(size(br%balno(i)%f),br%balno(i)%f,1)
				enddo
				write(13,"(6(a14,a1))") "轴向位移",tab,"径向位移",tab,"翻转角度(deg)",tab,"内圈阻力矩",tab,"外圈阻力矩",tab,"整体方程组余量",tab
				write(13,"(6(e23.15,a1))") br%dlta,tab,br%dltr,tab,rad2deg(br%theta),tab,tori,tab,toro,tab,dnrm2(size(br%f3),br%f3,1),tab
				write(13,*)
			enddo
		endif
		close(13)
		f_open=ShellExecute(0, "open", "temp.xls", NULL, NULL, 1)

!		call slv25_nleqslv()		!作为后续处理，可以达到极高的精度
		if(br%cmdcreep) then
			call wrttm(1)  !输出拖动力、拖动力矩
			call optslu()	!输出求解数据，用于论文图表
			call optlog()  !放置在程序的最后
!			call print_sfield()
!			call print_tfield()
!			call vector_combining()  !矢量合成x、y，并将结果存入前一变量
!			call print_field()
		endif
	!!!!!!!!!!!  test  !!!!!!!!!!!!
	!以下均为测试用数据
	elseif(t==.true.) then
!		open(6, file="test.dat")		!创建数据文件 excel文件
		open(6, file="test.xls")		!创建数据文件 excel文件
!!		call tt()	!方程求解算法测试程序
		write(*,*) 
		write(*,*) "!!!!!!!!!!!  test  !!!!!!!!!!!!"
		num=3
		if(num==2) then
		!		call ftest(1,0.0_dp,0.004_dp,0.0_dp) !测试程序heathcote=false hertz=true
			!将数据输出到TECPLOT格式的数据文件
		!	WRITE(6,*) "TITLE = ""DISTRIBUTION"""
		!	WRITE(6,*) "VARIABLES = ""X"",""Y"",""DISTRIBUTION"""
		!	WRITE(6,*) "ZONE T=""DISTRIBUTION"""
		!	WRITE(6,"(' I=',I5,', J=',I5,', K=1',', ZONETYPE=Ordered')") 41,41
		!	WRITE(6,*) "DATAPACKING=POINT"
		!	WRITE(6,*) "DT=(DOUBLE DOUBLE DOUBLE )"
		!		do i=0,40
		!		 do j=0,40
		!		 write(*,*) i,j
					vx=0.00000010_dp
					vy=0.00000000_dp
					vz=0.0000001_dp
		!			call ftest2(1._dp,1._dp,.00_dp,0.00000001_dp*i,0.00000_dp) !测试程序heathcote=false hertz=true
		!			call ftest2(1._dp,1._dp,0.000000005_dp*i,0.00000001_dp*j,0.0000004_dp) !测试程序heathcote=false hertz=true
					call ftest2(1._dp,1._dp,vx,vy,vz) !测n                 y7777i8&ZA/zz/.ZZz试程序heathcote=false hertz=true
		!			write(6,"(5(e23.15,a1))") br%tmp1,tab,br%tmp2,tab,abs(br%ai),tab,abs(br%ao),tab,br%aj,tab
					write(6,"(6(e23.15,a1))") vx,tab,br%tmp1,tab,vy,tab,br%tmp2,tab,vz,tab,br%aj,tab
		!			enddo
		!		enddo
		!		write(*,*) fs%tx,fs%ty
		!		call adh3(fs,br) !测试程序heathcote=false hertz=true
		!		call spiny(fs,br) !测试程序heathcote=false hertz=true
		!		call ellip(br)
				allocate(zeropos(-fs%ny-1:fs%ny+1,4))
				do i=fs%ny,-fs%ny,-1
					if(i==0) cycle
					do j=fs%nx,-fs%nx,-1
						if(j==0) cycle
						if(fs%sx(i,j)/=0.0_dp .or. fs%sy(i,j)/=0.0_dp) then
							zeropos(i,1)=fs%xc(i,j)
							zeropos(i,2)=fs%yc(i,j)
		!					write(6,"(2(e23.15,a1))") zeropos(i,1),tab,zeropos(i,2),tab
							exit
						endif
					enddo
				enddo
					call print_sfield()
					call print_tfield()
					call vector_combining()  !矢量合成x、y，并将结果存入前一变量
					call print_field()
				close(6)
		!		f_open=ShellExecute(0, "open", "test.xls", NULL, NULL, 1)
		elseif(num==3) then
			FS% HEATHCOTE	=.true.		!HEATHCOTE SLIP   FALSE
			br%hcmthd			=1				!1-johnson最简单;2-johnson+mzk;3-mzk最复杂
			FS% HERTZ		=.TRUE.			!Hertz Contact Model
			BR% HERTZ		=FS% HERTZ		!Hertz Contact Model
			BR% SIMPMOMENT	=.TRUE.			!采用简化方法计算作用于球心的力矩

			br%g1				=8.4e10_dp
			br%nu1			=0.25_dp
			br%e1				=br%g1*2.0_dp*(1.0_dp+br%nu1)
			br%g2				=br%g1
			br%nu2			=br%nu1
			br%e2				=br%e1
			BR%BALNO(:)%OUT%P=1.0E5_DP
			BR%BALNO(:)%INN%P=BR%BALNO(1)%OUT%P
			BR%MU=0.1_DP

			rb=1.0e-2_dp
			rx=huge(0.0_dp)
			ry=-1.05e-2_dp
			br%r=rb
			br%ocom%rx2	=rx
			br%ocom%ry2	=ry
			abk%rx1		=rb
			abk%ry1		=rb
			abk%rx2		=rx
			abk%ry2		=ry
			abk%kp		=br%kp
			abk%e1		=br%e1	!quasiidentical
			abk%e2		=br%e2	!quasiidentical
			abk%nu1		=br%nu1
			abk%nu2		=br%nu1
			
			call getabk(abk)
			br%ocom%ke(:)			=abk%ke
			br%ocom%and(:)		=abk%and
			br%ocom%bnd(:)		=abk%bnd
			br%ocom%xlarger		=abk%xlarger
			call getpchaxis_idx(1,0.0_dp,0.0_dp)
		!将数据输出到TECPLOT格式的数据文件
		!	WRITE(6,*) "TITLE = ""DISTRIBUTION"""
		!	WRITE(6,*) "VARIABLES = ""X"",""Y"",""DISTRIBUTION"""
		!	WRITE(6,*) "ZONE T=""DISTRIBUTION"""
		!	WRITE(6,"(' I=',I5,', J=',I5,', K=1',', ZONETYPE=Ordered')") 41,41
		!	WRITE(6,*) "DATAPACKING=POINT"
		!	WRITE(6,*) "DT=(DOUBLE DOUBLE DOUBLE )"
!				do i=0,10000,1
!				 do j=0,40
					vx=-0.03_dp
					vy=0.01_dp
					vz=20.0_dp
					call ftest3(vx,vy,vz) 
					write(6,"(8(e23.15,a1))") vx,tab,br%tmp1,tab,vy,tab,br%tmp2,tab,vz,tab,br%aj,tab,br%ai,tab,br%ao,tab
					write(6,"(1(e23.15,a1))") abk%ke,tab
		!			enddo
!				enddo
				allocate(zeropos(-fs%ny-1:fs%ny+1,4))
!				do i=fs%ny,-fs%ny,-1
!					if(i==0) cycle
!					do j=fs%nx,-fs%nx,-1
!						if(j==0) cycle
!						if(fs%sx(i,j)/=0.0_dp .or. fs%sy(i,j)/=0.0_dp) then
!							zeropos(i,1)=fs%xc(i,j)
!							zeropos(i,2)=fs%yc(i,j)
		!					write(6,"(2(e23.15,a1))") zeropos(i,1),tab,zeropos(i,2),tab
!							exit
!						endif
!					enddo
!				enddo
					call print_sfield()
					call print_tfield()
					call vector_combining()  !矢量合成x、y，并将结果存入前一变量
					call print_field()
				close(6)
				f_open=ShellExecute(0, "open", "test.xls", NULL, NULL, 1)
		endif
	endif
!!!!!!!!!!!  test  !!!!!!!!!!!!

	call freemem()	!释放内存
	i=etime(ta)
	call fdate(today)
	write(*,*) 
	write(*,*) 'program has used', i, 'seconds of cpu time.'
	write(*,*) 'today is ', today
end program sacab
