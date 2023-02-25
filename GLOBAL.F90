!仍然采用jj kalker的网格划分方法，对于计算而言，此种方法对计算要求量最低
!均匀网格的弊端在于边界的处理较为繁琐，需要进行更多的计算。
!边界点的处理比较特殊，不许进行单独处理
module fsbr
	use nrtype
	implicit none
	logical			:: pntxo
	logical			:: maxit
	logical			:: embed
	real(dp)		:: mprec,sqmprec
	type:: fastdata    !fastsim算法数据结构
 		integer		:: nx,ny		!x,y方向半轴等分段数，nx=ny
		real(dp)	:: hy		!y方向网格边长
		real(dp)	:: ux,uy,phx,phy   !无量纲蠕滑量，采用kalker量纲分析方法
		real(dp)	:: uh   !无量纲蠕滑量，采用kalker量纲分析方法
		real(dp)	:: spole(1,2)	!旋转中心  
		logical		:: hertz						!hertz压力分布控制字，true:hertz；false:non-hertz
		logical		:: heathcote				!heathcote slip
		real(dp),allocatable:: hx(:)										!x方向网格边长
		real(dp),allocatable:: area(:)									!单元面积
		real(dp),allocatable:: xc(:,:),yc(:,:)					!法向压力点坐标
		real(dp),allocatable:: tx(:,:),ty(:,:)	!无量钢拖动力数组
		real(dp),allocatable:: sx(:,:),sy(:,:)	!相对滑移速度，无量钢 relative velocity
		real(dp)	:: tfx,tfy
	end type fastdata

	!以下为轴承系统、滚动体、解除对数据结构
	type:: patchabk
		real(dp)	:: rx1,ry1,rx2,ry2
		real(dp)	:: rx,ry
		real(dp)	:: and,bnd,appnd !k为接触刚度系数 and,bnd 为不包括法向力的半轴长度，并非是无量纲量
		real(dp)	:: ei1,ei2
		real(dp)	:: ke   !椭圆率
		real(dp)	:: kp   !椭圆率计算精度
		real(dp)	:: e1
		real(dp)	:: e2
		real(dp)	:: nu1
		real(dp)	:: nu2
		logical		:: xlarger				!
	end type patchabk

	type comm  !共有参数
		real(dp),allocatable	:: and(:),bnd(:),appnd(:)  !无压力半轴、趋近量
		real(dp),allocatable	:: ke(:)
		real(dp),allocatable	:: ei1(:),ei2(:)
		real(dp),allocatable	:: rx(:),ry(:)	!综合曲率半径、曲率和
		real(dp),allocatable	:: k(:)	!接触刚度系数
		real(dp)	:: rx2,ry2			!滚道曲率半径
		logical		:: xlarger				!
		real(dp)	:: c11,c22,c23,c32,c33 !kalker系数
		REAL(DP)	:: RC		!Common Radius of Curvature 
	endtype comm

	type contact
		real(dp)	:: ax,by,app,ab ! 半轴、趋近量
		real(dp)	:: rdmd ! 接触中心变形后半径
		real(dp)	:: p		!接触对正压力
		real(dp)	:: alfst !静态接触角contact angle 包括受离心力的情况
		real(dp)	:: alf !contact angle
		real(dp)	:: xix,xiy,psi	!蠕滑率
		real(dp)	:: ux,uh,uy,phx,phy   !无量纲蠕滑量，采用kalker量纲分析方法
		real(dp)	:: tx,ty !滚动体受到的集中拖动力
		real(dp)	:: mx,my,mz !滚动体受到的力矩
		real(dp)	:: rd			!????????
		real(dp)	:: rbd		!滚动体变形后无量纲半径
		real(dp)	:: spn2rol  !spin to roll ratio
	endtype contact

	type roller
		type(contact) :: inn,out
		real(dp)	:: e	!滚动体节圆半径
		real(dp)	:: e0_ictl,e0_octl	!内、外控制无摩擦（理想）状态滚动体节圆半径
		real(dp)	:: est	!静平衡无摩擦状态滚动体节圆半径
		real(dp)	:: mgx,mgz !滚动体受到的陀螺力矩
		real(dp)	:: fc	!滚动体受到的离心力
		real(dp)	:: posang !
		real(dp)	:: omgc		!每个滚动体公转角速度
		real(dp)	:: omgb,omgx,omgy,omgz		!滚动体自转角速度
		real(dp),allocatable	:: f(:)  !方程组最终余量
	endtype roller
	
	type:: bearsysdata
		type(comm)	:: icom,ocom
		logical		:: simpmoment			!简化计算
 		logical		:: sym			!是否对称载荷
		logical		:: hertz						!hertz压力分布控制字，true:hertz；false:non-hertz
		logical		:: OPTIN
		logical		:: unscl	!unscale coordinates? TRUE:UNSCALE;FALSE:NO
		logical		:: rcm		!raceway control method
		logical		:: absolute		!绝对运动,所有角速度均相对于惯性坐标系而言
		logical		:: cmdcreep	!联合载荷作用下的蠕滑模型
		integer		:: s2rmthd	!旋滚比计算方法
		integer		:: hcmthd  !heathcote计算方法
		integer		:: z				!滚动体个数
		integer		:: cyc
		integer		:: inert	!0:no inert;1:fc only;2:fc+mgx;3:fc+mgz;4:fc+mgx+mgz
		real(dp),allocatable	:: kn(:)  !联合刚度系数，静力平衡时
		real(dp)	:: d,r			!diameter of balls直径、半径
		real(dp)	:: ri,ro		!inner and outer raceway groove curvature radius
		real(dp)	:: d_i,d_o	!inner and outer raceway diameter
		real(dp)	:: alf0			!free contact angle
		real(dp)	:: a,ai,ao  !distance between raceway groove curvature centers
		real(dp)	:: rx1,ry1	!滚动体曲率半径
		real(dp)	:: dm				!mean diameter
		real(dp)	:: fo,fi		!raceway curvature/ball dia
		real(dp)	:: pd				!bearing diametral clearance
		real(dp)	:: kp,kpst	!椭圆率、静态接触角计算精度
		real(dp)	:: e1,g1		!滚动体弹性、剪切模量
		real(dp)	:: e2,g2		!套圈弹性、剪切模量
		real(dp)	:: m1				!滚动体质量
		real(dp)	:: ip				!滚动体转动惯量
		real(dp)	:: mu				!摩擦系数
		real(dp)	:: nu1			!滚动体弹性模量
		real(dp)	:: nu2			!滚动体弹性模量
		real(dp)	:: rho1			!密度
		real(dp)	:: rho2			!密度
		real(dp)	:: np		!n prime  详见kalker fastsim算法文献

		real(dp)	:: aload		!轴承轴向载荷
		real(dp)	:: rload		!轴承径向载荷
		real(dp)	:: mload		!轴承径向载荷
		real(dp)	:: xtran=0.0_dp	!平易x坐标，便于tecplot，该值为椭圆主轴的倍数，
		type(roller),allocatable :: balno(:)
		real(dp),allocatable	:: x(:),x0(:)
		real(dp)	:: dlta,dltr,theta
		real(dp)	:: f3(3)
		real(dp)	:: dltast,dlta0 !静平衡纯轴向位移，不考虑摩擦拟静力学轴向位移
		real(dp)	:: dlta0_ictl,dlta0_octl !内、外控动平衡轴向位移
		real(dp)	:: s2rinn	!外滚道控制，内滚道旋滚比
		real(dp)	:: s2rout	!内滚道控制，外滚道旋滚比
		real(dp)	:: omgi		!内圈转速
		real(dp)	:: omgo		!外圈转速
		real(dp)	:: omgo_ictl,omgo_octl !内外孔外圈转速
		real(dp)	:: omgc		!保持架角速度
		real(dp)	:: omgb_ictl,omgb_octl	!滚动体转速初值
		real(dp)	:: iconang_ictl,oconang_ictl !内控接触角
		real(dp)	:: iconang_octl,oconang_octl !内控接触角
		real(dp)	:: startang !起始滚动体的方位角
		real(dp)	:: pchang !滚动体俯仰角
		real(dp)	:: pchang_octl	!外滚道控制俯仰角
		real(dp)	:: pchang_ictl	!内滚道控制俯仰角

		real(dp)	:: aj,gama	!传递参数
		real(dp)	:: tmp1,tmp2	!传递参数

		integer		:: m=6,n=6  !方程、未知数的个数
		
		integer		:: nsec
		real(dp),allocatable	:: alf(:)
		
		integer		::	index=0
		integer		::	intmthd	!1:接触角->变形量;2:变形量->接触角 忽略摩擦力
		integer		::	crpmthd	!1:接触角->变形量;2:变形量->接触角　考虑摩擦力
		logical		::	init
	end type bearsysdata

	type(fastdata)	 ,target	:: fs
	type(bearsysdata),target	:: br
end module fsbr

subroutine freefs()
	use fsbr
	implicit none
	logical				::saveslip

	deallocate(fs%hx,fs%area,fs%xc,fs%yc,fs%tx,fs%ty,fs%sx,fs%sy)
	deallocate(br%balno)
end subroutine freefs

subroutine freemem()
	use fsbr
	implicit none
	logical				::saveslip

	deallocate(fs%hx,fs%area,fs%xc,fs%yc,fs%tx,fs%ty,fs%sx,fs%sy)
	deallocate(br%balno)
end subroutine freemem
