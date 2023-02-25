subroutine getpchaxis(alfi,alfo)
	!计算静态接触角
	use fsbr
	use nr
	use intface,only:pow
	implicit none
	integer	 :: i
	real(dp)	::alfi,alfo
	real(dp)	:: p13,p23
	real(dp)	::and,bnd,appnd,dy
	
	!在给定载荷下，计算接触椭圆半轴、趋近量
	!内接触
	p13=pow(br%balno(1)%inn%p,1.0_dp/3.0_dp)
	p23=pow(br%balno(1)%inn%p,2.0_dp/3.0_dp)
	call ratint(br%alf,br%icom%and	,alfi,and,dy)
	call ratint(br%alf,br%icom%bnd	,alfi,bnd,dy)
	call ratint(br%alf,br%icom%appnd,alfi,appnd,dy)
	br%balno(:)%inn%ax  =and	*p13
	br%balno(:)%inn%by  =bnd	*p13
	br%balno(:)%inn%ab  =br%balno(1)%inn%ax*br%balno(1)%inn%by
	br%balno(:)%inn%app =appnd	*p23
	!外接触
	p13=pow(br%balno(1)%out%p,1.0_dp/3.0_dp)
	p23=pow(br%balno(1)%out%p,2.0_dp/3.0_dp)
	call ratint(br%alf,br%ocom%and	,alfo,and,dy)
	call ratint(br%alf,br%ocom%bnd	,alfo,bnd,dy)
	call ratint(br%alf,br%ocom%appnd,alfo,appnd,dy)
	br%balno(:)%out%ax  =and		*p13
	br%balno(:)%out%by  =bnd		*p13
	br%balno(:)%out%ab  =br%balno(1)%out%ax*br%balno(1)%out%by
	br%balno(:)%out%app =appnd		*p23

	select case(br%hcmthd)
	case(1:2) !johnson & johnson+mzk
		br%balno(:)%inn%rd	=br%r/br%balno(1)%inn%by				!rc/b  axis b of contact patch 用来计算herthcote slip
		br%balno(:)%out%rd	=br%r/br%balno(1)%out%by				!rc/b  axis b of contact patch 用来计算herthcote slip
	case(3) !mzk
		br%balno(:)%inn%rd	=br%icom%rc/br%balno(1)%inn%by				!rc/b  axis b of contact patch 用来计算herthcote slip
		br%balno(:)%out%rd	=br%ocom%rc/br%balno(1)%out%by				!rc/b  axis b of contact patch 用来计算herthcote slip
		br%balno(:)%inn%rbd	=(br%r-.5_dp*br%balno(1)%inn%app)/br%balno(1)%inn%by
		br%balno(:)%out%rbd	=(br%r-.5_dp*br%balno(1)%out%app)/br%balno(1)%out%by
	endselect
endsubroutine getpchaxis

subroutine getpchaxis_idx(idx,alfi,alfo)
	!计算静态接触角
	use fsbr
	use nr
	use intface,only:pow
	implicit none
	integer	 :: idx
	real(dp)	::alfi,alfo
	real(dp)	:: p13,p23
	real(dp)	::and,bnd,appnd,dy
	
	!在给定载荷下，计算接触椭圆半轴、趋近量
	!内接触
	p13=pow(br%balno(idx)%inn%p,1.0_dp/3.0_dp)
	p23=pow(br%balno(idx)%inn%p,2.0_dp/3.0_dp)
	call ratint(br%alf,br%icom%and	,alfi,and		,dy)
	call ratint(br%alf,br%icom%bnd	,alfi,bnd		,dy)
	call ratint(br%alf,br%icom%appnd,alfi,appnd	,dy)
	br%balno(idx)%inn%ax  =and	*p13
	br%balno(idx)%inn%by  =bnd	*p13
	br%balno(idx)%inn%ab  =br%balno(idx)%inn%ax*br%balno(idx)%inn%by
	br%balno(idx)%inn%app =appnd	*p23
	!外接触
	p13=pow(br%balno(idx)%out%p,1.0_dp/3.0_dp)
	p23=pow(br%balno(idx)%out%p,2.0_dp/3.0_dp)
	call ratint(br%alf,br%ocom%and	,alfo,and		,dy)
	call ratint(br%alf,br%ocom%bnd	,alfo,bnd		,dy)
	call ratint(br%alf,br%ocom%appnd,alfo,appnd	,dy)
	br%balno(idx)%out%ax  =and		*p13
	br%balno(idx)%out%by  =bnd		*p13
	br%balno(idx)%out%ab  =br%balno(idx)%out%ax*br%balno(idx)%out%by
	br%balno(idx)%out%app =appnd		*p23

	select case(br%hcmthd)
	case(1:2) !johnson & johnson+mzk
		br%balno(idx)%inn%rd	=br%r/br%balno(idx)%inn%by				!rc/b  axis b of contact patch 用来计算herthcote slip
		br%balno(idx)%out%rd	=br%r/br%balno(idx)%out%by				!rc/b  axis b of contact patch 用来计算herthcote slip
	case(3) !mzk
		br%balno(idx)%inn%rd	=br%icom%rc/br%balno(idx)%inn%by				!rc/b  axis b of contact patch 用来计算herthcote slip
		br%balno(idx)%out%rd	=br%ocom%rc/br%balno(idx)%out%by				!rc/b  axis b of contact patch 用来计算herthcote slip
		br%balno(idx)%inn%rbd	=(br%r-.5_dp*br%balno(idx)%inn%app)/br%balno(idx)%inn%by
		br%balno(idx)%out%rbd	=(br%r-.5_dp*br%balno(idx)%out%app)/br%balno(idx)%out%by
	endselect
endsubroutine getpchaxis_idx