subroutine nedflt(analjac,restart,maxtaken,&
									lshmsd,global,method,termcode,&
									itncount,itnlimit,itconsec,consecmax,nmt,&
									ftol,xtol,mtol,btol,etak,stepmx,dlt,&
									x0,sx,sy)
	use nrtype
	implicit none
	logical::analjac,restart,maxtaken,nmtlsh
	integer::lshmsd,global,method,termcode
	integer::itncount,itnlimit,itconsec,consecmax,nmt
	real(dp)::ftol,xtol,mtol,btol,etak,stepmx,dlt
	real(dp)::x0(:),sx(:),sy(:)
		
	!defaul setting: 
	!								factored broyden update
	!								line search global strategy
	analjac		=.false.
	restart		=.true.
	maxtaken	=.false.

	!1-   monotone line search
	!2-nonmonotone line search by limited avaeraged
	!3-nonmonotone line search by limited maximum 
	!4-nonmononone line search by convex combination
	lshmsd		=1     
	global			=1
	method			=2
	termcode		=0

	itncount		=0
	itnlimit		=200
	itconsec		=0
	consecmax		=5
	nmt					=5

	ftol			=1.0e-8_dp
	xtol			=1.0e-8_dp
	mtol		=1.0e-8_dp
	btol			=1.0e-8_dp
	etak			=.85_dp
	stepmx			=1.0_dp
	dlt					=1.0_dp

!	x0(:)=8.8888888888888_dp
	sx(:)=1.0_dp
	sy(:)=1.0_dp
endsubroutine nedflt