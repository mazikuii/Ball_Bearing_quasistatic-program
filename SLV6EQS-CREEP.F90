subroutine slv6eqs_creep(x)
	use	nrtype
	use neintface
	USE FSBR
	implicit none
		
	integer,parameter:: n=6
	real(dp)::x(:)
	real(dp),allocatable::fp(:),gp(:),typx(:),typy(:)
	logical::analjac,restart,maxtaken
	integer::lshmsd,global,method,termcode
	integer::itncount,itnlimit,itconsec,consecmax,nmt,nfnt,njnt
	real(dp)::ftol,xtol,mtol,btol,eta0,stepmx,dlt

	integer::i
	real(dp):: fvec(N)
	logical(lgt)::pnt

	allocate(typx(n),typy(n),fp(n),gp(n))
	call nedflt(analjac,restart,maxtaken,&
							lshmsd,global,method,termcode,&
							itncount,itnlimit,itconsec,consecmax,nmt,&
							ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
							x,typx,typy)
	ftol			=1.e-16_dp
	xtol			=1.e-16_dp
	mtol			=1.e-16_dp
	btol			=1.e-16_dp
	typy			=1.0e-0_dp
	typx			=1.0e-0_dp
	stepmx		=001.000_dp
	method		=2	!1-newton;2-broyden
	nmt				=6
	dlt				=.10000_dp

	global		=1	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=3	!1-mono;2-averged;3-max;4-convex
	eta0			=1._dp

	typy(1)	=1.0e-12_dp
	typy(2)	=1.0e-11_dp		!11
	typy(3)	=1.0e-11_dp		!11
	typy(4)	=1.0e-14_dp
	typy(5)	=1.0e-14_dp
	typy(6)	=1.0e-14_dp

	if(br%crpmthd==1) then
		typx(5)	=1.0_dp
		typx(6)	=1.0_dp
	elseif(br%crpmthd==2) then
		typx(5)	=1.0e-10_dp		!这里是怎么回事，到底该不该要，取值多少．
		typx(6)	=1.0e-10_dp
	endif

	consecmax	=120
	itnlimit	=200
	pnt=.false.
	
	i=0
	do while(i==0)
		call neslv(x,fp,gp,fv25,jac24,typx,typy,&
							 analjac,restart,maxtaken,&
							 lshmsd,global,method,termcode,&
							 itncount,itnlimit,itconsec,consecmax,nmt,&
							 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
							 nfnt,njnt)
		call slv6eqs_idx(x,fvec,.true.,br%index,br%crpmthd)  !重新计算，保存必要的中间变量
		if(sqrt(dot_product(fvec,fvec))<=1.0e-9_dp) i=1

	enddo
	global=2
	call neslv(x,fp,gp,fv25,jac24,typx,typy,&
						 analjac,restart,maxtaken,&
						 lshmsd,global,method,termcode,&
						 itncount,itnlimit,itconsec,consecmax,nmt,&
						 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
						 nfnt,njnt)
	call slv6eqs_idx(x,fvec,.true.,br%index,br%crpmthd)  !重新计算，保存必要的中间变量
	deallocate(typx,typy,fp)
endsubroutine slv6eqs_creep

function fv25(x)
	use nrtype
	use fsbr
	use neintface,only:slv6eqs_idx
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	fv25(size(x(:)))
	real(dp)::	fvec(size(x(:)))

	call slv6eqs_idx(x,fvec,.false.,br%index,br%crpmthd)
	call dcopy(size(x),fvec,1,fv25,1)
endfunction fv25