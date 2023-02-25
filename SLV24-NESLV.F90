subroutine slv24_neslv()
	use	nrtype
	use neintface
	USE FSBR
	implicit none
		
	real(dp),allocatable::x(:),fp(:),gp(:),typx(:),typy(:),fvec(:)
	logical::analjac,restart,maxtaken
	integer::lshmsd,global,method,termcode
	integer::itncount,itnlimit,itconsec,consecmax,nmt,nfnt,njnt
	real(dp)::ftol,xtol,mtol,btol,eta0,stepmx,dlt

	integer	:: n
	integer::i,mxit
	logical(lgt)::pnt

	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************    BY NESLV METHOD     ************'
	WRITE(*,*)'            ***************************************************'
	
	n=br%n
	allocate(x(n),typx(n),typy(n),fp(n),gp(n),fvec(n))
	call nedflt(analjac,restart,maxtaken,&
							lshmsd,global,method,termcode,&
							itncount,itnlimit,itconsec,consecmax,nmt,&
							ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
							x,typx,typy)
	ftol			=1.e-16_dp
	xtol			=1.e-16_dp
	mtol			=1.e-16_dp
	btol			=1.e-16_dp
	stepmx		=1.0_dp
!	typy			=1.0e-17_dp
	typx			=1.0e-0_dp
	method		=2	!1-newton;2-broyden
	global		=1	!1-line search;2-double dogleg;3-simplified dogleg;4-dogleg mzk
	lshmsd		=3	!1-mono;2-averged;3-max;4-convex
	eta0			=.7_dp
	nmt				=6
	dlt				=1.0_dp

	consecmax	=120
	itnlimit	=200
	pnt=.false.
	
	x(:)=0.0_dp
	call x0vct(x)
	call printx0(x)
	
	mxit=7
	do i=1,mxit
		if(i>5) then
!			ftol			=ftol*1_dp
!			xtol			=ftol
!			mtol			=ftolhttps://securelogin.arubanetworks.com/cgi-bin/login?cmd=logout
!			btol			=ftol
!			typx=x
!			br%mu=br%mu*.995_dp
!			nmt=nmt-2
!			typy			=fvec
!	global		=2
!			eta0=.0_dp
		endif
if(i==7) global=2  !计算精度最高
		call neslv(x,fp,gp,fv24,jac24,typx,typy,&
							 analjac,restart,maxtaken,&
							 lshmsd,global,method,termcode,&
							 itncount,itnlimit,itconsec,consecmax,nmt,&
							 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
							 nfnt,njnt)
		write(*,"(a7,i4,5x,a9,i3)") 'global:',global,'termcode:',termcode
		write(*,"(a5,i4,5x,a5,i3)") 'nfnt:',nfnt,'njnt:',njnt
		call savexf(x,fp)
		if(n==6) then
			call slv6eqs(x,fvec,.true.)  !重新计算，保存必要的中间变量 
		elseif(n==8) then
			call slv8eqs(x,fvec,.true.)  !重新计算，保存必要的中间变量 
		endif
		call printxf(x,fvec)
		call savexf(x,fvec)
		if(pnt) call neexpln(termcode,consecmax)
	enddo
	deallocate(x,typx,typy,fp)
endsubroutine slv24_neslv

function fv24(x)
	use nrtype
	use neintface,only:slv6eqs,slv8eqs
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	fv24(size(x(:)))
	real(dp)::	fvec(size(x(:)))

	if(size(x)==6)then
		call slv6eqs(x,fvec,.false.)
	elseif(size(x)==8) then
		call slv8eqs(x,fvec,.false.)
	endif
	call dcopy(size(x),fvec,1,fv24,1)
endfunction fv24	

subroutine jac24(x,jc)
	use nrtype
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp),intent(out)::	jc(:,:)
	
	jc(1,1)=1.0_dp
	jc(1,2)=-1.0_dp
	jc(2,1)=1.0_dp
	jc(2,2)=1.0_dp
end subroutine jac24	
