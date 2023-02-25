SUBROUTINE TT()
	use	nrtype
	use neintface
	implicit none
		
	integer,parameter:: n=2
	real(dp),allocatable::x(:),fp(:),gp(:),typx(:),typy(:)
	logical::analjac,restart,maxtaken
	integer::lshmsd,global,method,termcode
	integer::itncount,itnlimit,itconsec,consecmax,nmt,nfnt,njnt
	real(dp)::ftol,xtol,mtol,btol,etak,stepmx,dlt,macheps

	integer::i,mxit
	real(dp):: fvec(N)


	WRITE(*,*)'            ***************************************************'
	WRITE(*,*)'            ***************    BY NESLV METHOD     ************'
	WRITE(*,*)'            ***************************************************'
	

	allocate(x(n),typx(n),typy(n),fp(n),gp(n))
	call nedflt(analjac,restart,maxtaken,&
							lshmsd,global,method,termcode,&
							itncount,itnlimit,itconsec,consecmax,nmt,&
							ftol,xtol,mtol,btol,etak,stepmx,dlt,&
							x,typx,typy)
	xtol=1.0d-16
	ftol=1.0d-16
	btol=1.0d-16
	mtol=1.0d-6
	stepmx=100.0d0
	analjac=.false.
	global=1
	dlt=-1.0_dp
	consecmax=180
	itnlimit=200
  
  x(1)  = 01.0000d-5
  x(2)  = 7.0d0
	mxit=3
	do i=1,mxit
		call neslv(x,fp,gp,fvtt,jactt,typx,typy,&
							 analjac,restart,maxtaken,&
							 lshmsd,global,method,termcode,&
							 itncount,itnlimit,itconsec,consecmax,nmt,&
							 ftol,xtol,mtol,btol,etak,stepmx,dlt,&
							 nfnt,njnt)
		write(*,"(15x,a9,i3)") 'termcode:',termcode
		write(*,*) 'nfnt:',nfnt,'njnt:',njnt
		write(*,*) 'x :',x
		write(*,*) 'fp:',fp
		write(*,*) 'gp:',gp
	enddo
	
	deallocate(x,typx,typy,fp)
ENDSUBROUTINE TT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function fvtt(x)
	use nrtype
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp)::	fvtt(size(x(:)))
	real(dp)::	fvec
	
	fvtt(1) = 1.0D4 * X(2)* X(1)-1.0D0
	fvtt(2) = EXP(-X(1))+EXP(-X(2))-1.001D0
endfunction fvtt	

subroutine jactt(x,jc)
	use nrtype
	implicit none
	real(dp),intent(in)	::	x(:)
	real(dp),intent(out)::	jc(:,:)

	jc(1,1)=1.0D4 * X(2)
	jc(1,2)=1.0D4 * X(1)
	jc(2,1)=-X(1)*EXP(-X(1))
	jc(2,2)=-X(2)*EXP(-X(2))
end subroutine jactt	
