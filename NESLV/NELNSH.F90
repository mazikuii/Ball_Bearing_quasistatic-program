subroutine nelnsh(xold,fold,g,p,x,f,fvec,stpmax,btol,sx,sy,&
									maxtaken,lshmsd,eta0,nmt,retcode,fmin,func,nfnt)
	use nrtype
	use nrutil, only : assert_eq,nrerror,vabs
	use neintface,only:svfm
	implicit none
	real(dp):: xold(:),fold,g(:),p(:),x(:),f,fvec(:),stpmax,btol,sx(:),sy(:)
	logical	:: maxtaken
	integer :: lshmsd
	real(dp):: eta0
	integer	:: nmt,retcode,nfnt
	interface
		function fmin(x,func,fvec,sy)
			use nrtype
			implicit none
			real(dp),dimension(:), intent(in) :: x,sy
			real(dp),dimension(:)::fvec
			real(dp)::fmin
			interface
				function func(x)
					use nrtype
					implicit none
					real(dp),intent(in)	::	x(:)
					real(dp)::	func(size(x(:)))
				endfunction func	
			end interface
		endfunction fmin
		function func(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	func(size(x(:)))
		endfunction func	
	end interface

	real(dp), parameter :: alf=1.0e-4_dp
	integer :: ndum
	real(dp) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	real(dp):: fm(nmt),qkp
	integer ::it

	it=1
	fm(:)=0.0_sp
	fm(1)=fold
	qkp=1.0_dp
	retcode=2
	maxtaken=.false.
	ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	pabs=vabs(sx(:)*p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	pabs=stpmax
	slope=dot_product(g,p)
	if (slope >= 0.0) then
		if(lshmsd==1) then
			call nrerror('roundoff problem in lnsrch')		
		else
			p(:)=-p(:)
			slope=-slope
		endif
	endif
	alamin=btol/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp/sx(:)))
	alam=1.0_dp
	do
		x(:)=xold(:)+alam*p(:)
		f=fmin(x,func,fvec,sy)
		nfnt=nfnt+1
		if(lshmsd==1) then
			!monotone line search
			fold=f			!not changed
		elseif(lshmsd==2) then
			!nonmonotone line search by limited averaged 
			call svfm(f,fm)
			it=it+1
			fold=sum(fm)/it	!mean value
		elseif(lshmsd==3) then
			!nonmonotone line search by limited maximum
			call svfm(f,fm)
			fold=maxval(fm)	!by grippo
		else !(lshmsd==4)
			!nonmonotone line search by convex combination
			fold=(eta0*qkp*fold+f)/(eta0*qkp+1.0_dp)
			qkp=eta0*qkp+1.0_dp
		endif
		if (alam < alamin) then
			x(:)=xold(:)
			retcode=1
			return
		else if (f <= fold+alf*alam*slope) then
			retcode=0
			if(alam==1.0_dp .and. pabs>0.99_dp*stpmax) maxtaken=.true. 
			return
		else
			if (alam == 1.0_dp) then
				tmplam=-slope/(2.0_dp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2_dp)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2_dp)/&
					(alam-alam2)
				if (a == 0.0_dp) then
					tmplam=-slope/(2.0_dp*b)
				else
					disc=b*b-3.0_dp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_dp*alam
					else if (b <= 0.0_dp) then
						tmplam=(-b+sqrt(disc))/(3.0_dp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_dp*alam)
	end do
end subroutine nelnsh

subroutine svfm(f,fm)
	use nrtype
	implicit none
	real(dp)::f,fm(:)

	integer i,n

	n=size(fm)
	do i=n,2,-1
		fm(i)=fm(i-1)
	enddo
	fm(i)=f
end subroutine
