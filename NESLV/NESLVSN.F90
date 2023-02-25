subroutine neslvsn(qt,r,fvec,p,g,qtf,sx,sy,macheps)
	use nrtype; use nrutil, only : get_diag,vabs
	use nr, only : rsolv
	implicit none
	real(dp), dimension(:,:) :: qt,r
	real(dp), dimension(:)	 :: fvec,p,g,qtf,sx,sy
	real(dp)								 :: macheps

	integer(i4b) :: n,i,ierr,j
	real(dp), dimension(size(qt(1,:)))::  d,y,rcdwrk,icdwrk
	real(dp)	:: rcond,mu
		
	real(dp), dimension(size(qt(1,:)))::  wa,w

	logical	:: dogleg
	real(dp):: stpmx

	n=size(qt(1,:))

	dogleg=.true.
	if(dogleg) then
		stpmx=1333._dp
		qtf(:)=-matmul(qt(:,:),sy(:)*fvec(:))  !注意与下面代码的符号差别
		p(:)=-qtf(:)
		do i=1,n
			g(i)=dot_product(r(1:i,i),p(1:i))
		end do
		!可以替代QR分解方法求解牛顿方程
		call nedogleg(n,r,sx,qtf,stpmx,p)  !dogleg from minpack for newton direction
	else
		ierr=0
		call NEcndjac(n,r,n,macheps,rcond,y,rcdwrk,icdwrk,ierr,mu) 
		qtf(:)=matmul(qt(:,:),sy(:)*fvec(:))
		if(ierr==0) then
			p(:)=-qtf(:)
			do i=1,n
				g(i)=-dot_product(r(1:i,i),p(1:i))
			end do
			d=get_diag(r(:,:))
			call rsolv(r,d,p)
		else
			call dcopy(n,sx,1,wa,1)
			call dscal(n,mu,wa,1)
			call liqrev(n,r,n,wa,qtf,p,y,w)
		endif
	endif
end subroutine neslvsn