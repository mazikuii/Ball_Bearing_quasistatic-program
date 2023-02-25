subroutine neqrudt(qt,r,x,xold,f,fold,fvec,fvcold,sx,sy,stpmax,macheps)
	use nrtype; use nrutil, only : nrerror,get_diag
	use nr, only : qrupdt
	implicit none
	real(dp), dimension(:,:) :: qt,r
	real(dp), dimension(:) :: x,xold
	real(dp)::f,fold
	real(dp), dimension(:) :: fvec,fvcold,sx,sy
	real(dp)::stpmax,macheps

	INTEGER :: i,n
	REAL(DP), DIMENSION(size(x)) :: d,s,t,w

	n=size(x)
	s(:)=x(:)-xold(:)
	do i=1,n
		t(i)=dot_product(r(i,i:n),s(i:n))
	end do
	w(:)=sy(:)*(fvec(:)-fvcold(:))-matmul(t(:),qt(:,:))
	where (abs(w(:)) < macheps*sy(:)*(abs(fvec(:))+abs(fvcold(:)))) &
		w(:)=0.0_DP
	if (any(w(:) /= 0.0_DP)) then
		t(:)=matmul(qt(:,:),w(:))
		s(:)=sx(:)*sx(:)*s(:)/dot_product(sx(:)*s,sx(:)*s)
		call qrupdt(r,qt,t,s)
		d(:)=get_diag(r(:,:))
!		if (any(d(:) == 0.0_DP)) &			!真奇怪，去掉这条语句，精度反而更高！
!			call nrerror('r singular in neqrudt')
	end if
endsubroutine neqrudt