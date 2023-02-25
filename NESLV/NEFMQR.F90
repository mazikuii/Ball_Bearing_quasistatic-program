subroutine nefmqr(qt,r)
	use nrtype
	use nrutil, only : nrerror,put_diag,unit_matrix, &
		outerprod,lower_triangle
	use nr, only : qrdcmp
	implicit none
	real(dp), dimension(:,:) :: qt,r

	real(dp), dimension(size(r(1,:)))	 :: c,d
	logical		:: sing
	integer :: k,n

	n=size(r(1,:))
	call qrdcmp(r,c,d,sing)
	if (sing) call nrerror('singular jacobian in broydn')
	call unit_matrix(qt)
	do k=1,n-1
		if (c(k) /= 0.0_dp) then
			qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),&
				matmul(r(k:n,k),qt(k:n,:)))/c(k)
		end if
	end do
	where (lower_triangle(n,n)) r(:,:)=0.0_dp
	call put_diag(d(:),r(:,:))
endsubroutine nefmqr