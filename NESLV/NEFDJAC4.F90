subroutine nefdjac4(x,func,fvec,df,sx)
	use nrtype; use nrutil, only : assert_eq
	implicit none
	real(dp), dimension(:), intent(inout) :: x
	real(dp), dimension(size(x)), intent(in) :: fvec,sx
	real(dp), dimension(size(x),size(x)), intent(out) :: df
	interface
		function func(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	func(size(x(:)))
		endfunction func	
	end interface

	real(dp) :: eps
	integer :: j,n
	real(dp), dimension(size(x)) :: xsav,xph,h,fv
	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'nefdjac')
	j=1
!	eps=1.0e-4_dp	!效果很差，很难收敛
	eps=sqrt(epsilon(0.0_dp))		!效果很好
!	eps=1.0e-12_dp	!效果很好
!	eps=epsilon(x)

	call dcopy(n,x(:),1,xsav(:),1)
	h(:)=eps*max(abs(xsav(:)),1.0_dp/sx(:))
	where (h(:) == 0.0_dp) h(:)=eps
	xph(:)=xsav(:)+h(:)
	h=xph-xsav
	do j=1,n
		x(j)=xph(j)
		fv=func(x)
		df(:,j)=(fv-fvec(:))/h(j)
		x(j)=xsav(j)
	end do

end subroutine nefdjac4