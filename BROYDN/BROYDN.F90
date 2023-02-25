	SUBROUTINE broydn(x,check)
	USE FSBR	
	USE nrtype; USE nrutil, ONLY : get_diag,lower_triangle,nrerror,&
		outerprod,put_diag,unit_matrix,vabs
	USE nr, ONLY : fdjac,lnsrch,qrdcmp,qrupdt,rsolv
	USE fminln
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=200
	REAL(DP), PARAMETER ::	EPS=epsilon(x),				&
													TOLF=1.0e-12_Dp,			&
													TOLMIN=1.0e-12_Dp,		&
													TOLX=1.0e-12_Dp,			&
													STPMX=100.0000000_DP
	INTEGER(I4B) :: i,its,k,n
	REAL(DP) :: f,fold,stpmax
	REAL(DP), DIMENSION(size(x)), TARGET :: fvec
	REAL(DP), DIMENSION(size(x)) :: c,d,fvcold,g,p,s,t,w,xold
	REAL(DP), DIMENSION(size(x),size(x)) :: qt,r
	LOGICAL :: restrt,sing

	fmin_fvecp=>fvec
	n=size(x)
	f=fmin(x)
	if (maxval(abs(fvec(:))) < 0.01_Dp*TOLF) then
		check=.false.
		RETURN
	end if
	stpmax=STPMX*max(vabs(x(:)),real(n,Dp))
	restrt=.true.	
	do its=1,MAXITS
		if (restrt) then
			call fdjac(x,fvec,r)
			call qrdcmp(r,c,d,sing)
			if (sing) call nrerror('singular Jacobian in broydn')
			call unit_matrix(qt)
			do k=1,n-1
				if (c(k) /= 0.0_DP) then
					qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),&
						matmul(r(k:n,k),qt(k:n,:)))/c(k)
				end if
			end do
			where (lower_triangle(n,n)) r(:,:)=0.0_DP
			call put_diag(d(:),r(:,:))
		else
			s(:)=x(:)-xold(:)
			do i=1,n
				t(i)=dot_product(r(i,i:n),s(i:n))
			end do
			w(:)=fvec(:)-fvcold(:)-matmul(t(:),qt(:,:))
			where (abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) &
				w(:)=0.0_DP
			if (any(w(:) /= 0.0_DP)) then
				t(:)=matmul(qt(:,:),w(:))
				s(:)=s(:)/dot_product(s,s)
				call qrupdt(r,qt,t,s)
				d(:)=get_diag(r(:,:))
				if (any(d(:) == 0.0_DP)) &
					call nrerror('r singular in broydn')
			end if
		end if
		p(:)=-matmul(qt(:,:),fvec(:))
		do i=1,n
			g(i)=-dot_product(r(1:i,i),p(1:i))
		end do
		xold(:)=x(:)
		fvcold(:)=fvec(:)
		fold=f
		call rsolv(r,d,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
		if (maxval(abs(fvec(:))) < TOLF) then
			check=.false.
			RETURN
		end if
		if (check) then
			if (restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
				1.0_Dp)/max(f,0.5_Dp*n)) < TOLMIN) RETURN
			restrt=.true.
		else
			restrt=.false.
			if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
				1.0_Dp)) < TOLX) RETURN
		end if
	end do
	MAXIT=.TRUE.
	call nrerror('MAXITS exceeded in broydn')
	END SUBROUTINE broydn