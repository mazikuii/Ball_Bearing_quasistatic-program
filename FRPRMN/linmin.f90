MODULE f1dim_mod
	USE nrtype
	INTEGER(I4B) :: ncom
	REAL(DP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
!BL
	FUNCTION f1dim(x)
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: f1dim
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP), DIMENSION(:), ALLOCATABLE :: xt
	allocate(xt(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	f1dim=func(xt)
	deallocate(xt)
	END FUNCTION f1dim
END MODULE f1dim_mod

	SUBROUTINE linmin(p,xi,fret)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : mnbrak,brent
	USE f1dim_mod
	IMPLICIT NONE
	REAL(DP), INTENT(OUT) :: fret
	REAL(DP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
	REAL(DP), PARAMETER :: TOL=1.0e-4_Dp
	REAL(DP) :: ax,bx,fa,fb,fx,xmin,xx
	ncom=assert_eq(size(p),size(xi),'linmin')
	pcom=>p
	xicom=>xi
	ax=0.0_Dp
	xx=1.0_Dp
	call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=brent(ax,xx,bx,f1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
	END SUBROUTINE linmin
