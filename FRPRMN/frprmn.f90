	SUBROUTINE frprmn(p,ftol,iter,fret)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : linmin
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(IN) :: ftol
	REAL(DP), INTENT(OUT) :: fret
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(p)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: p
		REAL(DP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(p)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: p
		REAL(DP), DIMENSION(size(p)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=2000
	REAL(DP), PARAMETER :: EPS=1.0e-10_Dp
	INTEGER(I4B) :: its
	REAL(DP) :: dgg,fp,gam,gg
	REAL(DP), DIMENSION(size(p)) :: g,h,xi
	fp=func(p)
	xi=dfunc(p)
	g=-xi
	h=g
	xi=h
	do its=1,ITMAX
		iter=its
		call linmin(p,xi,fret)
		if (2.0_Dp*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
		fp=fret
		xi=dfunc(p)
		gg=dot_product(g,g)
!		dgg=dot_product(xi,xi)
		dgg=dot_product(xi+g,xi)
		if (gg == 0.0) RETURN
		gam=dgg/gg
		g=-xi
		h=g+gam*h
		xi=h
	end do
	call nrerror('frprmn: maximum iterations exceeded')
	END SUBROUTINE frprmn
