	SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc)
	USE nrtype; USE nrutil, ONLY : nrerror,outerprod,unit_matrix,vabs
	USE nr, ONLY : lnsrch
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(IN) :: gtol
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
	INTEGER(I4B), PARAMETER :: ITMAX=200
	REAL(DP), PARAMETER :: STPMX=100.0_dp,EPS=epsilon(p),TOLX=4.0_dp*EPS
	INTEGER(I4B) :: its
	LOGICAL :: check
	REAL(DP) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
	REAL(dP), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
	REAL(DP), DIMENSION(size(p),size(p)) :: hessin
	fp=func(p)
	g=dfunc(p)
	call unit_matrix(hessin)
	xi=-g
	stpmax=STPMX*max(vabs(p),real(size(p),sp))
	do its=1,ITMAX
		iter=its
		call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
		fp=fret
		xi=pnew-p
		p=pnew
		if (maxval(abs(xi)/max(abs(p),1.0_dp)) < TOLX) RETURN
		dg=g
		g=dfunc(p)
		den=max(fret,1.0_Dp)
		if (maxval(abs(g)*max(abs(p),1.0_dp)/den) < gtol) RETURN
		dg=g-dg
		hdg=matmul(hessin,dg)
		fac=dot_product(dg,xi)
		fae=dot_product(dg,hdg)
		sumdg=dot_product(dg,dg)
		sumxi=dot_product(xi,xi)
		if (fac > sqrt(EPS*sumdg*sumxi)) then
			fac=1.0_dp/fac
			fad=1.0_dp/fae
			dg=fac*xi-fad*hdg
			hessin=hessin+fac*outerprod(xi,xi)-&
				fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
		end if
		xi=-matmul(hessin,g)
	end do
	call nrerror('dfpmin: too many iterations')
	END SUBROUTINE dfpmin
