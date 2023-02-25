	SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(DP), INTENT(IN) :: fold,stpmax
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x
	REAL(DP), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	INTERFACE
		FUNCTION func(x)
		USE NRTYPE
		IMPLICIT NONE
		REAL(DP) :: func
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(DP), PARAMETER :: ALF=1.0e-4_Dp,TOLX=epsilon(x)
	INTEGER(I4B) :: ndum
	REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	slope=dot_product(g,p)
	if (slope >= 0.0_DP) call nrerror('roundoff problem in lnsrch')
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_Dp))
	alam=1.0_DP
	do
		x(:)=xold(:)+alam*p(:)
		f=func(x)
		if (alam < alamin) then
			x(:)=xold(:)
			check=.true.
			RETURN
		else if (f <= fold+ALF*alam*slope) then
			RETURN
		else
			if (alam == 1.0_DP) then
				tmplam=-slope/(2.0_Dp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2_Dp)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2_Dp)/&
					(alam-alam2)
				if (a == 0.0_Dp) then
					tmplam=-slope/(2.0_Dp*b)
				else
					disc=b*b-3.0_Dp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_Dp*alam
					else if (b <= 0.0_Dp) then
						tmplam=(-b+sqrt(disc))/(3.0_Dp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_Dp*alam) tmplam=0.5_Dp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_Dp*alam)
	end do
	END SUBROUTINE lnsrch
