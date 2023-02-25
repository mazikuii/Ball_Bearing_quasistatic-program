	SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
	USE nrtype
	USE nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(DP), INTENT(IN) :: tolx,tolf
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
		END SUBROUTINE usrfun
	END INTERFACE
	INTEGER(I4B) :: i
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(DP) :: d
	REAL(DP), DIMENSION(size(x)) :: fvec,p
	REAL(DP), DIMENSION(size(x),size(x)) :: fjac
	do  i=1,ntrial
		call usrfun(x,fvec,fjac)
		if (sum(abs(fvec)) <= tolf) RETURN
		p=-fvec
		call ludcmp(fjac,indx,d)
		call lubksb(fjac,indx,p)
		x=x+p
		if (sum(abs(p)) <= tolx) RETURN
	end do
	END SUBROUTINE mnewt
