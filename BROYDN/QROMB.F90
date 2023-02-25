	FUNCTION qromb(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP) :: qromb,KE
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(DP), PARAMETER :: EPS=1.0e-8_dp
	REAL(DP), DIMENSION(JMAXP) :: h,s
	REAL(DP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_dp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
