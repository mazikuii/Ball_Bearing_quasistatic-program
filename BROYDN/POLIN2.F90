	SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : polint
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: ya
	REAL(DP), INTENT(IN) :: x1,x2
	REAL(DP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: j,m,ndum
	REAL(DP), DIMENSION(size(x1a)) :: ymtmp
	REAL(DP), DIMENSION(size(x2a)) :: yntmp
	m=assert_eq(size(x1a),size(ya,1),'polin2: m')
	ndum=assert_eq(size(x2a),size(ya,2),'polin2: ndum')
	do j=1,m
		yntmp=ya(j,:)
		call polint(x2a,yntmp,x2,ymtmp(j),dy)
	end do
	call polint(x1a,ymtmp,x1,y,dy)
	END SUBROUTINE polin2
