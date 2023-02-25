SUBROUTINE GETCREEP2()
	USE FSBR
	IMPLICIT NONE
	INTEGER		:: I
	REAL(DP)	:: P,AX,BY,XIX,XIY,PSI,RD
	REAL(DP)	:: TEMP,AB

		!外接触
		P =BR%BALNO(1)%OUT%P
		AX=BR%BALNO(1)%OUT%AX
		BY=BR%BALNO(1)%OUT%BY
		AB=BR%BALNO(1)%OUT%AB
		XIX=BR%BALNO(1)%OUT%XIX
		XIY=BR%BALNO(1)%OUT%XIY
		PSI=BR%BALNO(1)%OUT%PSI

		BR%BALNO(:)%OUT%UX=0.375_DP*AB*BR%NP*BR%G1*BR%OCOM%C11*XIX/BR%MU/P
		BR%BALNO(:)%OUT%UY=0.375_DP*AB*BR%NP*BR%G1*BR%OCOM%C22*XIY/BR%MU/P
		BR%BALNO(:)%OUT%PHX=BY**2*SQRT(AB)*BR%G1*BR%OCOM%C23*PSI/BR%MU/P
		BR%BALNO(:)%OUT%PHY=   AB*SQRT(AB)*BR%G1*BR%OCOM%C23*PSI/BR%MU/P

ENDSUBROUTINE GETCREEP2

SUBROUTINE GETCREEP3(inn21,out21)
	USE FSBR
	IMPLICIT NONE
	real(dp),intent(in)::inn21,out21
	INTEGER		:: I
	REAL(DP)	:: P,AX,BY,Z0,XIX,XIY,PSI,RD
	REAL(DP)	:: TEMP,AB
	REAL(DP)	:: L1,L2,L3,C11,C22,C23,G,MU
	real(dp)	:: rbd

		!外接触
		P  =BR%BALNO(1)%OUT%P
		AX =BR%BALNO(1)%OUT%AX
		BY =BR%BALNO(1)%OUT%BY
		AB =BR%BALNO(1)%OUT%AB
		XIX=BR%BALNO(1)%OUT%XIX
		XIY=BR%BALNO(1)%OUT%XIY
		PSI=BR%BALNO(1)%OUT%PSI

		Z0=P/AB/BR%NP
		C11=BR%OCOM%C11
		C22=BR%OCOM%C22
		C23=BR%OCOM%C23
		G	 =BR%G1
		MU =BR%MU
		L1=          8.0_DP*AX/(3.0_DP*C11*G)
		L2=          8.0_DP*AX/(3.0_DP*C22*G)
		L3=PI_D*AX*SQRT(AX/BY)/(4.0_DP*C23*G)
		BR%BALNO(:)%OUT%UX =   AX*XIX/(MU*Z0*L1)
		BR%BALNO(:)%OUT%UY =   AX*XIY/(MU*Z0*L2)
		BR%BALNO(:)%OUT%PHX=   AB*PSI/(MU*Z0*L3)
		BR%BALNO(:)%OUT%PHY=AX*AX*PSI/(MU*Z0*L3)

		RD=BR%BALNO(1)%OUT%RD 
		selectcase(br%hcmthd) !精确计算HEATHCOTE
		case(1)
			!johnson
			BR%BALNO(:)%OUT%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)   !用来计算HERTHCOTE SLIP
		case(2)
			!refined johnson
			BR%BALNO(:)%OUT%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)*(1.0_dp-out21)   !用来计算HERTHCOTE SLIP
		case(3)
			!mzk
			rbd=br%BALNO(1)%OUT%rbd
			BR%BALNO(:)%OUT%UH=ax*(1.0_DP/(2.0_DP*RD*rbd))/(MU*Z0*L1)*(1.0_dp-out21)   !用来计算HERTHCOTE SLIP
		endselect
	  !内接触
		P  =BR%BALNO(1)%INN%P
		AX =BR%BALNO(1)%INN%AX
		BY =BR%BALNO(1)%INN%BY
		AB =BR%BALNO(1)%INN%AB
		XIX=BR%BALNO(1)%INN%XIX
		XIY=BR%BALNO(1)%INN%XIY
		PSI=BR%BALNO(1)%INN%PSI

		Z0 =P/AB/BR%NP
		C11=BR%ICOM%C11
		C22=BR%ICOM%C22
		C23=BR%ICOM%C23
		G	 =BR%G1
		MU =BR%MU
		L1=          8.0_DP*AX/(3.0_DP*C11*G)
		L2=          8.0_DP*AX/(3.0_DP*C22*G)
		L3=PI_D*AX*SQRT(AX/BY)/(4.0_DP*C23*G)
		BR%BALNO(:)%INN%UX =   AX*XIX/(MU*Z0*L1)
		BR%BALNO(:)%INN%UY =   AX*XIY/(MU*Z0*L2)
		BR%BALNO(:)%INN%PHX=   AB*PSI/(MU*Z0*L3)
		BR%BALNO(:)%INN%PHY=AX*AX*PSI/(MU*Z0*L3)

		RD=BR%BALNO(1)%INN%RD	
		selectcase(br%hcmthd) !精确计算HEATHCOTE
		case(1)
			!johnson
			BR%BALNO(:)%INN%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)    !用来计算HERTHCOTE SLIP
		case(2)
			!refined johnson
			BR%BALNO(:)%INN%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)*(1.0_dp+inn21)    !用来计算HERTHCOTE SLIP
		case(3)
			!mzk
			rbd=br%BALNO(1)%INN%rbd
			BR%BALNO(:)%INN%UH=ax*(1.0_DP/(2.0_DP*RD*rbd))/(MU*Z0*L1)*(1.0_dp+inn21)    !用来计算HERTHCOTE SLIP
		endselect
ENDSUBROUTINE GETCREEP3

SUBROUTINE GETCREEP_idx(inn21,out21,idx)
	USE FSBR
	IMPLICIT NONE
	real(dp),intent(in)::inn21,out21
	integer		:: idx
	INTEGER		:: I
	REAL(DP)	:: P,AX,BY,Z0,XIX,XIY,PSI,RD
	REAL(DP)	:: TEMP,AB
	REAL(DP)	:: L1,L2,L3,C11,C22,C23,G,MU
	real(dp)	:: rbd

	!外接触
	P  =BR%BALNO(idx)%OUT%P
	AX =BR%BALNO(idx)%OUT%AX
	BY =BR%BALNO(idx)%OUT%BY
	AB =BR%BALNO(idx)%OUT%AB
	XIX=BR%BALNO(idx)%OUT%XIX
	XIY=BR%BALNO(idx)%OUT%XIY
	PSI=BR%BALNO(idx)%OUT%PSI

	Z0=P/AB/BR%NP
	C11=BR%OCOM%C11
	C22=BR%OCOM%C22
	C23=BR%OCOM%C23
	G	 =BR%G1
	MU =BR%MU
	L1=          8.0_DP*AX/(3.0_DP*C11*G)
	L2=          8.0_DP*AX/(3.0_DP*C22*G)
	L3=PI_D*AX*SQRT(AX/BY)/(4.0_DP*C23*G)
	BR%BALNO(idx)%OUT%UX =   AX*XIX/(MU*Z0*L1)
	BR%BALNO(idx)%OUT%UY =   AX*XIY/(MU*Z0*L2)
	BR%BALNO(idx)%OUT%PHX=   AB*PSI/(MU*Z0*L3)
	BR%BALNO(idx)%OUT%PHY=AX*AX*PSI/(MU*Z0*L3)

	RD=BR%BALNO(idx)%OUT%RD 
	selectcase(br%hcmthd) !精确计算HEATHCOTE
	case(1)
		!johnson
		BR%BALNO(idx)%OUT%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)   !用来计算HERTHCOTE SLIP
	case(2)
		!refined johnson
		BR%BALNO(idx)%OUT%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)*(1.0_dp-out21)   !用来计算HERTHCOTE SLIP
	case(3)
		!mzk
		rbd=br%BALNO(idx)%OUT%rbd
		BR%BALNO(idx)%OUT%UH=ax*(1.0_DP/(2.0_DP*RD*rbd))/(MU*Z0*L1)*(1.0_dp-out21)   !用来计算HERTHCOTE SLIP
	endselect

	!内接触
	P  =BR%BALNO(idx)%INN%P
	AX =BR%BALNO(idx)%INN%AX
	BY =BR%BALNO(idx)%INN%BY
	AB =BR%BALNO(idx)%INN%AB
	XIX=BR%BALNO(idx)%INN%XIX
	XIY=BR%BALNO(idx)%INN%XIY
	PSI=BR%BALNO(idx)%INN%PSI

	Z0 =P/AB/BR%NP
	C11=BR%ICOM%C11
	C22=BR%ICOM%C22
	C23=BR%ICOM%C23
	G	 =BR%G1
	MU =BR%MU
	L1=          8.0_DP*AX/(3.0_DP*C11*G)
	L2=          8.0_DP*AX/(3.0_DP*C22*G)
	L3=PI_D*AX*SQRT(AX/BY)/(4.0_DP*C23*G)
	BR%BALNO(idx)%INN%UX =   AX*XIX/(MU*Z0*L1)
	BR%BALNO(idx)%INN%UY =   AX*XIY/(MU*Z0*L2)
	BR%BALNO(idx)%INN%PHX=   AB*PSI/(MU*Z0*L3)
	BR%BALNO(idx)%INN%PHY=AX*AX*PSI/(MU*Z0*L3)

	RD=BR%BALNO(idx)%INN%RD	
	selectcase(br%hcmthd) !精确计算HEATHCOTE
	case(1)
		!johnson
		BR%BALNO(idx)%INN%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)    !用来计算HERTHCOTE SLIP
	case(2)
		!refined johnson
		BR%BALNO(idx)%INN%UH=ax*(1.0_DP/(2.0_DP*RD**2))/(MU*Z0*L1)*(1.0_dp+inn21)    !用来计算HERTHCOTE SLIP
	case(3)
		!mzk
		rbd=br%BALNO(idx)%INN%rbd
		BR%BALNO(idx)%INN%UH=ax*(1.0_DP/(2.0_DP*RD*rbd))/(MU*Z0*L1)*(1.0_dp+inn21)    !用来计算HERTHCOTE SLIP
	endselect
ENDSUBROUTINE GETCREEP_idx