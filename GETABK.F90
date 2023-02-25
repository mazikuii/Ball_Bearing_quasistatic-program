subroutine getabk(pch)
	use fsbr
	use nr, only:qromb
	use intface,only:ellip_new,pow,firstinte,secondinte
	implicit none
	type(patchabk),intent(inout),target	:: pch
	real(dp)	::newke,ke,ei1,ei2,rdiff,rx,ry,rsum,re
	common /ke/ ke

	!计算x、y方向上的有效半径（effective radius)
	rx=pch%rx1*pch%rx2/(pch%rx1+pch%rx2)
	ry=pch%ry1*pch%ry2/(pch%ry1+pch%ry2)
	rsum=rx*ry/(rx+ry)
	rdiff=abs(ry-rx)/(rx+ry)!计算球综合主曲率差
	re=2.0_dp/((1.0_dp-pch%nu1**2)/pch%e1+(1.0_dp-pch%nu2**2)/pch%e2)
	if(rx==ry) then
		ei1=pio2_d
		ei2=pio2_d
		ke=1.0_dp
		pch%xlarger=.true.
	else
		!叠代求解接触椭圆率
		ke=1.0_dp/sqrt(1.0_dp-(min(rx,ry)/max(rx,ry))**(4.0_dp/3.0_dp))	!初值
		newke=ellip_new(rdiff)
		do while(abs(newke-ke)>pch%kp)
			ke=newke
			newke=ellip_new(rdiff)
		end do
		ke=newke
		ei1=qromb(firstinte ,0.0_dp,pio2_d)
		ei2=qromb(secondinte,0.0_dp,pio2_d)
		!为了计算微积分方程系数的需要，必须使椭圆率在0～1之间
		if(rx<ry) then
			ke=1.0_dp/ke
			pch%xlarger=.false.
			!重新计算第一类和第二类椭圆积分值
			ei1=qromb(firstinte ,0.0_dp,pio2_d)
			ei2=qromb(secondinte,0.0_dp,pio2_d)
		else
			pch%xlarger=.true.
		end if	
	endif
	pch%ke=ke
	!计算干接触工况下接触椭圆半轴长度
	pch%and=pow(6.0_dp*(ke**2)*ei2*rsum/pi_d/re,1.0_dp/3.0_dp)
	pch%bnd=pow(6.0_dp*ei2*rsum/pi_d/ke/re,1.0_dp/3.0_dp)
	!计算干接触工况下趋近量
	pch%appnd=ei1*pow((9.0_dp/2.0_dp/ei2/rsum)*pow(1.0_dp/pi_d/ke/re,2.0_dp),1.0_dp/3.0_dp) 
	pch%rx=rx
	pch%ry=ry
	pch%ei1=ei1
	pch%ei2=ei2
endsubroutine getabk