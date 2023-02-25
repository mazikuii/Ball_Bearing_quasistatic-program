subroutine nestop(x,xold,fvec,f,g,sx,sy,ftol,xtol,mtol,btol,retcode,&
									analjac,itncount,itnlimit,maxtaken,itconsec,consecmax,termcode)
	use nrtype
	implicit none
	real(dp),dimension(:)::x,xold,fvec,g,sx,sy
	real(dp)::f
	real(dp)::ftol,xtol,mtol,btol
	integer ::retcode,itncount,itnlimit,itconsec,consecmax,termcode
	logical(lgt) analjac,maxtaken
	integer ::n

	n=size(x(:))
	termcode=0

	if(retcode==1) then
		termcode=3
	elseif(maxval(abs(sy(:)*fvec(:))) < ftol) then
		termcode=1
	elseif(maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
			1.0_dp/sx(:))) < xtol) then
		termcode=2
	elseif(itncount>=itnlimit) then
		termcode=4
	elseif(maxtaken) then
		itconsec=itconsec+1
		if(itconsec==consecmax) termcode=5
	else
		itconsec=0
		if(analjac) then
			if(maxval(abs(g(:))*max(abs(x(:)),1.0_dp/sx(:))/max(f,0.5_dp*n))<mtol) then
				termcode=6
			endif
		endif
	endif
end subroutine nestop