recursive subroutine neslv4(x,fp,gp,func,fjac,typx,typy,&
								 analjac,restart,maxtaken,&
								 lshmsd,global,method,termcode,&
								 itncount,itnlimit,itconsec,consecmax,nmt,&
								 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
								 nfnt,njnt)
	use nrtype
	use nrutil,only:vabs,put_diag,unit_matrix
	use neintface
	implicit none
	real(dp),dimension(:)::x,fp,gp
	interface
		function func(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	func(size(x(:)))
		endfunction func
		subroutine fjac(x,jc)
			use nrtype
			implicit none
			real(dp),intent(in)::x(:)
			real(dp)::	jc(:,:)
		end subroutine fjac	
	end interface
	real(dp),dimension(:)::typx,typy
	logical::analjac,restart,maxtaken
	integer::lshmsd,global,method,termcode
	integer::itncount,itnlimit,itconsec,consecmax,nmt
	real(dp)::ftol,xtol,mtol,btol,eta0,stepmx,dlt
	integer::nfnt,njnt

	real(dp),dimension(size(x),size(x))::jc,qt,r
	real(dp),dimension(size(x))::xold,fvec,fvcold,g,p,sx,sy,qtf
	real(dp) ::f,fold,stpmax,dlt0,macheps
	integer ::nstp
	integer	:: n,retcode
	!for ddlg
	real(dp)::d(size(x)),ssd(size(x)),v(size(x)),wa(size(x))
	integer	::priter=0
	real(dp)::sigma=.4_dp
	real(dp)::ared,pred,rk
	real(dp):: DNRM2

	n=size(x(:))
	where(typx==0.0_dp) typx=1.0_dp
	where(typy==0.0_dp) typy=1.0_dp
	sx(:)=1.0_dp/abs(typx(:))
	sy(:)=1.0_dp/abs(typy(:))
	nfnt=0
	njnt=0
	call neinchk(ftol,xtol,mtol,btol,eta0,itnlimit,sx,sy,&
							 restart,maxtaken,&
							 lshmsd,global,method,termcode,&
							 itncount,itconsec,consecmax,nmt,&
							 macheps)
	if(termcode<0) return	!error, subroutine return
	f=nefmin(x,func,fvec,sy)	!output f and fvec
	nfnt=nfnt+1
	call nestop0(fvec,sy,ftol,itconsec,termcode)
	if(termcode==1) return	!solved,无需进一步求解

	stpmax=stepmx*max(vabs(x(:)*sx),real(n,dp))
	dlt0=dlt
	do while(termcode==0)
		itncount=itncount+1
		if(analjac) then
			call fjac(x,jc) !analytical jacobian
			njnt=njnt+1
			fvec(:)=func(x)
			r(:,:)=jc(:,:)
			call nefmqr(qt,r)
			call unit_matrix(jc(:,:))
			call put_diag(sy(:),jc(:,:))
			g(:)=matmul(matmul(qt(:,:),jc(:,:)),fvec(:))
		else
			if(restart .or. method==1) then
				fvec(:)=func(x)
				nfnt=nfnt+1
				call nefdjac4(x,func,fvec,jc,sx) !diff jac
				nfnt=nfnt+n
				njnt=njnt+1
				call unit_matrix(qt(:,:))
				call put_diag(sy(:),qt(:,:))
				r=matmul(qt(:,:),jc(:,:))
				call nefmqr(qt,r)
			else
				call neqrudt(qt,r,x,xold,f,fold,fvec,fvcold,sx,sy,stpmax,macheps)	
			endif
		endif
		call dcopy(n,x,1,xold,1)
		call dcopy(n,fvec,1,fvcold,1)
		fold=f
		call neslvsn(qt,r,fvec,p,g,qtf,sx,sy,macheps)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!global strategy
		select case(global)
		case(1) 
			!retcode=0	0:search made;1:search failed;
			!monotone/nonmonotone line search
			call nelnsh(xold,fold,g,p,x,f,fvec,stpmax,btol,sx,sy,&	!fold changed
									maxtaken,lshmsd,eta0,nmt,retcode,nefmin,func,nfnt)
		case(2)
			!double dogleg method
			stpmax=.1_dp
			call neddlg(n,r,n,p,g,xold,fold,stpmax,&
									btol,maxtaken,dlt,qtf,sx,sy,&
									func,d,x,ssd,v,wa,fvec,&
									x,fvec,f,retcode,nstp,priter,itncount)
			nfnt=nfnt+nstp
		case(3)
			!simplified dogleg method
			stpmax=100._dp
			call nedogleg(n,r,sx,qtf,stpmax,p)
			x=x+p
			f=nefmin(x,func,fvec,sy)
		case(4)
			!dogleg method
			dlt=abs(dlt)
			call nedogleg(n,r,sx,qtf,dlt,p)
			x=x+p
			f=nefmin(x,func,fvec,sy)
			ared=DNRM2(n,fvcold,1)**2-DNRM2(n,fvec,1)**2
			call nefdjac4(xold,func,fvcold,jc,sx) !计算JACOBIAN
			ssd=fvcold+matmul(jc,p)
			pred=DNRM2(n,fvcold,1)**2-DNRM2(n,ssd,1)**2
			rk=ared/pred
			if(rk>=0.05_dp) then	
				dlt=dlt*2.0_dp		!需要对此处进行调整
			elseif(rk>=0.001_dp) then
				dlt=dlt
			else
!				dlt=.5_dp*dlt	!需要对此处进行调整
			endif
		endselect
		call nestop(x,xold,fvec,f,g,sx,sy,ftol,xtol,mtol,btol,retcode,&
									analjac,itncount,itnlimit,maxtaken,itconsec,consecmax,termcode)
		if(termcode==3 .and. not(restart) .and. not(analjac)) then
			termcode=0
			dlt=dlt0 !disregard it if global=1
			restart=.true.		!
			cycle
		elseif(termcode>0) then		!程序的唯一正常出口
			call dcopy(n,fvec,1,fp,1)
			call dcopy(n,x,1,x,1)
			call dcopy(n,g,1,gp,1)	!前一点的梯度，如需输出精确梯度需重新计算
		else
			restart=.false.
		endif
	enddo
endsubroutine neslv4