
module neintface
	interface
		subroutine nedflt(analjac,restart,maxtaken,&
											lshmsd,global,method,termcode,&
											itncount,itnlimit,itconsec,consecmax,nmt,&
											ftol,xtol,mtol,btol,etak,stepmx,dlt,&
											x0,sx,sy)
			use nrtype
			implicit none
			logical::analjac,restart,maxtaken
			integer::lshmsd,global,method,termcode
			integer::itncount,itnlimit,itconsec,consecmax,nmt
			real(dp)::ftol,xtol,mtol,btol,etak,stepmx,dlt
			real(dp)::x0(:),sx(:),sy(:)
		endsubroutine nedflt
	end interface

	interface
		function fv24(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	fv24(size(x(:)))
		endfunction fv24	
		function fv25(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	fv25(size(x(:)))
		endfunction fv25	
	end interface

	interface
		subroutine jac24(x,jc)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp),intent(out)::	jc(:,:)
		end subroutine jac24	
	end interface

	interface
		function nefmin(x,func,fvec,sy)
			use nrtype
			implicit none
			real(dp),dimension(:), intent(in) :: x,sy
			real(dp),dimension(:)::fvec
			real(dp)::nefmin
			interface
				function func(x)
					use nrtype
					implicit none
					real(dp),intent(in)	::	x(:)
					real(dp)::	func(size(x(:)))
				endfunction func	
			end interface
		endfunction nefmin
	end interface

	interface
		subroutine nestop0(fvec,typy,ftol,itconsec,termcode)
			use nrtype
			real(dp),dimension(:)::fvec,typy
			real(dp)::ftol
			integer ::itconsec,termcode
		endsubroutine nestop0
	end interface

	interface
		subroutine nefdjac(x,func,fvec,df,sx)
			use nrtype; use nrutil, only : assert_eq
			implicit none
			real(dp), dimension(:), intent(inout) :: x
			real(dp), dimension(size(x)), intent(in) :: fvec,sx
			real(dp), dimension(size(x),size(x)), intent(out) :: df
			interface
				function func(x)
					use nrtype
					implicit none
					real(dp),intent(in)	::	x(:)
					real(dp)::	func(size(x(:)))
				endfunction func	
			end interface
		end subroutine nefdjac	
		subroutine nefdjac4(x,func,fvec,df,sx)
			use nrtype; use nrutil, only : assert_eq
			implicit none
			real(dp), dimension(:), intent(inout) :: x
			real(dp), dimension(size(x)), intent(in) :: fvec,sx
			real(dp), dimension(size(x),size(x)), intent(out) :: df
			interface
				function func(x)
					use nrtype
					implicit none
					real(dp),intent(in)	::	x(:)
					real(dp)::	func(size(x(:)))
				endfunction func	
			end interface
		end subroutine nefdjac4	
	end interface

	interface
		subroutine neslv(x,fp,gp,func,fjac,typx,typy,&
										 analjac,restart,maxtaken,&
										 lshmsd,global,method,termcode,&
										 itncount,itnlimit,itconsec,consecmax,nmt,&
										 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
										 nfnt,njnt)
			use nrtype
			use nrutil,only:vabs,put_diag,unit_matrix
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
		endsubroutine
		subroutine neslv4(x,fp,gp,func,fjac,typx,typy,&
										 analjac,restart,maxtaken,&
										 lshmsd,global,method,termcode,&
										 itncount,itnlimit,itconsec,consecmax,nmt,&
										 ftol,xtol,mtol,btol,eta0,stepmx,dlt,&
										 nfnt,njnt)
			use nrtype
			use nrutil,only:vabs,put_diag,unit_matrix
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
		endsubroutine
	end interface

	interface
		subroutine nefmqr(qt,r)
			use nrtype
			implicit none
			real(dp), dimension(:,:) :: qt,r
		end subroutine nefmqr
	end interface

	interface
		subroutine neslvsn(qt,r,fvec,p,g,qtf,sx,sy,macheps)
			use nrtype
			implicit none
			real(dp), dimension(:,:) :: qt,r
			real(dp), dimension(:)	 :: fvec,p,g,qtf,sx,sy
			real(dp)								 :: macheps
		end subroutine neslvsn
	end interface

	interface
		subroutine nelnsh(xold,fold,g,p,x,f,fvec,stpmax,btol,sx,sy,&
											maxtaken,lshmsd,etak,nmt,retcode,fmin,func,nfnt)
			use nrtype
			implicit none
			real(dp):: xold(:),fold,g(:),p(:),x(:),f,fvec(:),stpmax,btol,sx(:),sy(:),etak
			logical	:: maxtaken
			integer :: lshmsd,nmt,retcode,nfnt
			interface
				function fmin(x,func,fvec,sy)
					use nrtype
					implicit none
					real(dp),dimension(:), intent(in) :: x,sy
					real(dp),dimension(:)::fvec
					real(dp)::fmin
					interface
						function func(x)
							use nrtype
							implicit none
							real(dp),intent(in)	::	x(:)
							real(dp)::	func(size(x(:)))
						endfunction func	
					end interface
				endfunction fmin
				function func(x)
					use nrtype
					implicit none
					real(dp),intent(in)	::	x(:)
					real(dp)::	func(size(x(:)))
				endfunction func	
			end interface
		end subroutine nelnsh
	end interface

	interface
		subroutine nestop(x,xold,fvec,f,g,sx,sy,ftol,xtol,mtol,btol,retcode,&
									analjac,itncount,itnlimit,maxtaken,itconsec,consecmax,termcode)
			use nrtype
			implicit none
			real(dp),dimension(:)::x,xold,fvec,g,sx,sy
			real(dp)::f
			real(dp)::ftol,xtol,mtol,btol
			logical(lgt) analjac,maxtaken
			integer ::n,retcode,itncount,itnlimit,itconsec,consecmax,termcode
		end subroutine nestop
	end interface

	interface
		subroutine neinchk(ftol,xtol,mtol,btol,etak,itnlimit,sx,sy,&
											 restart,maxtaken,&
											 lshmsd,global,method,termcode,&
											 itncount,itconsec,consecmax,nmt,&
											 macheps)
			use nrtype
			implicit none
			real(dp)::ftol,xtol,mtol,btol,etak
			integer::itnlimit
			real(dp)::sx(:),sy(:)
			logical::restart,maxtaken
			integer::lshmsd,global,method,termcode
			integer::itncount,itconsec,consecmax,nmt
			real(dp)::macheps
		endsubroutine neinchk
	end interface

	interface
		subroutine neqrudt(qt,r,x,xold,f,fold,fvec,fvcold,sx,sy,stpmax,macheps)
			use nrtype; use nrutil, only : nrerror,get_diag
			use nr, only : qrupdt
			implicit none
			real(dp), dimension(:,:) :: qt,r
			real(dp), dimension(:) :: x,xold
			real(dp)::f,fold
			real(dp), dimension(:) :: fvec,fvcold,sx,sy
			real(dp)::stpmax,macheps
		endsubroutine neqrudt	
	end interface

	interface
		function fvtt(x)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	fvtt(size(x(:)))
		endfunction fvtt	
		subroutine jactt(x,jc)
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp),intent(out)::	jc(:,:)
		end subroutine jactt	
	end interface

	interface
		subroutine neexpln ( itrmcd,consecmax)
			implicit none
			integer itrmcd,consecmax
		endsubroutine neexpln
	end interface

	interface
		subroutine svfm(f,fm)
			use nrtype
			implicit none
			real(dp)::f,fm(:)
		end subroutine
	end interface

	interface
		subroutine conangs(alfi,alfo)
			use	nrtype
			implicit none
			real(dp),intent(inout):: alfi,alfo
		endsubroutine conangs
	end interface

	interface
		subroutine slv3eqs_st(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	f(:)
			logical	:: sv
		endsubroutine slv3eqs_st
		subroutine slv3eqs_creep(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	f(:)
			logical	:: sv
		endsubroutine slv3eqs_creep
		recursive subroutine slv3eqs_inert(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	f(:)
			logical	:: sv
		endsubroutine slv3eqs_inert
	end interface

	interface
		subroutine slv6eqs(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(:), intent(out) :: f
			logical				::sv
		endsubroutine slv6eqs
		subroutine slv8eqs(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(:), intent(out) :: f
			logical				::sv
		endsubroutine slv8eqs
		subroutine slv6eqs_idx(x,f,sv,index,mthd)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(:), intent(out) :: f
			logical				::sv
			integer				::index,mthd
		endsubroutine slv6eqs_idx
		subroutine slv6eqs_creep(x)
			use	nrtype
			USE FSBR
			implicit none
			real(dp)::x(:)
		endsubroutine slv6eqs_creep
	end interface

	interface
		subroutine slv2eqs_sym(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:):: x
			real(dp), dimension(:):: f
			logical	:: sv
		endsubroutine slv2eqs_sym
		subroutine slv2eqs_asym(x,f,sv)
			use fsbr
			use nrtype
			implicit none
			real(dp), dimension(:):: x
			real(dp), dimension(:):: f
			logical	:: sv
		endsubroutine slv2eqs_asym
	end interface

	interface
		function qs2eqs_sym(x)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	qs2eqs_sym(size(x(:)))
		endfunction qs2eqs_sym
		function qs2eqs_asym(x)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	qs2eqs_asym(size(x(:)))
		endfunction qs2eqs_asym
	end interface

	interface
		function qs3eqs_st(x)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	qs3eqs_st(size(x(:)))
		endfunction qs3eqs_st
		function qs3eqs_creep(x)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	qs3eqs_creep(size(x(:)))
		endfunction qs3eqs_creep
		function qs3eqs_inert(x)
			use fsbr
			use nrtype
			implicit none
			real(dp),intent(in)	::	x(:)
			real(dp)::	qs3eqs_inert(size(x(:)))
		endfunction qs3eqs_inert
	end interface

	interface
		subroutine x0vct(x)
			use nrtype
			implicit none
			real(dp), dimension(:) :: x
		end subroutine x0vct
	endinterface


endmodule neintface