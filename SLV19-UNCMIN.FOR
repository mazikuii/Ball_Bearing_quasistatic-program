      SUBROUTINE SLV19_UNCMIN()
      implicit none

      integer nr
      parameter (nr=6)

      double precision a(nr,nr)
      external d1fcn_u,d2fcn_u
      double precision dlt
      external f19
      double precision fpls
      double precision fscale
      double precision gpls(nr)
      double precision gradtl
      integer i
      integer iagflg
      integer iahflg
      integer iexp
      integer iprint
      integer itnlim
      integer itrmcd
      integer lounit
      integer method
      integer msg
      integer n
      integer ndigit
      double precision stepmx
      double precision steptl
      double precision typsiz(nr)
      double precision wrk(nr,8)
      double precision x(nr)
      double precision xpls(nr)

	WRITE(*,*)'  ***************************************************'
	WRITE(*,*)'  ***************  BY UNCMIN METHOD  ****************'
	WRITE(*,*)'  ***************************************************'
        n = nr

	CALL X0VCT(x)
	CALL PRINTX0(x)					
c
c  Set the variables to default values.
c
        call dfault_u(n,typsiz,fscale,method,iexp,msg,ndigit,itnlim,
     &  iagflg,iahflg,iprint,lounit,dlt,gradtl,stepmx,steptl)
c
c  Change some values from their defaults.
c
        gradtl = 1.D-12
	do i=1,n
!        typsiz(i)= 1.0D+00
        typsiz(i)= abs(x(i))
	enddo
        fscale = 1.0D-12
        method = 2
        iexp = 0
		msg=0
		ndigit=-1
        itnlim = 1000
        iagflg = 0
        iahflg = 0
        iprint = 0
		dlt=-1
        steptl = 1.0D-12
        stepmx = 1.0D+02
c
c  Call OPTIF9, the full interface to UNCMIN.
c

        call optif9_u(nr,n,x,f19,d1fcn_u,d2fcn_u,typsiz,fscale,method,
     & iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,dlt,gradtl,
     &  stepmx,steptl,xpls,fpls,gpls,itrmcd,a,wrk)
c
	CALL SLV6EQS(X,gpls,.false.)
	CALL SAVEXf(X,gpls)
	CALL PRINTXF(X,gpls)

      endSUBROUTINE SLV19_UNCMIN

      subroutine f19 ( n, x, f )
 	USE INTFACE,ONLY:func
      implicit none
      integer n
      double precision f
      double precision x(n)
      f = func(X)
      endsubroutine      



