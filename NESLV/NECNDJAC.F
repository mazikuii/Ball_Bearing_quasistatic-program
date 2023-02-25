      subroutine NEcndjac(n,r,ldr,epsm,rcond,y,rcdwrk,icdwrk,ierr,mu)
	implicit none
      integer n,ldr,icdwrk(*),ierr
      double precision epsm,rcond,mu,r(ldr,*),y(*),rcdwrk(*)

c---------------------------------------------------------------------
c
c     Check r for singularity and/or ill conditioning
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       r       Real(ldr,*)     upper triangular R from QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       epsm    Real            machine precision
c     Out      rcond   Real            inverse condition  of r
c     Wk       y       Real(*)         workspace
c     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
c     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c     Out      mu      Real            0 if ierr == 0
c                                      small positive number when ierr > 0
c                                      to make trans(R)*R+mu*I non singular
c
c---------------------------------------------------------------------

      integer i,info
      logical rsing
      double precision Rzero,R2d3
      parameter(Rzero=0.0d0, R2d3=2.0d0/3.0d0)

      mu = Rzero
      ierr = 0

      rsing = .false.
      do 10 i=1,n
         if( r(i,i) .eq. Rzero ) then
             rsing = .true.
         endif
   10 continue

      if( rsing ) then
         ierr = 2
         rcond = Rzero
      else
         call dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
         if( rcond .lt. epsm**R2d3 ) then
             ierr = 1
         endif
      endif

      if( ierr .gt. 0 ) then
         call necompmu(r,ldr,n,epsm,mu,y)
      endif

      return
      end
