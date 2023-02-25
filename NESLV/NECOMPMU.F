      subroutine necompmu(r,ldr,n,epsm,mu,y)
	implicit none
      integer ldr,n
      double precision r(ldr,*),epsm,mu,y(*)

c-------------------------------------------------------------------------
c
c     Compute a small perturbation mu for the (almost) singular matrix R.
c     mu is used in the computation of the levenberg-marquardt step.
c
c     Arguments
c
c     In       R       Real(ldr,*)     upper triangular matrix from QR
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         column dimension of R
c     In       epsm    Real            machine precision
c     Out      mu      Real            sqrt(l1 norm of R * infinity norm of R
c                                      * n * epsm * 100) designed to make
c                                        trans(R)*R + mu * I not singular
c     Wk       y       Real(*)         workspace for dlange
c
c-------------------------------------------------------------------------

      double precision  aifnrm,al1nrm
      double precision dlantr

      double precision Rhund
      parameter(Rhund=100d0)

c     get the infinity norm of R
c     get the l1 norm of R

      aifnrm = dlantr('I','U','N',n,n,r,ldr,y)
      al1nrm = dlantr('1','U','N',n,n,r,ldr,y)

c     compute mu

      mu = sqrt(aifnrm*al1nrm*n*epsm*Rhund)

      return
      end