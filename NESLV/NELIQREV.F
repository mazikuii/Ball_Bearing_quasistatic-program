      subroutine liqrev(n,r,ldr,diag,b,x,sdiag,wrk)
      integer n,ldr
      double precision  r(ldr,*),diag(*),b(*),x(*),sdiag(*),wrk(*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c       In     n      Integer         order of R.
c       Inout  R      Real(ldr,*)     upper triangular matrix R from QR
c                                     unaltered
c                                     strict lower triangle contains 
c                                        transposed strict upper triangle of the upper
c                                        triangular matrix S.
c
c       In     diag   Real(*)         vector with diagonal elements of matrix D
c                                     all elements > 0
c
c       In     ldr    Integer         leading dimension of the array R.
c       In     b      Real(*)         vector of size n
c
c       Out    x      Real(*)         vector of size n
c                                     on output contains least squares solution
c                                     of the system R*x = b, D*x = 0.
c
c       Out    sdiag  Real(*)         vector of size n, containing the
c                                     diagonal elements of the upper
c                                     triangular matrix S.
c
c       Out    wrk    Real(*)         workspace of size n.
c
c     Description
c
c     Given an n by n upper triangular matrix R, a diagonal matrix D with positive entries
c     and an n-vector b, determine an x which solves the system
c
c         |R*x| = |b|
c         |D*x| = |0|
c
c     in the least squares sense.
c     The routine is used when the matrix R from the QR decomposition of a Jacobian
c     is ill-conditioned (or singular). Then it is difficult to calculate a Newton step
c     accurately (Dennis and Schnabel). D&S advise perturbing trans(J)*J with a positive
c     diagonal matrix. The idea is then to solve (J^T * J + mu*I)x=b where mu
c     is a small positive number. Using a QR decomposition of J solving this system
c     is equivalent solving (R^T*R + mu*I)x=b, where R comes from the QR decomposition.
c     Solving this system is equivalent to solving the above least squares problem.
c     On output the routine also provides an upper triangular matrix S such that
c     (see description of arguments above for the details)
c
c         (trans(R)*R + D*D) = trans(S)*S .
c
c     Method used here is described in
c     Nocedal and Wright, 2006, Numerical Optimization, Springer, ISBN 978-0-387-30303-1
c     page 259--261
c-----------------------------------------------------------------------------

      integer j,k
      double precision  bj,c,s,sum,temp
      double precision  ddot
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

c     copy R and b to preserve input and initialise S.
c     Save the diagonal elements of R in wrk.
c     Beware: the algorithm operates on an upper triangular matrix,
c     which is stored in lower triangle of R.
c     
      do 10 j=1,n
         call dcopy(n-j+1,r(j,j),ldr,r(j,j),1)
         wrk(j) = r(j,j)
   10 continue
      call dcopy(n,b,1,x,1)

c     eliminate the diagonal matrix D using givens rotations.
c     Nocedal method: start at the bottom right
c     after 100 loop had finished R contains diagonal of S
c     save in sdiag and restore original diagonal of R

      do 100 j=n,1,-1

c        initialise the row of D to be eliminated

         call nuvset(n-j+1,Rzero,sdiag(j),1)
         sdiag(j) = diag(j)

c        the transformations to eliminate the row of D

         bj = Rzero
         do 50 k=j,n

c           determine a givens rotation which eliminates the
c           appropriate element in the current row of D.
c           accumulate the transformation in the row of S.

c           eliminate the diagonal element in row j of D
c           this generates fill-in in columns [j+1 .. n] of row j of D
c           successively eliminate the fill-in with givens rotations
c           for R[j+1,j+1] and D[j,j+1].
c           rows of R have been copied into the columns of R in 10 loop
c           perform all operations on those columns to preserve the original R

            if (sdiag(k) .ne. Rzero) then

               call nuvgiv(r(k,k),sdiag(k),c,s)
               if( k .lt. n ) then
                   call drot(n-k,r(k+1,k),1,sdiag(k+1),1,c,s)
               endif

c              compute the modified element of (b,0).

               temp =  c*x(k) + s*bj
               bj   = -s*x(k) + c*bj
               x(k) = temp

            endif

   50    continue

  100 continue

c     retrieve diagonal of S from diagonal of R
c     restore original diagonal of R

      do 150 k=1,n
          sdiag(k) = r(k,k)
          r(k,k) = wrk(k)
  150 continue

c     x now contains modified b
c     solve trans(S)*x = x

      x(n) = x(n) / sdiag(n)
      do 200 j=n-1,1,-1
         sum  = ddot(n-j,r(j+1,j),1,x(j+1),1)
         x(j) = (x(j) - sum)/sdiag(j)
  200 continue

      return
      end
