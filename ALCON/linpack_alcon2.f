C
C*    Group  Linear Solver subroutines from LINPACK (incl. BLAS)
C
      SUBROUTINE DSPFA(AP,N,KPVT,INFO)
      IMPLICIT LOGICAL(A-Z)
      INTEGER N,KPVT(1),INFO
      DOUBLE PRECISION AP(1)
C
C     sspfa factors a real symmetric matrix stored in
C     packed form by elimination with symmetric pivoting.
C
C     to solve  a*x = b , follow sspfa by sspsl.
C     to compute  inverse(a)*c , follow sspfa by sspsl.
C     to compute  determinant(a) , follow sspfa by sspdi.
C     to compute  inertia(a) , follow sspfa by sspdi.
C     to compute  inverse(a) , follow sspfa by sspdi.
C
C     on entry
C
C        ap      real (n*(n+1)/2)
C                the packed form of a symmetric matrix  a .  the
C                columns of the upper triangle are stored sequentially
C                in a one-dimensional array of length  n*(n+1)/2 .
C                see comments below for details.
C
C        n       integer
C                the order of the matrix  a .
C
C     output
C
C        ap      a block diagonal matrix and the multipliers which
C                were used to obtain it stored in packed form.
C                the factorization can be written  a = u*d*trans(u)
C                where  u  is a product of permutation and unit
C                upper triangular matrices , trans(u) is the
C                transpose of  u , and  d  is block diagonal
C                with 1 by 1 and 2 by 2 blocks.
C
C        kpvt    integer(n)
C                an integer vector of pivot indices.
C
C        info    integer
C                = 0  normal value.
C                = k  if the k-th pivot block is singular. this is
C                     not an error condition for this subroutine,
C                     but it does indicate that sspsl or sspdi may
C                     divide by zero if called.
C
C     packed storage
C
C          the following program segment will pack the upper
C          triangle of a symmetric matrix.
C
C                k = 0
C                do 20 j = 1, n
C                   do 10 i = 1, j
C                      k = k + 1
C                      ap(k)  = a(i,j)
C             10    continue
C             20 continue
C
C     linpack. this version dated 08/14/78 .
C     james bunch, univ. calif. san diego, argonne nat. lab.
C
C     subroutines and functions
C
C     blas DAXPY,DSWAP,IDAMAX
C     fortran abs,amax1,sqrt
C
C     internal variables
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IDAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
C
C     initialize
C
C     alpha is used in choosing pivot block size.
      ALPHA = (1.0D0 + SQRT(17.0D0))/8.0D0
C
      INFO = 0
C
C     main loop on k, which goes from n to 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        leave the loop if k=0 or k=1.
C
C     ...exit
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (AP(1) .EQ. 0.0E0) INFO = 1
C     ......exit
            GO TO 200
   20    CONTINUE
C
C        this section of code determines the kind of
C        elimination to be performed.  when it is completed,
C        kstep will be set to the size of the pivot block, and
C        swap will be set to .true. if an interchange is
C        required.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = ABS(AP(KK))
C
C        determine the largest off-diagonal element in
C        column k.
C
         IMAX = IDAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = ABS(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           determine the largest off-diagonal element in
C           row imax.
C
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,ABS(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = DMAX1(ROWMAX,ABS(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (ABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
C
C           column k is zero.  set info and iterate the loop.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 x 1 pivot block.
C
            IF (.NOT.SWAP) GO TO 120
C
C              perform an interchange.
C
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = AP(JK)
                  AP(JK) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           perform the elimination.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = MULK
               CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           set the pivot array.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 x 2 pivot block.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              perform an interchange.
C
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = AP(JKM1)
                  AP(JKM1) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           perform the elimination.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/AP(KM1K)
               DENOM = 1.0D0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/AP(KM1K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = MULKM1
                  CALL DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           set the pivot array.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END

Caveat receptor.  (Jack) dongarra@anl-mcs, (Eric Grosse) research!ehg
Compliments of netlib   Mon Dec  8 17:47:13 CST 1986
      SUBROUTINE DSPSL(AP,N,KPVT,B)
      IMPLICIT LOGICAL(A-Z)
      INTEGER N,KPVT(1)
      DOUBLE PRECISION AP(1),B(1)
C
C     ssisl solves the real symmetric system
C     a * x = b
C     using the factors computed by sspfa.
C
C     on entry
C
C        ap      real(n*(n+1)/2)
C                the output from sspfa.
C
C        n       integer
C                the order of the matrix  a .
C
C        kpvt    integer(n)
C                the pivot vector from sspfa.
C
C        b       real(n)
C                the right hand side vector.
C
C     on return
C
C        b       the solution vector  x .
C
C     error condition
C
C        a division by zero may occur if  sspco  has set rcond .eq. 0.0
C        or  sspfa  has set info .ne. 0  .
C
C     to compute  inverse(a) * c  where  c  is a matrix
C     with  p  columns
C           call sspfa(ap,n,kpvt,info)
C           if (info .ne. 0) go to ...
C           do 10 j = 1, p
C              call sspsl(ap,n,kpvt,c(1,j))
C        10 continue
C
C     linpack. this version dated 08/14/78 .
C     james bunch, univ. calif. san diego, argonne nat. lab.
C
C     subroutines and functions
C
C     blas DAXPY,DDOT
C     fortran iabs
C
C     internal variables.
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     loop backward applying the transformations and
C     d inverse to b.
C
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 x 1 pivot block.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 interchange.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              apply the transformation.
C
               CALL DAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           apply d inverse.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 x 2 pivot block.
C
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 interchange.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              apply the transformation.
C
               CALL DAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           apply d inverse.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/AP(KM1K)
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/AP(KM1K)
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     loop forward applying the transformations.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 x 1 pivot block.
C
            IF (K .EQ. 1) GO TO 110
C
C              apply the transformation.
C
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 interchange.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 x 2 pivot block.
C
            IF (K .EQ. 1) GO TO 140
C
C              apply the transformation.
C
               B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + DDOT(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 interchange.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
