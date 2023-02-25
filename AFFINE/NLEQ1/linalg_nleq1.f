      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN 
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH      
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum 
* 
      SFMIN = DLAMCH('S')  
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN 
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN 
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 ) 
               ELSE 
                 DO 20 I = 1, M-J 
                    A( J+I, J ) = A( J+I, J ) / A( J, J ) 
   20            CONTINUE 
               END IF 
            END IF 
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U, because of fill-in resulting from the row
*  interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in.
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Gaussian elimination with partial pivoting
*
*     Set fill-in elements in columns KU+2 to KV to zero.
*
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     JU is the index of the last column affected by the current stage
*     of the factorization.
*
      JU = 1
*
      DO 40 J = 1, MIN( M, N )
*
*        Set fill-in elements in column J+KV to zero.
*
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
*
*        Find pivot and test for singularity. KM is the number of
*        subdiagonal elements in the current column.
*
         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
*
*           Apply interchange to columns J to JU.
*
            IF( JP.NE.1 )
     $         CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     $                     AB( KV+1, J ), LDAB-1 )
*
            IF( KM.GT.0 ) THEN
*
*              Compute multipliers.
*
               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
*
*              Update trailing submatrix within the band.
*
               IF( JU.GT.J )
     $            CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     $                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     $                       LDAB-1 )
            END IF
         ELSE
*
*           If pivot is zero, set INFO to the index of the pivot
*           unless a zero pivot has already been found.
*
            IF( INFO.EQ.0 )
     $         INFO = J
         END IF
   40 CONTINUE
      RETURN
*
*     End of DGBTF2
*
      END
      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTRF computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U because of fill-in resulting from the row interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, K2, KM, KV, NB, NW
      DOUBLE PRECISION   TEMP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX, ILAENV
      EXTERNAL           IDAMAX, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL,
     $                   DSWAP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
*
*     The block size must not exceed the limit set by the size of the
*     local arrays WORK13 and WORK31.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
*
*        Use unblocked code
*
         CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
*
*        Use blocked code
*
*        Zero the superdiagonal elements of the work array WORK13
*
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
*
*        Zero the subdiagonal elements of the work array WORK31
*
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
*
*        Gaussian elimination with partial pivoting
*
*        Set fill-in elements in columns KU+2 to KV to zero
*
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
*
*        JU is the index of the last column affected by the current
*        stage of the factorization
*
         JU = 1
*
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
*
*           The active part of the matrix is partitioned
*
*              A11   A12   A13
*              A21   A22   A23
*              A31   A32   A33
*
*           Here A11, A21 and A31 denote the current block of JB columns
*           which is about to be factorized. The number of rows in the
*           partitioning are JB, I2, I3 respectively, and the numbers
*           of columns are JB, J2, J3. The superdiagonal elements of A13
*           and the subdiagonal elements of A31 lie outside the band.
*
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
*
*           J2 and J3 are computed after JU has been updated.
*
*           Factorize the current block of JB columns
*
            DO 80 JJ = J, J + JB - 1
*
*              Set fill-in elements in column JJ+KV to zero
*
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
*
*              Find pivot and test for singularity. KM is the number of
*              subdiagonal elements in the current column.
*
               KM = MIN( KL, M-JJ )
               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
*
*                    Apply interchange to columns J to J+JB-1
*
                     IF( JP+JJ-1.LT.J+KL ) THEN
*
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
*
*                       The interchange affects columns J to JJ-1 of A31
*                       which are stored in the work array WORK31
*
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     $                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
*
*                 Compute multipliers
*
                  CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                        1 )
*
*                 Update trailing submatrix within the band and within
*                 the current block. JM is the index of the last column
*                 which needs to be updated.
*
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     $               CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
*
*                 If pivot is zero, set INFO to the index of the pivot
*                 unless a zero pivot has already been found.
*
                  IF( INFO.EQ.0 )
     $               INFO = JJ
               END IF
*
*              Copy current column of A31 into the work array WORK31
*
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
*
*              Apply the row interchanges to the other blocks.
*
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
*
*              Use DLASWP to apply the row interchanges to A12, A22, and
*              A32.
*
               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     $                      IPIV( J ), 1 )
*
*              Adjust the pivot indices.
*
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
*
*              Apply the row interchanges to A13, A23, and A33
*              columnwise.
*
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
*
*              Update the relevant part of the trailing submatrix
*
               IF( J2.GT.0 ) THEN
*
*                 Update A12
*
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A22
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A32
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
*
               IF( J3.GT.0 ) THEN
*
*                 Copy the lower triangle of A13 into the work array
*                 WORK13
*
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
*
*                 Update A13 in the work array
*
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A23
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A33
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
*
*                 Copy the lower triangle of A13 back into place
*
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
*
*              Adjust the pivot indices.
*
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
*
*           Partially undo the interchanges in the current block to
*           restore the upper triangular form of A31 and copy the upper
*           triangle of A31 back into place
*
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
*
*                 Apply interchange to columns J to JJ-1
*
                  IF( JP+JJ-1.LT.J+KL ) THEN
*
*                    The interchange does not affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
*
*                    The interchange does affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
*
*              Copy the current column of A31 back into place
*
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
*
      RETURN
*
*     End of DGBTRF
*
      END
      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general band matrix A using the LU factorization computed
*  by DGBTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     $                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
*
         DO 20 I = 1, NRHS
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     $                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
*
      ELSE
*
*        Solve A'*X = B.
*
         DO 30 I = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     $                  LDAB, B( 1, I ), 1 )
   30    CONTINUE
*
*        Solve L'*X = B, overwriting B with X.
*
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     $                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGBTRS
*
      END
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*    .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*    ..
*
* Purpose
* =======
*
* DGER performs the rank 1 operation
*
*    A := alpha*x*y' + A,
*
* where alpha is a scalar, x is an m element vector, y is an n element
* vector and A is an m by n matrix.
*
* Arguments
* ==========
*
* M       - INTEGER.
*             On entry, M specifies the number of rows of the matrix A.
*             M must be at least zero.
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the number of columns of the matrix A.
*             N must be at least zero.
*             Unchanged on exit.
*
* ALPHA - DOUBLE PRECISION.
*             On entry, ALPHA specifies the scalar alpha.
*             Unchanged on exit.
*
* X       - DOUBLE PRECISION array of dimension at least
*             ( 1 + ( m - 1 )*abs( INCX ) ).
*             Before entry, the incremented array X must contain the m
*             element vector x.
*             Unchanged on exit.
*
* INCX - INTEGER.
*             On entry, INCX specifies the increment for the elements of
*             X. INCX must not be zero.
*             Unchanged on exit.
*
* Y       - DOUBLE PRECISION array of dimension at least
*             ( 1 + ( n - 1 )*abs( INCY ) ).
*             Before entry, the incremented array Y must contain the n
*             element vector y.
*             Unchanged on exit.
*
* INCY - INTEGER.
*             On entry, INCY specifies the increment for the elements of
*             Y. INCY must not be zero.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*             Before entry, the leading m by n part of the array A must
*             contain the matrix of coefficients. On exit, A is
*             overwritten by the updated matrix.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. LDA must be at least
*             max( 1, m ).
*             Unchanged on exit.
*
*
* Level 2 Blas routine.
*
* -- Written on 22-October-1986.
*    Jack Dongarra, Argonne National Lab.
*    Jeremy Du Croz, Nag Central Office.
*    Sven Hammarling, Nag Central Office.
*    Richard Hanson, Sandia National Labs.
*
*
*    .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX
*    ..
*
*    Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
         INFO = 1
      ELSE IF (N.LT.0) THEN
         INFO = 2
      ELSE IF (INCX.EQ.0) THEN
         INFO = 5
      ELSE IF (INCY.EQ.0) THEN
         INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = 9
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGER ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*    Start the operations. In this version the elements of A are
*    accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
         JY = 1
      ELSE
         JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
         DO 20 J = 1,N
            IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                        A(I,J) = A(I,J) + X(I)*TEMP
10             CONTINUE
            END IF
            JY = JY + INCY
20    CONTINUE
      ELSE
         IF (INCX.GT.0) THEN
            KX = 1
         ELSE
            KX = 1 - (M-1)*INCX
         END IF
         DO 40 J = 1,N
            IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
30             CONTINUE
            END IF
            JY = JY + INCY
40    CONTINUE
      END IF
*
      RETURN
*
*    End of DGER .
*
      END


      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
*    .. Scalar Arguments ..
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*    ..
*
* Purpose
* =======
*
* DTBSV solves one of the systems of equations
*
*    A*x = b, or A'*x = b,
*
* where b and x are n element vectors and A is an n by n unit, or
* non-unit, upper or lower triangular band matrix, with ( k + 1 )
* diagonals.
*
* No test for singularity or near-singularity is included in this
* routine. Such tests must be performed before calling this routine.
*
* Arguments
* ==========
*
* UPLO - CHARACTER*1.
*             On entry, UPLO specifies whether the matrix is an upper or
*             lower triangular matrix as follows:
*
*             UPLO = 'U' or 'u' A is an upper triangular matrix.
*
*             UPLO = 'L' or 'l' A is a lower triangular matrix.
*
*             Unchanged on exit.
*
* TRANS - CHARACTER*1.
*             On entry, TRANS specifies the equations to be solved as
*             follows:
*
*             TRANS = 'N' or 'n' A*x = b.
*
*             TRANS = 'T' or 't' A'*x = b.
*
*             TRANS = 'C' or 'c' A'*x = b.
*
*             Unchanged on exit.
*
* DIAG - CHARACTER*1.
*             On entry, DIAG specifies whether or not A is unit
*             triangular as follows:
*
*             DIAG = 'U' or 'u' A is assumed to be unit triangular.
*
*             DIAG = 'N' or 'n' A is not assumed to be unit
*                                     triangular.
*
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the order of the matrix A.
*             N must be at least zero.
*             Unchanged on exit.
*
* K       - INTEGER.
*             On entry with UPLO = 'U' or 'u', K specifies the number of
*             super-diagonals of the matrix A.
*             On entry with UPLO = 'L' or 'l', K specifies the number of
*             sub-diagonals of the matrix A.
*             K must satisfy 0 .le. K.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*             Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*             by n part of the array A must contain the upper triangular
*             band part of the matrix of coefficients, supplied column by
*             column, with the leading diagonal of the matrix in row
*             ( k + 1 ) of the array, the first super-diagonal starting at
*             position 2 in row k, and so on. The top left k by k triangle
*             of the array A is not referenced.
*             The following program segment will transfer an upper
*             triangular band matrix from conventional full matrix storage
*             to band storage:
*
*                   DO 20, J = 1, N
*                      M = K + 1 - J
*                      DO 10, I = MAX( 1, J - K ), J
*                         A( M + I, J ) = matrix( I, J )
*             10 CONTINUE
*             20 CONTINUE
*
*             Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*             by n part of the array A must contain the lower triangular
*             band part of the matrix of coefficients, supplied column by
*             column, with the leading diagonal of the matrix in row 1 of
*             the array, the first sub-diagonal starting at position 1 in
*             row 2, and so on. The bottom right k by k triangle of the
*             array A is not referenced.
*             The following program segment will transfer a lower
*             triangular band matrix from conventional full matrix storage
*             to band storage:
*
*                   DO 20, J = 1, N
*                      M = 1 - J
*                      DO 10, I = J, MIN( N, J + K )
*                         A( M + I, J ) = matrix( I, J )
*             10 CONTINUE
*             20 CONTINUE
*
*             Note that when DIAG = 'U' or 'u' the elements of the array A
*             corresponding to the diagonal elements of the matrix are not
*             referenced, but are assumed to be unity.
*             Unchanged on exit.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. LDA must be at least
*             ( k + 1 ).
*             Unchanged on exit.
*
* X       - DOUBLE PRECISION array of dimension at least
*             ( 1 + ( n - 1 )*abs( INCX ) ).
*             Before entry, the incremented array X must contain the n
*             element right-hand side vector b. On exit, X is overwritten
*             with the solution vector x.
*
* INCX - INTEGER.
*             On entry, INCX specifies the increment for the elements of
*             X. INCX must not be zero.
*             Unchanged on exit.
*
*
* Level 2 Blas routine.
*
* -- Written on 22-October-1986.
*    Jack Dongarra, Argonne National Lab.
*    Jeremy Du Croz, Nag Central Office.
*    Sven Hammarling, Nag Central Office.
*    Richard Hanson, Sandia National Labs.
*
*
*    .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOUNIT
*    ..
*    .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
*    ..
*
*    Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
         INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +       .NOT.LSAME(TRANS,'C')) THEN
         INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
         INFO = 7
      ELSE IF (INCX.EQ.0) THEN
         INFO = 9
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DTBSV ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOUNIT = LSAME(DIAG,'N')
*
*    Set up the start point in X if the increment is not unity. This
*    will be ( N - 1 )*INCX too small for descending loops.
*
      IF (INCX.LE.0) THEN
         KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
         KX = 1
      END IF
*
*    Start the operations. In this version the elements of A are
*    accessed by sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*       Form x := inv( A )*x.
*
         IF (LSAME(UPLO,'U')) THEN
            KPLUS1 = K + 1
            IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                        IF (X(J).NE.ZERO) THEN
                              L = KPLUS1 - J
                              IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                              TEMP = X(J)
                              DO 10 I = J - 1,MAX(1,J-K),-1
                                 X(I) = X(I) - TEMP*A(L+I,J)
10                         CONTINUE
                        END IF
20             CONTINUE
            ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                        KX = KX - INCX
                        IF (X(JX).NE.ZERO) THEN
                              IX = KX
                              L = KPLUS1 - J
                              IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                              TEMP = X(JX)
                              DO 30 I = J - 1,MAX(1,J-K),-1
                                 X(IX) = X(IX) - TEMP*A(L+I,J)
                                 IX = IX - INCX
30                         CONTINUE
                        END IF
                        JX = JX - INCX
40             CONTINUE
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                        IF (X(J).NE.ZERO) THEN
                              L = 1 - J
                              IF (NOUNIT) X(J) = X(J)/A(1,J)
                              TEMP = X(J)
                              DO 50 I = J + 1,MIN(N,J+K)
                                 X(I) = X(I) - TEMP*A(L+I,J)
50                         CONTINUE
                        END IF
60             CONTINUE
            ELSE
                  JX = KX
                  DO 80 J = 1,N
                        KX = KX + INCX
                        IF (X(JX).NE.ZERO) THEN
                              IX = KX
                              L = 1 - J
                              IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                              TEMP = X(JX)
                              DO 70 I = J + 1,MIN(N,J+K)
                                 X(IX) = X(IX) - TEMP*A(L+I,J)
                                 IX = IX + INCX
70                         CONTINUE
                        END IF
                        JX = JX + INCX
80             CONTINUE
            END IF
         END IF
      ELSE
*
*       Form x := inv( A')*x.
*
         IF (LSAME(UPLO,'U')) THEN
            KPLUS1 = K + 1
            IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                        TEMP = X(J)
                        L = KPLUS1 - J
                        DO 90 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(I)
90                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                        X(J) = TEMP
100             CONTINUE
            ELSE
                  JX = KX
                  DO 120 J = 1,N
                        TEMP = X(JX)
                        IX = KX
                        L = KPLUS1 - J
                        DO 110 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX + INCX
110                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                        X(JX) = TEMP
                        JX = JX + INCX
                        IF (J.GT.K) KX = KX + INCX
120             CONTINUE
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                        TEMP = X(J)
                        L = 1 - J
                        DO 130 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(I)
130                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(1,J)
                        X(J) = TEMP
140             CONTINUE
            ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                        TEMP = X(JX)
                        IX = KX
                        L = 1 - J
                        DO 150 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX - INCX
150                   CONTINUE
                        IF (NOUNIT) TEMP = TEMP/A(1,J)
                        X(JX) = TEMP
                        JX = JX - INCX
                        IF ((N-J).GE.K) KX = KX - INCX
160             CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*    End of DTBSV .
*
      END

      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*    .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*    ..
*    .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*    ..
*
* Purpose
* =======
*
* DGEMM performs one of the matrix-matrix operations
*
*    C := alpha*op( A )*op( B ) + beta*C,
*
* where op( X ) is one of
*
*    op( X ) = X or op( X ) = X',
*
* alpha and beta are scalars, and A, B and C are matrices, with op( A )
* an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
*
* Arguments
* ==========
*
* TRANSA - CHARACTER*1.
*             On entry, TRANSA specifies the form of op( A ) to be used in
*             the matrix multiplication as follows:
*
*             TRANSA = 'N' or 'n', op( A ) = A.
*
*             TRANSA = 'T' or 't', op( A ) = A'.
*
*             TRANSA = 'C' or 'c', op( A ) = A'.
*
*             Unchanged on exit.
*
* TRANSB - CHARACTER*1.
*             On entry, TRANSB specifies the form of op( B ) to be used in
*             the matrix multiplication as follows:
*
*             TRANSB = 'N' or 'n', op( B ) = B.
*
*             TRANSB = 'T' or 't', op( B ) = B'.
*
*             TRANSB = 'C' or 'c', op( B ) = B'.
*
*             Unchanged on exit.
*
* M       - INTEGER.
*             On entry, M specifies the number of rows of the matrix
*             op( A ) and of the matrix C. M must be at least zero.
*             Unchanged on exit.
*
* N       - INTEGER.
*             On entry, N specifies the number of columns of the matrix
*             op( B ) and the number of columns of the matrix C. N must be
*             at least zero.
*             Unchanged on exit.
*
* K       - INTEGER.
*             On entry, K specifies the number of columns of the matrix
*             op( A ) and the number of rows of the matrix op( B ). K must
*             be at least zero.
*             Unchanged on exit.
*
* ALPHA - DOUBLE PRECISION.
*             On entry, ALPHA specifies the scalar alpha.
*             Unchanged on exit.
*
* A       - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*             k when TRANSA = 'N' or 'n', and is m otherwise.
*             Before entry with TRANSA = 'N' or 'n', the leading m by k
*             part of the array A must contain the matrix A, otherwise
*             the leading k by m part of the array A must contain the
*             matrix A.
*             Unchanged on exit.
*
* LDA - INTEGER.
*             On entry, LDA specifies the first dimension of A as declared
*             in the calling (sub) program. When TRANSA = 'N' or 'n' then
*             LDA must be at least max( 1, m ), otherwise LDA must be at
*             least max( 1, k ).
*             Unchanged on exit.
*
* B       - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*             n when TRANSB = 'N' or 'n', and is k otherwise.
*             Before entry with TRANSB = 'N' or 'n', the leading k by n
*             part of the array B must contain the matrix B, otherwise
*             the leading n by k part of the array B must contain the
*             matrix B.
*             Unchanged on exit.
*
* LDB - INTEGER.
*             On entry, LDB specifies the first dimension of B as declared
*             in the calling (sub) program. When TRANSB = 'N' or 'n' then
*             LDB must be at least max( 1, k ), otherwise LDB must be at
*             least max( 1, n ).
*             Unchanged on exit.
*
* BETA - DOUBLE PRECISION.
*             On entry, BETA specifies the scalar beta. When BETA is
*             supplied as zero then C need not be set on input.
*             Unchanged on exit.
*
* C       - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*             Before entry, the leading m by n part of the array C must
*             contain the matrix C, except when beta is zero, in which
*             case C need not be set on entry.
*             On exit, the array C is overwritten by the m by n matrix
*             ( alpha*op( A )*op( B ) + beta*C ).
*
* LDC - INTEGER.
*             On entry, LDC specifies the first dimension of C as declared
*             in the calling (sub) program. LDC must be at least
*             max( 1, m ).
*             Unchanged on exit.
*
*
* Level 3 Blas routine.
*
* -- Written on 8-February-1989.
*    Jack Dongarra, Argonne National Laboratory.
*    Iain Duff, AERE Harwell.
*    Jeremy Du Croz, Numerical Algorithms Group Ltd.
*    Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*    .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*    ..
*    .. External Subroutines ..
      EXTERNAL XERBLA
*    ..
*    .. Intrinsic Functions ..
      INTRINSIC MAX
*    ..
*    .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
*    ..
*    .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*    ..
*
*    Set NOTA and NOTB as true if A and B respectively are not
*    transposed and set NROWA, NCOLA and NROWB as the number of rows
*    and columns of A and the number of rows of B respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*    Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     + (.NOT.LSAME(TRANSA,'T'))) THEN
         INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     +       (.NOT.LSAME(TRANSB,'T'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGEMM ',INFO)
         RETURN
      END IF
*
*    Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     + (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*    And if alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 20 J = 1,N
                  DO 10 I = 1,M
                        C(I,J) = ZERO
10             CONTINUE
20       CONTINUE
         ELSE
            DO 40 J = 1,N
                  DO 30 I = 1,M
                        C(I,J) = BETA*C(I,J)
30             CONTINUE
40       CONTINUE
         END IF
         RETURN
      END IF
*
*    Start the operations.
*
      IF (NOTB) THEN
         IF (NOTA) THEN
*
*             Form C := alpha*A*B + beta*C.
*
            DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                        DO 50 I = 1,M
                              C(I,J) = ZERO
50                   CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                        DO 60 I = 1,M
                              C(I,J) = BETA*C(I,J)
60                   CONTINUE
                  END IF
                  DO 80 L = 1,K
                        IF (B(L,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(L,J)
                              DO 70 I = 1,M
                                 C(I,J) = C(I,J) + TEMP*A(I,L)
70                         CONTINUE
                        END IF
80             CONTINUE
90       CONTINUE
         ELSE
*
*             Form C := alpha*A'*B + beta*C
*
            DO 120 J = 1,N
                  DO 110 I = 1,M
                        TEMP = ZERO
                        DO 100 L = 1,K
                              TEMP = TEMP + A(L,I)*B(L,J)
100                   CONTINUE
                        IF (BETA.EQ.ZERO) THEN
                              C(I,J) = ALPHA*TEMP
                        ELSE
                              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                        END IF
110             CONTINUE
120       CONTINUE
         END IF
      ELSE
         IF (NOTA) THEN
*
*             Form C := alpha*A*B' + beta*C
*
            DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                        DO 130 I = 1,M
                              C(I,J) = ZERO
130                   CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                        DO 140 I = 1,M
                              C(I,J) = BETA*C(I,J)
140                   CONTINUE
                  END IF
                  DO 160 L = 1,K
                        IF (B(J,L).NE.ZERO) THEN
                              TEMP = ALPHA*B(J,L)
                              DO 150 I = 1,M
                                 C(I,J) = C(I,J) + TEMP*A(I,L)
150                         CONTINUE
                        END IF
160             CONTINUE
170       CONTINUE
         ELSE
*
*             Form C := alpha*A'*B' + beta*C
*
            DO 200 J = 1,N
                  DO 190 I = 1,M
                        TEMP = ZERO
                        DO 180 L = 1,K
                              TEMP = TEMP + A(L,I)*B(J,L)
180                   CONTINUE
                        IF (BETA.EQ.ZERO) THEN
                              C(I,J) = ALPHA*TEMP
                        ELSE
                              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                        END IF
190             CONTINUE
200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*    End of DGEMM .
*
      END
