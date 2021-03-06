*DECK DGEMM
      SUBROUTINE DGEMM_PRT (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B,
     $   LDB, BETA, C, LDC)
C***BEGIN PROLOGUE  DGEMM
C***PURPOSE  Perform one of the matrix-matrix operations.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B6
C***TYPE      DOUBLE PRECISION (SGEMM-S, DGEMM-D, CGEMM-C)
C***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S. (NAG)
C***DESCRIPTION
C
C  DGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X',
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = A'.
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = B'.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
C                 A set of level 3 basic linear algebra subprograms.
C                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
C***ROUTINES CALLED  XERBLA
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  DGEMM
C     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, INDEX
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      LOGICAL            NOTA, NOTB
      CHARACTER*6        NCTnct
      DATA NCTnct /'NCTnct'/
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C***FIRST EXECUTABLE STATEMENT  DGEMM
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
C     and  columns of  A  and the  number of  rows  of  B  respectively.
C
      NOTA  = (TRANSA.EQ.'N').OR.(TRANSA.EQ.'n')
      NOTB  = (TRANSB.EQ.'N').OR.(TRANSB.EQ.'n')
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF(      INDEX(NCTnct,TRANSA).EQ.0 ) THEN
         INFO = 1
      ELSE IF( INDEX(NCTnct,TRANSB).EQ.0 ) THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     And if  alpha.eq.zero.
C
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  WRITE(*,*) 'C(',I,',',J,') = ',ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  WRITE(*,*) 'C(',I,',',J,') *= ',BETA
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF( NOTB )THEN
         IF( NOTA )THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     WRITE(*,*) 'C(',I,',',J,') = ',ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     WRITE(*,*) 'C(',I,',',J,') *= ',BETA
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = B( L, J )
                     DO 70, I = 1, M
                        WRITE(*,*) 'C(',I,',',J,') += ',ALPHA,' * ',
     &                             TEMP,' * ',A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 120, J = 1, N
               DO 110, I = 1, M
                  IF( BETA.EQ.ZERO )THEN
                     WRITE(*,*) 'C(',I,',',J,') = ',ZERO
                  ELSE
                     WRITE(*,*) 'C(',I,',',J,') *= ',BETA
                  END IF
                  DO 100, L = 1, K
                     WRITE(*,*) 'C(',I,',',J,') += ',ALPHA,' * ',
     &                          A( L, I ),' * ',B( L, J )
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
C
C           Form  C := alpha*A*B' + beta*C
C
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     WRITE(*,*) 'C(',I,',',J,') = ',ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     WRITE(*,*) 'C(',I,',',J,') *= ',BETA
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = B( J, L )
                     DO 150, I = 1, M
                        WRITE(*,*) 'C(',I,',',J,') += ',ALPHA,' * ',
     &                             TEMP,' * ',A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 200, J = 1, N
               DO 190, I = 1, M
                  IF( BETA.EQ.ZERO )THEN
                     WRITE(*,*) 'C(',I,',',J,') = ',ZERO
                  ELSE
                     WRITE(*,*) 'C(',I,',',J,') *= ',BETA
                  END IF
                  DO 180, L = 1, K
                     WRITE(*,*) 'C(',I,',',J,') += ',ALPHA,' * ',
     &                          A( L, I ),' * ',B( J, L )
  180             CONTINUE
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of DGEMM .
C
      END
