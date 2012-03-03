      SUBROUTINE SP2SY(N, A, LDA, UPLO, IERR)
C      
C     PURPOSE
C     
C     To expand a symmetric matrix stored in packed storage mode to
C     standard (quadratic) storage mode using the same work array.       
C
C     ARGUMENT LIST
C       ARGUMENTS IN
C
C          N - INTEGER.
C              The order of the matrix A.
C      
C          A - DOUBLE PRECISION array of DIMENSION at least LDA*N.
C              The first  N*(N+1)/2  elements of this array contain the
C              symmetric matrix A in (upper or lower) packed storage
C              mode.   
C              NOTE that this array is overwritten.
C
C        LDA - INTEGER.
C              The leading dimension of array A. After calling SP2SY, A 
C              is stored like an array declared  A(LDA,N) containing a
C              symmetric matrix in its leading N-by-N part.
C      
C     ARGUMENTS OUT
C    
C          A - DOUBLE PRECISION array of DIMENSION at least LDA*N.
C              This array contains the matrix A stored like an array
C              declared  A(LDA,N) containing a symmetric matrix in its 
C              leading N-by-N part.
C
C     MODE PARAMETERS
C
C       UPLO - CHARACTER*1.
C              Indicates whether on input, A is stored in packed lower
C              or packed upper storage mode, i.e. if the columns of the
C              lower or upper trinagle of the matrix A are stored
C              consecutively.
C              UPLO = 'U' or 'u' : Upper triangle is stored in A.
C              UPLO = 'L' or 'l' : Lower triangle is stored in A.
C              NOTE that if UPLO contains any other character, A is
C              treated as packed in lower symmetric storage mode.
C      
C     ERROR INDICATOR
C     
C       IERR - INTEGER. 
C              Unless the routine detects an error (see next section),
C              IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C      IERR = 1 : (N .LT. 1) or (LDA. LT. N).
C
C     METHOD
C
C     Let
C               ( A(1,1) ... A(1,N) )
C           A = (  ...        ...   ).
C               ( A(N,1) ... A(N,N) ) 
C
C     On input, A is stored in either upper or lower packed storage
C     mode as follows:
C     If UPLO = 'U' or 'u', then the columns of the upper triangle of A
C     are contained in A consecutively, i.e.      
C     
C       A = (A(1,1), A(1,2), A(2,2), A(1,3), ..., A(3,3), A(1,4), ...,
C            A(1,N), ..., A(N,N)).      
C
C     If UPLO = 'L' or 'l', then the columns of the lower triangle of A
C     are contained in A consecutively, i.e.      
C     
C       A = (A(1,1), A(2,1), ..., A(N,1), A(2,2), A(3,2), ..., A(N,2),
C            A(3,3), ..., A(N,N)).      
C   
C     On output, A contains the full symmetric matrix A, i.e. A may be
C     treated as a 2-dimensional array declared A(LDA,N) with the
C     leading N-by-N part containing the symmetric matrix A.
C      
C     REFERENCE
C                 
C     None.
C      
C     CONTRIBUTOR      
C
C     Peter Benner (Technische Universitaet Chemnitz-Zwickau, FRG)
C
C     REVISIONS
C
C     February 21, 1995.
C            
C***********************************************************************
C     
C     .. Scalar Arguments ..      
      INTEGER          IERR, LDA, N
      CHARACTER        UPLO
C     
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
C
C     .. Local Scalars ..
      INTEGER          I, SPPOS, SYPOS
C
C     .. External Subroutines ..
C     . BLAS, LAPACK .
      EXTERNAL         XCOPY
C     
C     .. Executable Statements ..
C
      IERR = 0
      IF ((N .LT. 1) .OR. (N .GT. LDA))  IERR = 1
      IF (IERR .NE. 0)  GOTO 2001
C      
      IF ((UPLO.EQ.'U').OR.(UPLO.EQ.'u')) THEN
        SPPOS = N*(N+1)/2 + 1
        SYPOS = N*LDA + 1
        DO 100  I = N, 2, -1
          SPPOS = SPPOS - I
          SYPOS = SYPOS - LDA
          CALL XCOPY(I, A(SPPOS), -1, A(SYPOS), -1)
100     CONTINUE       
        SYPOS = 1
        DO 150  I = N-1, 1, -1
          CALL XCOPY(I, A(SYPOS+LDA), LDA, A(SYPOS+1), 1)
          SYPOS = SYPOS + LDA + 1
150     CONTINUE
C       
      ELSE
        SPPOS  = N*(N+1)/2
        SYPOS  = (N-1)*LDA + N
        DO 200  I = 1, N - 1
          CALL XCOPY(I,   A(SPPOS),  -1, A(SYPOS),      -1)
          CALL XCOPY(I-1, A(SYPOS+1), 1, A(SYPOS+LDA), LDA)
          SPPOS = SPPOS - I - 1
          SYPOS = SYPOS - LDA - 1
200     CONTINUE
        CALL XCOPY(N-1, A(2), 1, A(LDA+1), LDA)
      END IF
C      
2001  CONTINUE
      RETURN
C *** Last Line of SP2SY      
      END
      

         
