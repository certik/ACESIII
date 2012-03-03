      SUBROUTINE SY2SP(N, A, LDA, UPLOI, UPLOO, IERR)
C      
C     PURPOSE
C     
C     To store a symmetric or triangular matrix given in packed storage
C     mode using the same work array.        
C
C     ARGUMENT LIST
C       ARGUMENTS IN
C
C          N - INTEGER.
C              The order of the matrix A.
C      
C          A - DOUBLE PRECISION array of DIMENSION (LDA,N).
C              On input, the leading N-by-N part contains the upper or
C              lower triangular part of A depending upon UPLOI (see
C              section MODE PARAMETERS).
C              NOTE that this array is overwritten.
C
C        LDA - INTEGER.
C              The leading dimension of array A. 
C      
C     ARGUMENTS OUT
C    
C          A - DOUBLE PRECISION array of DIMENSION (LDA,N).
C              On output, A contains an  N*(N+1)/2  vector containing
C              the symmetric matrix A in packed storage mode (lower or
C              upper, depending upon UPLOO (see section MODE
C              PARAMETERS)). 
C
C     MODE PARAMETERS
C
C      UPLOI - CHARACTER*1.
C              Indicates whether on input, A contains the upper or lower
C              triangular part of A.
C              UPLOI = 'U' or 'u' : Upper triangle is stored in A.
C              UPLOI = 'L' or 'l' : Lower triangle is stored in A.
C              NOTE that if UPLOI contains any other character, A is
C              considered as containing the lower triangular part of A.
C      
C      UPLOO - CHARACTER*1.
C              Indicates whether on output, A is packed in upper or
C              lower symmetric storage mode.
C              UPLOO = 'U' or 'u' : A is stored in packed upper
C                                   storage mode.
C              UPLOO = 'L' or 'l' : A is stored in packed lower
C                                   storage mode.
C              NOTE that if UPLOO contains any other character, A is
C              packed to lower storage mode.
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
C     On input, A contains either the upper or lower triangular part of
C     A depending upon the input parameter UPLOI. This triangular part
C     of A then is packed to either lower or upper storage mode as
C     follows:
C     If UPLOO = 'U' or 'u', then the columns of the upper triangle of A
C     are contained in A consecutively, i.e.      
C     
C       A = (A(1,1), A(1,2), A(2,2), A(1,3), ..., A(3,3), A(1,4), ...,
C            A(1,N), ..., A(N,N)).      
C
C     If UPLOO = 'L' or 'l', then the columns of the lower triangle of A
C     are contained in A consecutively, i.e.      
C     
C       A = (A(1,1), A(2,1), ..., A(N,1), A(2,2), A(3,2), ..., A(N,2),
C            A(3,3), ..., A(N,N)).      
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
      CHARACTER        UPLOI, UPLOO
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
      IF (((UPLOI.EQ.'U').OR.(UPLOI.EQ.'u')) .AND.
     &    ((UPLOO.EQ.'L').OR.(UPLOO.EQ.'l'))) THEN
        SYPOS = 1
        DO 100  I = N-1, 1, -1
          CALL XCOPY(I, A(SYPOS+LDA), LDA, A(SYPOS+1), 1)
          SYPOS = SYPOS + LDA + 1
100     CONTINUE
      ELSE IF (((UPLOI.EQ.'L').OR.(UPLOI.EQ.'l')) .AND.
     &         ((UPLOO.EQ.'U').OR.(UPLOO.EQ.'u'))) THEN
        SYPOS = 1
        DO 200  I = N-1, 1, -1
          CALL XCOPY(I, A(SYPOS+1), 1, A(SYPOS+LDA), LDA)
          SYPOS = SYPOS + LDA + 1
200     CONTINUE
      ENDIF
C      
      IF ((UPLOO.EQ.'U').OR.(UPLOO.EQ.'u')) THEN
        SPPOS = 2
        SYPOS = LDA + 1
        DO 300  I = 2, N
          CALL XCOPY(I, A(SYPOS), 1, A(SPPOS), 1)
          SPPOS = SPPOS + I
          SYPOS = SYPOS + LDA
300     CONTINUE       
C         
      ELSE
        SPPOS  = N + 1
        SYPOS  = LDA + 2
        DO 400  I = N-1, 1, -1
          CALL XCOPY(I, A(SYPOS), 1, A(SPPOS), 1)
          SPPOS = SPPOS + I 
          SYPOS = SYPOS + LDA + 1
400     CONTINUE
      END IF
C      
2001  CONTINUE
      RETURN
C *** Last Line of SY2SP      
      END
      

         
