
c NOTE: This routine was blatantly ripped from the SLATEC library.
c       It sorts the primary integer array and carries along the
c       operations on an auxiliary double precision array.
c       I, Anthony Yau, only modified the DSORT source making X an
c       integer array and removing options 1 and -1 from the KFLAG list
c       of possibilities. For the complete header of comments
c       pertaining to this subroutine, consult
c       http://www.netlib.org/slatec/src/dsort.f

      SUBROUTINE IDSORT (IX, DY, N, KFLAG)
C***PURPOSE  Sort an integer array and make the same interchanges in
C            an auxiliary double precision array. The array may be
c            sorted in increasing or decreasing order.
C            algorithm is used.
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION OF PARAMETERS
C
C      IX - array of values to be sorted (usually destinations)
C      DY - array to be carried along
C      N  - number of values in array IX to be sorted
C      KFLAG - control parameter
C            =  2  means sort IX in increasing order and carry DY along.
C            = -2  means sort IX in decreasing order and carry DY along.
C
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   000729  DATE MODIFIED FROM DSORT
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      INTEGER IX(*)
      DOUBLE PRECISION DY(*)
C     .. Local Scalars ..
      INTEGER T, TT
      DOUBLE PRECISION R, TY, TTY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C
C***FIRST EXECUTABLE STATEMENT IDSORT
      NN = N
      IF (NN.LT.1) THEN
         CALL XERMSG ('SLATEC', 'IDSORT',
     +      'The number of values to be sorted is not positive.', 1, 1)
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.2) THEN
         CALL XERMSG ('SLATEC', 'IDSORT',
     +      'The sort control parameter, K, is not 2 or -2.', 2, 1)
         RETURN
      ENDIF
C
C     Alter array IX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            IX(I) = -IX(I)
   10    CONTINUE
      ENDIF
C
C     Sort IX and carry DY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX(IJ)
      TY = DY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (IX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (IX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = IX(I+1)
      TY = DY(I+1)
      IF (IX(I) .LE. T) GO TO 170
      K = I
C
  180 IX(K+1) = IX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 180
      IX(K+1) = T
      DY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            IX(I) = -IX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
