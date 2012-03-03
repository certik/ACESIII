
C     A subroutine that calculates the determinant and inverse of
C          a matrix, as well as solving systems of linear equations.
C     Martin J. McBride.  11/25/85.
C     General Electric CRD, Information System Operation.
C
C  AB  - The matrix that need to inverted (or solved); It is done in place!
C   N  - The order of the matrix AB
C  LDA - Leading Dimension of the matrix AB
C  SCRATCH - scratch are of dimension n to hold pivot indices (integer)
C  DET - The determinat corresponding to AB matrix (output)
C  EPS - The tolerence to check the linear dependencies (1.0D-8 is recommended)
C  M   - Control whether you invert the matrix of solve linear equations AB*X=C
C      - 0 invert, 1 would solve linear eqns.
C MODE - is not used

C  THE COMMENTS ABOVE REFER TO AN OLD "SCIPORT" VERSION OF THIS ROUTINE,
C  PRESUMABLY WRITTEN BY MR. MC BRIDE.  UNFORTUNATELY, IT SEEMS THAT
C  HE IS A BAD PROGRAMMER AND THIS ROUTINE DID NOT WORK.  MINV IS NOW
C  A SIMPLE INTERFACE TO LINPACK ROUTINES WHICH *WORK*

c#include "aces.h"

      SUBROUTINE MINV(AB,N,LDA,SCRATCH,DET,EPS,M,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION AB(LDA,1),SCRATCH(1)

      INFO = 0
      CALL DGETRF(N,N,AB,LDA,SCRATCH,INFO)
      IF (INFO.NE.0) THEN
	IF (INFO.GT.0) THEN
	  WRITE(6,*) "@MINV Singular Matrix"
        ELSE
          WRITE(6,*) "@MINV dgetrf argument ",-INFO," is illegal."
        ENDIF
        CALL ERREX
      ENDIF
      IF(M.EQ.1) CALL DGETRS('N',N,1,AB,LDA,SCRATCH,AB(1,N+1),LDA,INFO)
      IF(M.GT.1)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@MINV-F, Only one system can be solved at a time.')
       CALL ERREX
      ENDIF

      IOFF=1+N
      IF(M.EQ.0) CALL DGETRI(N,AB,LDA,SCRATCH,SCRATCH(IOFF),N,INFO)
C
C THIS IS WRONG, BUT WHO CARES? I care, Ajith 12/99
C
c     DET=X(1)*10**X(2)
C      IF (DET .LE. EPS) Then
      IF (INFO.NE.0) THEN
        IF (INFO.GT.0) THEN 
          WRITE(6,*) "@MINV Singular Matrix"
        ELSE
          WRITE(6,*) "@MINV dgetri argument ",-INFO," is illegal."
        ENDIF
        CALL ERREX
      ENDIF
      RETURN
      END
