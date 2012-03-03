C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      LOGICAL FUNCTION QSAME(Q1,Q2)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: qsame.f,v 1.1.1.1 2003/04/02 19:21:47 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     qsame -- determine the equivalence of two cartesian coordinates
C
C SYNOPSIS
      DOUBLE PRECISION Q1(3),Q2(3)
C
C DESCRIPTION
C     True if the absolute difference in any component of the cartesian
c     coordinates is less than or equal to 1.0e-6.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer I
      Double precision Z
      QSAME=.TRUE.
      DO 10 I=1,3
       Z=ABS(Q1(I)-Q2(I))
       IF(Z.GT.1.D-6)QSAME=.FALSE.
10    CONTINUE
      RETURN
      END
