C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      LOGICAL FUNCTION QSAME(Q1,Q2)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: qsame.f,v 1.3 2010/06/30 15:58:20 ponton Exp $
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
