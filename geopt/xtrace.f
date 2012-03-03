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
      DOUBLE PRECISION FUNCTION XTRACE(A,N)
      IMPLICIT NONE
      INTEGER N, INC, I, ADD
      DOUBLE PRECISION A(N*N)
      XTRACE = A(1)
      INC = N + 1
      ADD = 1
      DO I = 1, N-1
         ADD = ADD + INC
         XTRACE = XTRACE + A(ADD)
      END DO
      RETURN
      END

