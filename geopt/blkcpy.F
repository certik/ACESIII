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

C THIS ROUTINE COPIES AN NROWFRM BY NCOLFRM BLOCK MATRIX (MATFRM)
C INTO AN NROWTAR BY NCOLTAR TARGET MATRIX (MATTAR) SUCH THAT 
C THE (1,1) ELEMENT IN MATFRM BECOMES THE IROWTAR,ICOLTAR ELEMENT 
C IN MATTAR.

      SUBROUTINE BLKCPY(MATFRM,NROWFRM,NCOLFRM,MATTAR,NROWTAR,NCOLTAR,
     &                  IROWTAR,ICOLTAR)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION MATTAR(NROWTAR*NCOLTAR),MATFRM(NROWFRM,NCOLFRM)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#ifdef _ASSERT
      if (ncolfrm.lt.0) then
         print *, '@BLKCPY: Assertion failed.'
         print *, '         ncolfrm = ',ncolfrm
         call errex
      end if
#endif /* _ASSERT */
      ISTART = IROWTAR + (ICOLTAR-1)*NROWTAR
      DO ICOL = 1, NCOLFRM
         CALL XCOPY(NROWFRM,MATFRM(1,ICOL),1,MATTAR(ISTART),1)
         ISTART = ISTART + NROWTAR
      END DO
      RETURN
      END
