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
      SUBROUTINE GETCOVLRADI(IATOMICNMBER, SMOFCOVRADI, NRATMS)
C
C Compute sum of covalent radii for each unquie pair of atoms. The
C covalent radii are look up from the table, and was provided by 
C Stefan Fau (1999/2000). 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MAXPERTABLE =86, SCALE = 1.30D00)
      DIMENSION COVLNTRADI(MAXPERTABLE), IATOMICNMBER(NRATMS), 
     &          SMOFCOVRADI(NRATMS, NRATMS)
C
      DATA COVLNTRADI /0.32, 0.93, 1.23, 0.90, 0.80, 0.77, 0.74, 0.73, 
     &                  0.72, 0.71, 1.54, 1.36, 1.18, 1.11, 1.06, 1.02,
     &                  0.99, 0.98, 2.03, 1.74, 1.44, 1.32, 1.22, 1.18,
     &                  1.17, 1.17, 1.16, 1.15, 1.17, 1.25, 1.26, 1.22,
     &                  1.20, 1.16, 1.14, 1.12, 2.16, 1.91, 1.62, 1.45,
     &                  1.34, 1.30, 1.27, 1.25, 1.25, 1.28, 1.34, 1.48,    
     &                  1.44, 1.41, 1.40, 1.36, 1.33, 1.31, 2.35, 1.98,
     &                  1.69, 1.65, 1.65, 1.64, 1.64, 1.68, 1.85, 1.61,
     &                  1.59, 1.59, 1.57, 1.57, 1.56, 1.70, 1.56, 1.44,
     &                  1.34, 1.30, 1.28, 1.26, 1.27, 1.30, 1.34, 1.49,      
     &                  1.48, 1.47, 1.46, 1.46, 1.45, 2.14/

      DATA AMTOBOHR /1.8897265/
C
      CALL ZERO(SMOFCOVRADI, NRATMS*NRATMS)
C
C Compute the Scaled covalents radii, The Scaling Factor of 1.3 is provided
C by Stefan Fau.
C
      DO 10 IBNDS = 2, NRATMS
C         DO 20 JBNDS = JBNDS, NRATMS
          DO 20 JBNDS = 1, (IBNDS - 1)
            IF (IBNDS .NE. JBNDS) THEN
               SMOFCOVRADI(JBNDS, IBNDS) = 
     &                                (COVLNTRADI(IATOMICNMBER(IBNDS))
     &                              + COVLNTRADI(IATOMICNMBER(JBNDS)))
     &                              * SCALE * AMTOBOHR
#ifdef _DEBUG_LVL_M1
            Print*, "The covalent radii assignments"
            Print*,  IBNDS, JBNDS, IATOMICNMBER(IBNDS),
     &               IATOMICNMBER(JBNDS)
            Print*,  SMOFCOVRADI(JBNDS, IBNDS)*0.529177249d0

#endif
            ENDIF
 20      CONTINUE
 10   CONTINUE
C
      RETURN
      END
