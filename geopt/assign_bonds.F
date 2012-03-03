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
      SUBROUTINE ASSIGN_BONDS(IBNDTO, IREDUNCO, TOTREDNCO, TOTNOFBND,
     &                        NRATMS, MAXREDUNCO)
C
C This routine setup the connectivity array to define bonds. This
C array is dimensioned to (4, MAXREDUNCO).
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER TOTREDNCO, TOTNOFBND
      DIMENSION IREDUNCO(4, MAXREDUNCO), IBNDTO(NRATMS, NRATMS)
C
      DATA IZERO /0/

      TOTNOFBND = IZERO
C
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Print*, "The connectivity assignments in assign bonds"
      Print*, ((ibndto(j, i), j=1, nratms),I=1, nratms)
      Write(6,*)
#endif
C
      DO 10 IBNDS = 2, NRATMS
         DO 20 JBNDS = 1, IBNDS - 1
            
            IF (IBNDTO(JBNDS, IBNDS) .NE. 0) THEN
C
               TOTNOFBND = TOTNOFBND + 1
C
               IREDUNCO(1, TOTREDNCO + TOTNOFBND) = IBNDS
               IREDUNCO(2, TOTREDNCO + TOTNOFBND) = JBNDS
               IREDUNCO(3, TOTREDNCO + TOTNOFBND) = IZERO
               IREDUNCO(4, TOTREDNCO + TOTNOFBND) = IZERO
C
            ENDIF
C
 20      CONTINUE
 10   CONTINUE
C
      TOTREDNCO = TOTREDNCO + TOTNOFBND
C
      IF (TOTREDNCO .GE. MAXREDUNCO) THEN
        WRITE(6, 99) 
 99     FORMAT(T3,'Maximum No. of redundent coordinates allowed',
     &         '/', 'exceeded')
        print *,'Required:',totrednco,' Current maximum: ',MAXREDUNCO
        CALL ERREX
      ENDIF
C
      RETURN
      END




