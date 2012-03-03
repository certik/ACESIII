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
      SUBROUTINE ASSIGN_DIHLS(IBNDTO, IREDUNCO, TOTREDNCO, TOTNOFBND, 
     &                        TOTNOFANG, NRATMS, TOTNOFDIH, MAXREDUNCO)
C
C This routine setup the connectivity array to define dihedral angles. This
C array is dimensioned to (4, MAXREDUNCO).
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER TOTREDNCO, TOTNOFBND, TOTNOFANG, TOTNOFDIH
      LOGICAL YESLINR1, YESLINR2
C
      DIMENSION IREDUNCO(4, MAXREDUNCO), IBNDTO(NRATMS, NRATMS)
C
      DATA IZERO /0/
C

#ifdef _DEBUG_LVLM1
      Write(6,*)
      Print*, "MAXREDUNCO, TOTNOFBND, TOTREDNCO"
      Print*, MAXREDUNCO, TOTNOFBND, TOTREDNCO
      Write(6,*)
      Print*, "The redundent bond and angle coor. assignments"
      Print*,  ((Iredunco(i, j), i=1,4), j=1, TOTREDNCO)
      Write(6,*)
#endif

      TOTNOFDIH = IZERO
      
      DO 10 ICON2 = 2, NRATMS
         DO 20 ICON3 = 1, (ICON2 - 1)
C
            IF (IBNDTO(ICON3, ICON2) .GT. 0) THEN
               DO 30 IBNDS = 1, TOTNOFBND
                  
                  IF (IREDUNCO(2, IBNDS) .GT. IZERO) THEN
C
                     ICON1 = IZERO
C
                     IF (IREDUNCO(1, IBNDS) .EQ. ICON2) 
     &                    ICON1 = IREDUNCO(2, IBNDS)
                     IF (IREDUNCO(2, IBNDS) .EQ. ICON2) 
     &                    ICON1 = IREDUNCO(1, IBNDS)
C
                     IF (ICON1 .NE. IZERO .AND. ICON1 .NE. ICON3) THEN
                        
                        DO 40 JBNDS =1 , TOTNOFBND
                           
                           IF (IREDUNCO(2, JBNDS) .GT. IZERO) THEN

                              ICON4 = IZERO

                              IF (IREDUNCO(1, JBNDS) .EQ. ICON3) 
     &                             ICON4 = IREDUNCO(2, JBNDS)
                              IF (IREDUNCO(2, JBNDS) .EQ. ICON3) 
     &                             ICON4 = IREDUNCO(1, JBNDS)
C
                              IF (ICON4 .NE. IZERO .AND. ICON4 .NE.
     &                            ICON1 .AND. ICON4 .NE. ICON2) THEN
C
                                 YESLINR1 = (IBNDTO(ICON3, ICON1) .EQ. 
     &                                      -99 .OR. 
     &                                      (IBNDTO(ICON3, ICON1) .NE.
     &                                       IZERO .AND.
     &                                       IBNDTO(ICON3, ICON2) .EQ.
     &                                      -99))
                                 YESLINR2 = (IBNDTO(ICON2, ICON4) .EQ. 
     &                                       -99 .OR. 
     &                                       IBNDTO(ICON2, ICON4) .NE.
     &                                       IZERO .AND.
     &                                       IBNDTO(ICON3, ICON2) .EQ.
     &                                       -99)
                              IF (ICON1 .LT. ICON3)
     &                           YESLINR1 = (IBNDTO(ICON1, ICON3) .EQ. 
     &                                      -99 .OR. 
     &                                      (IBNDTO(ICON1, ICON3) .NE.
     &                                       IZERO .AND.
     &                                       IBNDTO(ICON3, ICON2) .EQ.
     &                                      -99))
                              IF (ICON4 .LT. ICON2)
     &                           YESLINR2 = (IBNDTO(ICON4, ICON2) .EQ. 
     &                                       -99 .OR. 
     &                                       IBNDTO(ICON4, ICON2) .NE.
     &                                       IZERO .AND.
     &                                       IBNDTO(ICON3, ICON2) .EQ.
     &                                       -99)

                                 IF (.NOT. YESLINR1 .AND. .NOT. 
     &                               YESLINR2) THEN
C
                                    TOTNOFDIH = TOTNOFDIH + 1
                                    IREDUNCO(1, TOTREDNCO + TOTNOFDIH) =
     &                                       ICON1
                                    IREDUNCO(2, TOTREDNCO + TOTNOFDIH) =
     &                                       ICON2
                                    IREDUNCO(3, TOTREDNCO + TOTNOFDIH) =
     &                                       ICON3
                                    IREDUNCO(4, TOTREDNCO + TOTNOFDIH) =
     &                                       ICON4
                                  ENDIF
                                 ENDIF
C                              ENDIF
                           ENDIF
C
 40                     CONTINUE
C
                     ENDIF
                  ENDIF
C     
 30            CONTINUE
C
            ENDIF
C
 20      CONTINUE
 10   CONTINUE
C
      TOTREDNCO = TOTREDNCO + TOTNOFDIH
C
      IF (TOTREDNCO .GE. MAXREDUNCO) THEN
        WRITE(6, 99) 
 99     FORMAT(T3,'Maximum No. of redundent coordinates allowed',
     &         'exceeded.')
        print *,'Required:',TOTREDNCO,' Current Maximum:',MAXREDUNCO
        Write(6,*)
        CALL ERREX
      ENDIF
C
      RETURN
      END
