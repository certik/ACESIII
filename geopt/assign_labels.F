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
      SUBROUTINE ASSIGN_LABELS (IREDUNCO, TOTREDNCO, TOTNOFBND,
     &                          TOTNOFANG, TOTNOFDIH,INTLABEL)
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C Labels all the redundant coordinates according to whether they
C are bonds,angles or dihedrals
C
#include "mxatms.par"
C
      INTEGER TOTREDNCO,TOTNOFBND,TOTNOFANG,TOTNOFDIH
      DIMENSION IREDUNCO(4,TOTREDNCO)
      CHARACTER*4 INTLABEL(TOTREDNCO)
      CHARACTER*3 LABEL(MAXREDUNCO)
C
      DATA MONE/-1/, MTWO/-2/
      do i = 1, min(9,MAXREDUNCO)
         write(label(i),'(2i1)') 0, i
      end do
      do i = 10, min(99,MAXREDUNCO)
         write(label(i),'(i2)') i
      end do
      do i = 100, MAXREDUNCO
         write(label(i),'(i3)') i
      end do
C
      DO 10 I=1,TOTNOFBND
            IF (IREDUNCO(4,I).EQ.0) THEN
                 INTLABEL(I)='R'//LABEL(I)
            END IF
 10   CONTINUE
            
      DO 20 I=TOTNOFBND+1, (TOTREDNCO-TOTNOFDIH)
            IF (IREDUNCO(4,I).LE.MONE) THEN
                 INTLABEL(I)='A'//LABEL(I-TOTNOFBND)
            END IF
 20   CONTINUE
 
      DO 30 I=TOTNOFBND+TOTNOFANG+1,TOTREDNCO
            IF (IREDUNCO(4,I).GT.0) THEN
                INTLABEL(I)='D'//LABEL(I-TOTNOFBND-TOTNOFANG)
            END IF
 30   CONTINUE
  
C      We may avoid scanning 3 times the whole totrednco simply by getting
C      TOTNOFBND,TOTNOFANG,TOTNOFDIH (1-->TOTNOFBND,TOTNOFBND+1 -->...

      RETURN
      END
   
