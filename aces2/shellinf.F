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
      SUBROUTINE SHELLINF(NORBIT,NTOTATOM,NTOTSHEL,INUMSHEL,ISHLOFF,
     &                   ISHELTP,ISHELSZ)
C
C GENERATE SHELL INFORMATION FOR ATOMS AS ORDERED IN THE ZMAT
C FILE.
C
CJDW 10/28/96. Modifications to check value of IMOL and to check we do
C              not exceed fixed dimensions. For dummy atoms IMOL is set
C              to 999 and we have troubles. Thanks to Roger Edberg of ANU
C              for finding this error.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION INUMSHEL(NORBIT),ISHLOFF(NORBIT)
      DIMENSION ISHELTP(NTOTSHEL),ISHELSZ(NTOTSHEL)
C
#include <mxatms.par>
#include <baslims.par>
C
      DIMENSION IZMT2MOL(mxatms)
C
C I am going to set these two arrays to the maximum number of
C atoms times maximum number of shells. Ajith Perera, 11/07.
C
      DIMENSION ISHELTP2(mxatms*mxshel),ISHELSZ2(mxatms*mxshel)
C
      IF(NTOTATOM.GT.200)THEN
       WRITE(6,1000)
       CALL ERREX
      ENDIF
C
      CALL IGETREC(20,'JOBARC','ZMAT2MOL',NTOTATOM,IZMT2MOL)
#ifdef _DEBUG_LVLM1
      Print*, "ZMAT2MOL", (IZMT2MOL(I), I=1, NTOTATOM)
      Print*, "ISHELTP", (ISHELTP(I), I=1, NTOTSHEL)
      Print*, "ISHELSZ", (ISHELSZ(I), I=1, NTOTSHEL)
#endif 
C
C LOOP OVER ATOMS IN ZMAT
C
      IOFFZMAT=1
      DO 10 IZMAT=1,NTOTATOM
C
C GET MOL FILE POSITION
C
       IMOL=IZMT2MOL(IZMAT)
C
       IF(IMOL.GE.1 .AND. IMOL.LE.NORBIT)THEN
C
C GET SHELL INFORMATION FOR THIS ATOM
C
        NSHELL=INUMSHEL(IMOL)
        IOFFMOL=ISHLOFF(IMOL)
        if (IOFFZMAT+NSHELL.gt.500) then
           print *, '@SHELLINF: Assertion failed.'
           print *, '           maximum number of shells = 500'
           print *, '           require at least ',IOFFZMAT+NSHELL
           call errex
        end if
        CALL ICOPY(NSHELL,ISHELTP(IOFFMOL),1,ISHELTP2(IOFFZMAT),1)
        CALL ICOPY(NSHELL,ISHELSZ(IOFFMOL),1,ISHELSZ2(IOFFZMAT),1)
        IOFFZMAT=IOFFZMAT+NSHELL
C
       ELSE
C
        IF(IMOL.NE.999)THEN
         WRITE(6,1010) IMOL
         CALL ERREX
        ENDIF
C
       ENDIF
C
10    CONTINUE
C
      NSHELTOT=IOFFZMAT-1
C
      IF(NSHELTOT.GT.500)THEN
       WRITE(6,1020) NSHELTOT
       CALL ERREX
      ENDIF
C
#ifdef _DEBUG_LVLM1
      Print*, "+++++++++++++++++", NSHELTOT
      Print*, "ISHELTP", (ISHELTP2(I), I=1, NTOTSHEL)
      Print*, "ISHELSZ", (ISHELSZ2(I), I=1, NTOTSHEL)                          
#endif
C
      CALL IPUTREC(20,'JOBARC','FULSHLNM',1,NSHELTOT)
      CALL IPUTREC(20,'JOBARC','FULSHLTP',NSHELTOT,ISHELTP2)
      CALL IPUTREC(20,'JOBARC','FULSHLSZ',NSHELTOT,ISHELSZ2)
C
#ifdef _DEBUG_LVLM1
      Print*, "++++++++++++++++++"
      Print*, "ISHELTP", (ISHELTP2(I), I=1, NTOTSHEL)
      Print*, "ISHELSZ", (ISHELSZ2(I), I=1, NTOTSHEL)
#endif
C
      RETURN
 1000 FORMAT(' @SHELLINF-F, Too many atoms (over 200) ',I10)
 1010 FORMAT(' @SHELLINF-F, Invalid value of IMOL ',I10)
 1020 FORMAT(' @SHELLINF-F, Too many shells (over 500) ',I10)
      END
