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
      SUBROUTINE JSD1(JSD,JSD2,NREAL)
C
C COMPUTES THE SD CONTRIBUTION TOn INDIRECT SPIN-SPIN
C NMR COUPLING CONSTANTS
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL YESNO,PRINT
      DOUBLE PRECISION JSD,JSD2
C
#include "mxatms.par"
      DIMENSION JSD(3*NREAL,3),JSD2(NREAL,NREAL)
      DIMENSION NORD(MXATMS*2)
C
#include "cbchar.com"
#include "coord.com"
C
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
C
      DATA THIRD /0.333333333333333D0/
      DATA HALF  /0.5D0/
C
      DO 2000 ITYPE=1,2
C
      IF(ITYPE.EQ.1) THEN
       INQUIRE(FILE='JSD',EXIST=YESNO)
      ELSE
       INQUIRE(FILE='JSDSCF',EXIST=YESNO)
      ENDIF
C
      PRINT=.FALSE.
      IF(IFLAGS(1).GE.10)PRINT=.TRUE.
C
      IF(YESNO)THEN
C
C GET INFORMATION ABOUT ATOMIC ORDERING
C
       CALL IGETREC(20,'JOBARC','DUMSTRIP',NATOMS,NORD)
       CALL IGETREC(20,'JOBARC','MAP2ZMAT',NATOMS,NORD(NATOMS+1))

       IF(ITYPE.EQ.1) THEN
        OPEN(UNIT=90,FILE='JSD',FORM='FORMATTED',STATUS='OLD')
       ELSE
        OPEN(UNIT=90,FILE='JSDSCF',FORM='FORMATTED',STATUS='OLD')
       ENDIF
C
C READ IN THE CHEMICAL SHIELDING TENSOR
C
       JBOT=1
       DO 10 I=1,NREAL
C
        IZ=NORD(NATOMS+I)
        IZ0=NORD(IZ)
C
        READ(90,'((3F20.10))')(JSD(J,1),J=1,NREAL*3)
        READ(90,'((3F20.10))')(JSD(J,2),J=1,NREAL*3)
        READ(90,'((3F20.10))')(JSD(J,3),J=1,NREAL*3)
C
C COMPUTE FOR ALL ATOMS ISOTROPIC COUPLINGS
C
        I0=0
        DO 11 J=1,NREAL
C
        JZ=NORD(NATOMS+J)
        JZ0=NORD(JZ)
C
         JSD2(IZ0,JZ0)=THIRD*(JSD(I0+1,1)+JSD(I0+2,2)+JSD(I0+3,3))
         I0=I0+3
11      CONTINUE
C
10     CONTINUE
C
       WRITE(6,1000)
       IF(ITYPE.EQ.2) THEN
        WRITE(6,1001)
       ELSE
        IF(IFLAGS(2).EQ.0) THEN
         WRITE(6,1001)
        ELSE IF(IFLAGS(2).EQ.1) THEN
         WRITE(6,10001)
        ELSE IF(IFLAGS(2).EQ.2) THEN
         WRITE(6,10011)
         CALL ERREX
        ELSE IF(IFLAGS(2).EQ.10) THEN
         WRITE(6,10012)
        ELSE
         CALL ERREX
        ENDIF
       ENDIF
       WRITE(6,1000)
C
       JBOT=1
       NTIMES=1+(NREAL-1)/5
C
       DO 20 ICOUNT=1,NTIMES
        WRITE(6,*)
        WRITE(6,142)(ZSYM(ICN),ICN=JBOT,MIN(NREAL-1,JBOT+4))
        WRITE(6,144)(ICN,ICN=JBOT,MIN(NREAL-1,JBOT+4))
        DO 12 I=JBOT+1,NREAL
12      WRITE(6,143)ZSYM(I),I,(JSD2(I,J),J=JBOT
     & ,MIN(I-1,JBOT+4))
        JBOT=JBOT+5
20      CONTINUE
C
142    FORMAT(17X,A3,4(9X,A3))
143    FORMAT(T3,A3,'[',I2,']',5(2X,F10.5))
144    FORMAT(16X,:'[',I2,']',4(8X,:'[',I2,']'))

       WRITE(6,1000)
       CLOSE(UNIT=90,STATUS='KEEP')
      ELSE
       IF(ITYPE.EQ.1) THEN    
        WRITE(6,1004)
        CALL ERREX
       ENDIF
      ENDIF
C
2000  CONTINUE
      RETURN
C
1000  FORMAT(T3,67('-'))
1001  FORMAT(T5,'SCF Spin-Dipole NMR Coupling Constants (in Hz)')
10001 FORMAT(T3,'MBPT(2) Spin-Dipole NMR Coupling Constants (in Hz)')
10011 FORMAT(T3,'MBPT(3) Spin-Dipole NMR Coupling Constants (in Hz)')
10012 FORMAT(T3,'CCSD-EH Spin-Dipole NMR Coupling Constants (in Hz)')
1003  FORMAT(T7,I3,T20,I3,T33,F10.3,T51,F10.3)
1004  FORMAT(T3,'@JSD1-F, File JSD not found on disk.')
2001  FORMAT(T3,I3,3F20.10)
C
      END
