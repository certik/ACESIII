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

C Computes chemical shieldings and anisotropies from the chemical shielding
C tensors.

C Also Print the total shielding tensor and the orientation of the 
C shielding tensor. Added by Ajith 04/96 in accordence with many pepole
C who extensively use ACESII NMR capabilities (Don Chesnut, Jim Franz, 
C Dewey Barich and Jim Haw...).

      SUBROUTINE NMR1(CSHIFT,SHIFT,ANISOT, RANGE, NREAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL YESNO,PRINT
      DOUBLE PRECISION DWORK(8), DTMP(3)
      INTEGER ADD_I, ADD_J, NDX_I, NDX_J

#include "mxatms.par"
      DIMENSION CSHIFT(9,NREAL),SHIFT(NREAL),ANISOT(NREAL),
     &          RANGE(NREAL)
      DIMENSION NORD(MXATMS*2),BUF(9), TMP(3,3)
C     
#include "cbchar.com"
#include "coord.com"
C
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
C     
      DATA THIRD /0.333333333333333D0/
      DATA HALF  /0.5D0/
C     
C Print the orientation of the molecule again for convenient 
C interpretation of NMR shielding tensors.
C
      IONE = 1 
C
      CALL IGETREC(20, 'JOBARC', 'NREALATM', IONE, NREAL)
      CALL DUMPCORD(NREAL, Q, IATNUM)
C
      DO 2000 ITYPE = 1, 2
C     
         IF(ITYPE.EQ.1) THEN
            INQUIRE(FILE='CSHIFT', EXIST=YESNO)
         ELSE
            INQUIRE(FILE='CSHIFTSCF', EXIST=YESNO)
         ENDIF
C     
         PRINT=.FALSE.
         IF(IFLAGS(1).GE.10)PRINT=.TRUE.
C     
         IF (YESNO) THEN
C     
C GET INFORMATION ABOUT ATOMIC ORDERING
C     
            CALL IGETREC(20,'JOBARC','DUMSTRIP',NATOMS,NORD)
            CALL IGETREC(20,'JOBARC','MAP2ZMAT',NATOMS,NORD(NATOMS+1))

            IF(ITYPE.EQ.1) THEN
               OPEN(UNIT=90,FILE='CSHIFT',FORM='FORMATTED',
     &              STATUS='OLD')
            ELSE
               OPEN(UNIT=90,FILE='CSHIFTSCF',FORM='FORMATTED',
     &              STATUS='OLD')
            ENDIF
C     
C READ IN THE CHEMICAL SHIELDING TENSOR
C     
            JBOT=1
C
            DO 10 I = 1, NREAL
C
               JZ  = NORD(NATOMS + I)
               JZ0 = NORD(JZ)
C
               READ(90, '((3F20.10))') (CSHIFT(J,JZ0),J = 1, 9)
C     
C SYMMETRIZE MATRIX
C     
               CALL SYMMET2(CSHIFT(1,JZ0), 3)
C     
C COMPUTE ISOTROPIC CHEMICAL SHIELDING ( 1/3 TR)
C     
               SHIFT(JZ0) = THIRD*(CSHIFT(1,JZ0) + CSHIFT(5,JZ0) + 
     &                       CSHIFT(9,JZ0))
C     
C DIAGONALIZE TO GET ANISOTROPY
C
cYAU 20001127
c old:
c               CALL EIG2(BUF,CSHIFT(1,JZ0),3,3,-1)
c new:
               CALL DSYEV('V','L',3,CSHIFT(1,JZ0),3,BUF,DWORK,8,INFO)
               IF (INFO.NE.0) THEN
                  WRITE(*,*)
     & '@NMR1: There was a problem diagonalizing the matrix.'
                  CALL ERREX
               END IF
               DO NDX_I = 1, 2
                  DO NDX_J = NDX_I+1, 3
                     IF ( BUF(NDX_I) .LT. BUF(NDX_J) ) then
                        DTMP(1)    = BUF(NDX_I)
                        BUF(NDX_I) = BUF(NDX_J)
                        BUF(NDX_J) = DTMP(1)
                        ADD_I   = 3*NDX_I
                        DTMP(1) = CSHIFT(ADD_I-2,JZ0)
                        DTMP(2) = CSHIFT(ADD_I-1,JZ0)
                        DTMP(3) = CSHIFT(ADD_I  ,JZ0)
                        ADD_J               = 3*NDX_J
                        CSHIFT(ADD_I-2,JZ0) = CSHIFT(ADD_J-2,JZ0)
                        CSHIFT(ADD_I-1,JZ0) = CSHIFT(ADD_J-1,JZ0)
                        CSHIFT(ADD_I  ,JZ0) = CSHIFT(ADD_J  ,JZ0)
                        CSHIFT(ADD_J-2,JZ0) = DTMP(1)
                        CSHIFT(ADD_J-1,JZ0) = DTMP(2)
                        CSHIFT(ADD_J  ,JZ0) = DTMP(3)
                     END IF
                  END DO
               END DO
c :end
cYAU
C
               IF (ITYPE .EQ. 2) THEN
                  WRITE (6, 2005) 
               ELSE
                  WRITE (6, 2006)
               ENDIF           
C
C Write the principal values of the shielding tensor.
C     
               CALL ZERO(TMP, 9)
               DO 5 III = 1, 3
                  TMP(III, III) = BUF(III)
 5             CONTINUE
C
               CALL CSHIFTPRI(TMP(1,1), JZ0)
C     
               IF (ITYPE .EQ. 2) THEN
                  WRITE (6, 2008) 
               ELSE
                  WRITE (6, 2009)
               ENDIF           
C     
C Write orientation of the shielding tensor.
C
               CALL CSHIFTPRI(CSHIFT(1, JZ0), JZ0)
C     
               IF (PRINT) WRITE(6,2001) JZ0, BUF
C     
               ANISOT(JZ0) = BUF(1) - HALF*(BUF(2)+BUF(3))
               RANGE(JZ0) = BUF(1) - BUF(3)
C     
 10         CONTINUE
C     
            WRITE (6, 999)
C     
            WRITE (6,1000)
C
            IF (ITYPE.EQ.2) THEN
               WRITE(6,1001)
            ELSE
               IF(IFLAGS(2).EQ.0) THEN
                  WRITE(6,1001)
               ELSE IF(IFLAGS(2).EQ.1) THEN
                  WRITE(6,10001)
               ELSE IF(IFLAGS(2).EQ.2) THEN
                  WRITE(6,10011)
               ELSE IF(IFLAGS(2).EQ.3) THEN
                  WRITE(6,10013)
               ELSE IF(IFLAGS(2).EQ.4) THEN
                  WRITE(6,10014)
               ELSE IF(IFLAGS(2).EQ.8) THEN
                  WRITE(6,10015)
               ELSE IF(IFLAGS(2).EQ.23) THEN
                  WRITE(6,10016)
               ELSE IF(IFLAGS(2).EQ.10) THEN
                  WRITE(6,10017)
               ELSE
                  CALL ERREX
               ENDIF
            ENDIF

            WRITE(6,1000)
            WRITE(6,1002)
            WRITE(6,1000)
C
            IOFF=0
            DO 20 I=1,NREAL
C     
C SOMEWHAT A HASSLE TO GET THE CORRECT ATOMIC NUMBER, HOWEVER IT WORKS
C     
               WRITE(6,1003)I, ZSYM(I), SHIFT(I), ANISOT(I), RANGE(I)
C
 20         CONTINUE
C
            WRITE(6,1000)
            CLOSE(UNIT=90, STATUS='KEEP')
C
         ELSE
C
            IF (ITYPE .EQ. 1) THEN    
               WRITE(6,1004)
               CALL ERREX
            ENDIF
C
         ENDIF
C     
 2000 CONTINUE
C      
      RETURN
C     
 999  FORMAT(/)
 1000 FORMAT(T3,72('-'))
 1001 FORMAT(T6,'SCF Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10001 FORMAT(T4,'MBPT(2) Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10011 FORMAT(T4,'MBPT(3) Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10013 FORMAT(T4,'SDQ-MBPT(4) Nuclear Magnetic Resonance Shieldings ',
     &   'and Anisotropies')
10014 FORMAT(T4,'MBPT(4) Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10015 FORMAT(T4,'CCD Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10016 FORMAT(T4,'QCISD Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10017 FORMAT(T4,'CCSD Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
10012 FORMAT(T4,'CCSD-EH Nuclear Magnetic Resonance Shieldings and ',
     &   'Anisotropies')
 1002 FORMAT(T5,'Z-matrix',T19,'Atomic',T31,'Chemical Shielding',
     &   T53,'Anisotropy',T69, 'Range',/,T6,'Center',T19,'Symbol',
     &   T37,'(ppm)',T55,'(ppm)', T69, '(ppm)')
 1003 FORMAT(T7,I2,T21,A3,T33,F10.3,T51,F10.3, T65, F10.3)
 1004 FORMAT(T3,'@NMR1-F, File CSHIFT not found on disk.')
 2001 FORMAT(T3,I3,3F20.10)
C     
 2005 FORMAT(/,'  Principal values of the shielding tensor (SCF):',/)
 2006 FORMAT(/,'  Principal values of the shielding tensor (MBPT(2)):',
     &   /)
 2008 FORMAT(/,'  Orientation of the SCF shielding tensor:',/)
 2009 FORMAT(/,'  Orientation of the MBPT(2) shielding tensor:',/)
C     
      END
