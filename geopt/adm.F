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

C PRINTS LOWER DIAGONAL ATOMIC DISTANCE MATRIX FOR ACES2 PROGRAM SYSTEM.

      SUBROUTINE ADM
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
      PARAMETER (BA = 0.529177249d0)
      PARAMETER (NCOL = 6)
c      PARAMETER (TOL = 0.5)
      PARAMETER (TOL = 0.2)

C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
 
      DIMENSION V1(3),V2(3)

cJDW 1/6/98
c   Scratch array for holding up to five distances. Introduced
c   because previous code did not work on Suns with f90.
      DIMENSION SCR(5)
      INTEGER ITMP

C     Labels used throughout the program:
C     ZSYM    Atomic symbol given for each line of the Z-matrix
C     VARNAM  Symbols of all variable parameters
C     PARNAM  Symbols of all variables *and* (fixed) parameters
C
#include "cbchar.com"
#include "coord.com"
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT

      IQ(I)=3*I-2

c ----------------------------------------------------------------------

c      write(6,*) '  @ADM-I, Coordinates : '
c      write(6,'(3F15.10)') Q(1),Q(2),Q(3)
c      write(6,'(3F15.10)') Q(4),Q(5),Q(6)
c      write(6,'(3F15.10)') Q(7),Q(8),Q(9)

      iexit = 0

      WRITE(6,'(t3,a)') ' Interatomic distance matrix (Angstroms) '

      JBOT = 1
      NTIMES = 1 + (NATOMS-1)/5
      DO ICOUNT = 1, NTIMES
         WRITE(6,*)
         WRITE(6,142) (ZSYM(ICN),ICN=JBOT,MIN(NATOMS,JBOT+4))
 142     FORMAT(17X,A,4(9X,A))
         WRITE(6,144) (ICN,ICN=JBOT,MIN(NATOMS,JBOT+4))
 144     FORMAT(16X,:'[',I2,']',4(8X,:'[',I2,']'))
         DO I = JBOT, NATOMS
c      write(6,*) jbot,i,natoms,iq(1),iq(2),iq(3)
c      write(6,*) ' jbot,i,min(i,jbot+4) ',jbot,i,min(i,jbot+4)
c      write(6,*) '  @ADM-I, Coordinates : '
c      write(6,'(3F15.10)') Q(1),Q(2),Q(3)
c      write(6,'(3F15.10)') Q(4),Q(5),Q(6)
c      write(6,'(3F15.10)') Q(7),Q(8),Q(9)
            ITMP = MIN(I,JBOT+4)-JBOT+1
            IF (ITMP.GT.5) THEN
               WRITE(6,*) '@ADM: Loop limit exceeded. ',I,JBOT,ITMP
            END IF
            DO J = 1, ITMP
               SCR(J) = BA * DIST(Q(IQ(I)),Q(IQ(J+JBOT-1)))
            END DO
cJDW - old
c   This led to problems with f90 on Suns. The old code did not
c   work, but the new code seems to. Essential difference is we
c   calculate the MIN(I,JBOT+4)-JBOT+1 distances before the
c   write statement.
c            WRITE(6,143) ZSYM(I),I,(BA*DIST(Q(IQ(I)),Q(IQ(J))),J=JBOT,
c     &      MIN(I,JBOT+4))
cJDW - new
            WRITE(6,143) ZSYM(I),I,(SCR(J),J=1,ITMP)
cJDW - end
 143        FORMAT(T3,A,'[',I2,']',5(2X,F10.5))
         END DO
         JBOT = JBOT + 5
      END DO

c CHECK DISTANCES TO SEE IF ANY ARE TOO SHORT. CALL EXIT IF BELOW
c A CERTAIN TOLERANCE.

       dismax = 0.0
       DO 12 I = 1, NATOMS
       DO 12 J = I+1, NATOMS
       DIS=DIST(Q(IQ(I)),Q(IQ(J)))*BA
c-pr
c  get the longest distance, for a cavity radius in RFSCF calculations
       dismax = max(dis,dismax)
c-pr
       IF (DIS.LT.TOL) THEN
c         IF ONE IS A DUMMY (or a ghost), GO ON.
          IF(ATMASS(I).LT.1.D-3.OR.ATMASS(J).LT.1.D-3)GOTO 12
          IF(ATMASS(I).ge.100d0.OR.ATMASS(J).ge.100d0)GOTO 12
          WRITE(6,231)I,J,DIS,TOL
231       FORMAT(T3,' Atoms ',i2,' and ',i2,' are too close.',
     &/,T3,' Distance of ',f6.4,' Angstroms is below threshold of '
     &,f6.4,'.')
          IEXIT = 1
       END IF
c-pr
       dismax = 0.5 + 0.5*dismax
       call dputrec(20,'JOBARC','CAVITY',1,dismax)
c-pr
12     CONTINUE
C
C CRANK OUT ANGLES IF IPRNT IS SET TO A SUFFICIENTLY HIGH VALUE.
C
       IF (IPRNT.GE.600) THEN
       WRITE(6,702)
 702   FORMAT(/,T3,' Interatomic angles (degrees) ',/)
C
C LOOP OVER ALL ATOMS, PUTTING EACH AT THE ORIGIN OF A X-Y-Z ANGLE
C   PAIR.  INCLUDE THE DUMMY ATOMS TO MAKE CERTAIN USERS HAPPY.
C
       ANGTOL = 3.D0/BA
       ICNT = 0
       DO 15 I = 1, NATOMS
        DO 205 J = 1, NATOMS
         IF (J.EQ.I) GOTO 205
         DO K = 1, J-1
          IF (K.NE.I.AND.K.NE.J) THEN
C
C CHECK TO MAKE SURE THAT D(I,J) AND D(I,K) ARE BOTH LESS THAN 3 ANGSTROMS
C IF EITHER J OR K IS TOO FAR AWAY, GO ON.
C
          IF (DIST(Q(IQ(I)),Q(IQ(J))).LE.ANGTOL.AND.
     &        DIST(Q(IQ(I)),Q(IQ(K))).LE.ANGTOL) THEN
          CALL VEC(Q(IQ(I)),Q(IQ(J)),V1,1)
          CALL VEC(Q(IQ(I)),Q(IQ(K)),V2,1)
          Z=ANGLE(V1,V2,3)
          ICNT=ICNT+1
          IF(MOD(ICNT,2).EQ.1)THEN
           WRITE(6,400)ZSYM(K),K,ZSYM(I),I,ZSYM(J),J,Z
          ELSE
           WRITE(6,401)ZSYM(K),K,ZSYM(I),I,ZSYM(J),J,Z
          ENDIF
 400      FORMAT(/,T3,A2,'[',I2,']-',A2,'[',I2,']-',A2,'[',I2,']',
     &           2X,F9.5)
 401      FORMAT(T3,A2,'[',I2,']-',A2,'[',I2,']-',A2,'[',I2,']',
     &           2X,F9.5)
          END IF
          END IF
          END DO
 205     CONTINUE
 15    CONTINUE

       IF (MOD(ICNT,2).EQ.1) THEN
          WRITE(6,402)ICNT
 402   FORMAT(/ ,T3,I5,' interatomic angles printed.')
       ELSE
          WRITE(6,403)ICNT
 403   FORMAT(//,T3,I5,' interatomic angles printed.')
       END IF

       END IF

       if (iexit.eq.1) then
          write (6,*) '@ADM: Inspect distance matrix and correct input.'
          call errex
       end if

       RETURN
       END

