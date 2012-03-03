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
#include "flags.h"
      SUBROUTINE SUMMARY(SCRATCH, RFAMAT, GRD, SCALE, STPMAG,
     &                   IQFIX, NOPT, NX, NXM6, IBREAK, ICONTL,
     &                   LUOUT, NCYCLE, LUDONE, BMATRIX,
     &                   HES_INTACT, FSCR, VEC, STEP)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "cbchar.com"
#include "coord.com"
      integer ishell
      external ishell
      CHARACTER*5 a5, a5b
      CHARACTER*1 DISTANCE
      INTEGER TOTREDNCO
      DOUBLE PRECISION OPT_GEOMETRY
      LOGICAL bDonFil, skip, XYZIN, NWFINDIF, I_UNIQUE(MAXREDUNCO),
     &        Hessian_calcs, Geomopt_numrcl 
C
      DIMENSION SCRATCH(NX*NX), RFAMAT(NOPT, NOPT), VEC(NOPT),
     &          GRD(NXM6), BMATRIX(NXM6*NX), 
     &          HES_INTACT(NX*NX), FSCR(NX*NXM6),
     &          IREDUNCO(4, MAXREDUNCO), IQFIX(3*NATOMS, 3*NATOMS),
     &          STEP(NXM6) 
C
      COMMON /INPTYP/ XYZIN,NWFINDIF

C
C
      Common /Orient/ Orient(3,3)
c
#include "machsp.com"
C#define MONSTER_FLAGS /* use ioppar */
#include "jodaflags.com"
C
C Do some Printing
C
      WRITE(LUOUT, 88)STPMAG, SCALE
 88   FORMAT(T3,' Summary of Optimization Cycle: ',/,T3,
     &       ' The maximum unscaled ',
     &       'step is: ',F10.5,'.',/,T3,' Scale factor set to: ',
     &        F8.5,'.')
C
      WRITE(LUOUT, 93)
93    FORMAT(T3,' Forces are in hartree/bohr and hartree/radian.',
     &       /,T3,' Parameter values are in Angstroms and degrees.')
C
      WRITE(LUOUT, 81)
      WRITE(LuOut, 881)
      WRITE(LuOut, 81)
881   FORMAT(T3,'Parameter',T17,'dV/dR',T33,'Step',T47,'Rold',
     &       T63,'Rnew')
C
      DO 15 I=1, NOPT
C
         DO 23 J =1, I - 1
            IF (VARNAM(ISQUASH(NOPTI(I))) .EQ.
     &          VARNAM(ISQUASH(NOPTI(J)))) GOTO 15
 23      CONTINUE

         CNF = 180.0D0/DACOS(-1.D0)
         IF (.NOT. XYZIN) THEN
            IF (MOD(ISQUASH(NOPTI(I)),3) .EQ. 1) CNF = 0.529177249D0
         ELSE
            DISTANCE = VARNAM(NOPTI(I))
            IF (DISTANCE .EQ. "R")  CNF = 0.529177249D0
         ENDIF
         ZOP = SCRATCH(I) - SCRATCH(NOPT+I)

C
         WRITE(LUOUT, 84)VARNAM(ISQUASH(NOPTI(I))), GRD(NOPTI(I)),
     &                   SCRATCH(NOPT+I)*CNF, ZOP*CNF,
     &                   SCRATCH(I)*CNF
C
 15   CONTINUE
C
 84   FORMAT(T5,A,T10,4(F15.10,1X))
C
C Generate statistics based on gradients
C
      CALL ZERO(RFAMAT, NOPT*NOPT)
C
      DO 432 I =1, NOPT
         RFAMAT(I, 1) = GRD(NOPTI(I))
 432  CONTINUE
C
      CALL VSTAT(RFAMAT, SCRATCH(2*NOPT + 1), NOPT)
C
      WRITE(LUOUT,82) SCRATCH(2*NOPT + 2),SCRATCH(2*NOPT + 5)
 82   FORMAT(74('-'),/,T3,'Minimum force: ',F12.9,' / RMS force: ',
     &       F12.9)
 81   FORMAT(74('-'))
C
C see if the calcualtion has converged
C
      IF ( (SCRATCH(2*NOPT + 5) .LT. 10.D0**(-1*ICONTL) ) .AND.
     &     (IBREAK .EQ. 0)                                      ) THEN
C
         WRITE(6,1321)
 1321    FORMAT(80('-'))
         WRITE(6,1325)10.D0**(-1*ICONTL)
 1325    FORMAT(T3,' RMS gradient is below ',E10.5,'.')
         WRITE(6,1322)
 1322    FORMAT(T3,' Convergence criterion satisfied.  Optimization ',
     &      'completed.')
         WRITE(6,1321)

         CALL ADM

cmn
c write optimized values of coordinates in table
c for use with vibrational freq. calculation.
c
c The coordinates printed correspond to the geometry that generated the
c last set of energetics, i.e., the geometry whose gradient is below the
c convergence criterion.
C
         IF (Iflags2(h_IFLAGS2_geom_opt) .GE. 3) THEN
C
            write(6,*)
            write(6,*) '  Summary of optimized internal coordinates'
            write(6,*) '  (Angstroms and degrees)'
            write(6,*) (Iflags2(h_IFLAGS2_geom_opt))
            CALL IGETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)
            CALL IGETREC(20, 'JOBARC', 'UNIQUEDF', TOTREDNCO, I_UNIQUE)
            CALL IGETREC(20, 'JOBARC', 'CONTEVIT', 4*TOTREDNCO,
     &                  IREDUNCO)

            DO I=1,NOPT
            if (I_UNIQUE(NOPTI(I))) then

c            o remove spaces
               a5b=VARNAM(ISQUASH(NOPTI(I)))
               a5='     '
               l=0
               do k=1,5
                  if (a5b(k:k) .ne. ' ') then
                     l = l+1
                     a5(l:l) = a5b(k:k)
                  endif
               enddo

c >R01  :  2  1         =     1.0796431125
c >R03  :  4  1         =     1.3252448582
c >A01  :  2  1  3      =   116.8470115097
c >A02  :  2  1  4      =   121.5764336125
c >D01  :  5  4  1  2   =     0.0000000000
c >D02  :  5  4  1  3   =   180.0000000000

               INDEX = ISQUASH(NOPTI(I))
               DISTANCE = VARNAM(INDEX)
               OPT_GEOMETRY = SCRATCH(I) - SCRATCH(NOPT+I)
               IF (DISTANCE .EQ. "R") THEN
                  CNF = 0.529177249D0
                  WRITE(LUOUT, 9) a5, IREDUNCO(1, INDEX),
     &                            IREDUNCO(2, INDEX),
     &                            OPT_GEOMETRY*CNF
 9                FORMAT(A5,T6,':',2(1X,I2),9X,"=",F17.10)
               ELSE IF (DISTANCE .EQ. "A") THEN
                  CNF=180.D0/DACOS(-1.D0)
                  WRITE(LUOUT, 99) a5, IREDUNCO(1, INDEX),
     &                             IREDUNCO(2, INDEX),
     &                             IREDUNCO(3, INDEX),
     &                             OPT_GEOMETRY*CNF
 99               FORMAT(A5,T6,':',3(1X,I2),6X,'=',F17.10)
               ELSE
                  CNF=180.D0/DACOS(-1.D0)
                  WRITE(LUOUT, 999) a5, IREDUNCO(1, INDEX),
     &                              IREDUNCO(2, INDEX),
     &                              IREDUNCO(3, INDEX),
     &                              IREDUNCO(4, INDEX),
     &                              OPT_GEOMETRY*CNF
 999              FORMAT(A5,T6,':',4(1X,I2),3X,"=",F17.10)
               END IF

c           end if (I_UNIQUE(NOPTI(I)))
            end if
c           END DO I=1,NOPT
            END DO

c         o print out the Cartesian coordinate matrix (broken)
       write(6,*)
       write(6,*) ' Summary of optimized Cartesian coordinates (Bohr) '
       write(6,*)
C 
            CALL DGETREC(20, 'JOBARC', 'COORD_OP',3*NATOMS, Q)
C
            IOFF=1
            DO I=1,NATOMS
 1000          FORMAT(T6,A,T10,F14.8,T25,F14.8,T40,F14.8)
               WRITE(6,1000)ZSYM(I)(1:2),(Q(J),J=IOFF,IOFF+2)
               IOFF=IOFF+3
            END DO

c ELSE IF (INTERNAL)

         ELSE IF (Iflags2(h_IFLAGS2_geom_opt) .EQ. 1) THEN
C
            write(6,*)
            write(6,*) '  Summary of optimized internal coordinates'
            write(6,*) '  (Angstroms and degrees)'
            write(6,*) (Iflags2(h_IFLAGS2_geom_opt))

            DO I=1,NOPT
               OPT_GEOMETRY = SCRATCH(I) - SCRATCH(NOPT+I)
               skip = .false.
               do J=1,I-1
                  skip = skip.or.
     &                   (VARNAM(ISQUASH(NOPTI(I))).EQ.
     &                    VARNAM(ISQUASH(NOPTI(J)))    )
               end do
            if (.not.skip) then
c            o remove spaces
               a5b=VARNAM(ISQUASH(NOPTI(I)))
               a5='     '
               l=0
               do k=1,5
                  if (a5b(k:k) .ne. ' ') then
                     l = l+1
                     a5(l:l) = a5b(k:k)
                  endif
               enddo
               IF (MOD(ISQUASH(NOPTI(I)),3).EQ.1) THEN
                  CNF=0.529177249D0
               ELSE
                  CNF=180.D0/DACOS(-1.D0)
               END IF
               WRITE(LuOut,184) A5, OPT_GEOMETRY*CNF
 184           FORMAT(A5,T6,'=',F17.10)
c           end if (.not.skip)
            end if
c           END DO I=1,NOPT
            END DO

c ELES IF (Pure Cartesian)
C
         ELSE IF (Iflags2(h_IFLAGS2_geom_opt) .EQ. 2) THEN
C
c    o print out the Cartesian coordinate matrix (broken)
            write(6,*) 
            write(6,"(a,a)") ' Summary of optimized Cartesian',
     &                          ' coordinates (Bohr) '
            write(6,*) 
C      
            CALL DGETREC(20, 'JOBARC', 'COORD_OP',3*NATOMS, Q)
C           
            IOFF=1
            DO I=1,NATOMS
               WRITE(6,1000)ZSYM(I)(1:2),(Q(J),J=IOFF,IOFF+2)
               IOFF=IOFF+3
            END DO
C
C ENDIF (Iflags2(h_IFLAGS2_geom_opt))
C
         ENDIF
         WRITE(LUOUT,*)
cmn end

C      o Let's also update the archive file.
C
C Unpak the scratch array and regenerate the full R vector. Also
C unpak the step vector to avoid any problems with Hessian updates  
C in the next cycle.
C
         DO J = 1, NOPT
C
            INDEX = ISQUASH(NOPTI(J))
            DISTANCE = VARNAM(INDEX)
            IF (DISTANCE .EQ. "R") THEN
               CNF = 0.529177249D0
            ELSE IF (DISTANCE .EQ. "A") THEN
               CNF=180.D0/DACOS(-1.D0)
            ELSE
               CNF=180.D0/DACOS(-1.D0)
            END IF

            R(NOPTI(J))    = (SCRATCH(J) - SCRATCH(NOPT+J))
            STEP(NOPTI(J)) = SCRATCH(J + NOPT)
C
            DO K=1, NEQ(NOPTI(J))
C
               STEP(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
               R(IQFIX(NOPTI(J), K))    = (SCRATCH(J) - 
     &                                     SCRATCH(NOPT+J))
            ENDDO
         ENDDO
C
C The following calls do anything useful only when the optimizations
C are done numerically, and diretcly followed by frequency calculation
C with analytical or numerical second derivatives.
C What we need is the Cartesian coordinates of the optimized geometry
C in COORD jobarc record. In the case of analytical optimizations the
C previous step correspond to the optimized geometry and this call
C is simply a duplication. However, that is not the case for 
C numerical optimizations (the record contains the coordinates of the
C last dispalced point). Note to Anthony: this will effect the
C tests that use numerical gradients. 
C 
         Call igetrec(20, 'JOBARC', 'IFLAGS  ', 100, Iflags)
         Call igetrec(20, 'JOBARC', 'IFLAGS2 ', 500, Iflags2)
         Hessian_calcs   =  (Iflags(h_IFLAGS_vib) .NE. 0)
         Geomopt_numrcl  =  (Iflags2(h_IFLAGS2_grad_calc) .EQ. 2
     &                      .AND. Iflags2(h_IFLAGS2_geom_opt)
     &                      .NE. 0)
c
#ifdef _DEBUG_LVL0
         Write(6,*) "The internal coord:"
         Write(6,10) (R(I), I=1, NXM6)
         Write(6,*)
         Write(6,*) "The step size:"
         Write(6,*)
         Write(6,10)   (STEP(I), I=1, NXM6)
         Write(6,*)
         Print*, "The optimization cycle:", NCYCLE
         Print*, "Joda flag, grad_calc and geoom_opt",
     &            Iflags2(h_IFLAGS2_grad_calc),
     &            Iflags2(h_IFLAGS2_geom_opt) 
         Print*, "Internal flags; Hessian & geo. optimization:",
     &            Hessian_calcs,  Geomopt_numrcl 
  10     Format(5(1X,F10.6))
#endif
c
         If (Hessian_calcs .OR. Geomopt_numrcl) Then 
            Call GMetry(.TRUE., .FALSE.)
            IStart = 1
            Iscr   = Istart + NX*3
            IQtmp  = Iscr   + NX
            IQnew  = IQtmp  + NX

            Call Symmetry(Scratch(Istart), Scratch(Iscr), 
     &                    Scratch(IQtmp), .FALSE.)
            Call dputrec(20, 'JOBARC', 'COORD   ', NX, Q)
            Call dputrec(20, 'JOBARC', 'CORD_INT', NX, R)
#ifdef _DEBUG_LVL0
         Write(6,*)
         Print*, "The variable from summary" 
         Write(6,20) "The Cartesian coord:", (Q(i), i=1, NX)
   20    Format(A20, 3(1X, F10.6))

#endif
         Endif
c
c         call iputrec(20, 'JOBARC', 'CYCL2CON', 1, NCYCLE)
c         call archive(0.0D0, GRD, HES_INTACT, STEP, 0, VEC)
c
C Marshall's request to provide the frequencies of the updated Hessian
C might be a good idea. The frequencies of the updated Hessian at the
C last step might give a clue about the character of the stationary
C point. If this approx. test passes, then the expensive "exact"
C frequency calculation should be done. Ajith Perera, 05/2005.

         CALL GEN_APPROX_FREQ(SCRATCH, BMATRIX, FSCR, HES_INTACT)

c      o signal the end of the calculation
         call igetrec(1,'JOBARC','JODADONE',1,i)
         if (i.ne.0) then
            print *, '@SUMMARY: JODADONE is already true.'
            print *, '          Something is very wrong.'
         end if
         call iputrec(1,'JOBARC','JODADONE',1,1)
         call iputrec(1,'JOBARC','HAVEGEOM',1,1)

c      o remove the backup directory
         call rm_backup

c         CALL ACES_JA_FIN
c         STOP

      END IF

      RETURN
      END

