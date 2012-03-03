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
      SUBROUTINE PROJPLD(NATOM,NIRREP,IORDER,ATMASS,PLD,PLDSCR,
     &                   SYOP,IPTR,NBFATM,ILCATM,SCR)
C
C THIS ROUTINE PROJECTS THE TOTALLY SYMMETRIC COMPONENT FROM A
C "PETITE" POLARIZABILITY DERIVATIVE MATRIX USING THE SYMMETRY
C OPERATIONS OF THE FULL POINT GROUP, Ajith and John 10/1998
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL YESNO
      DIMENSION PLD(3*NATOM,9),SCR(84*NATOM),PLDSCR(3*NATOM,9)
      DIMENSION SYOP(9*IORDER),IPTR(NATOM,IORDER),NBFATM(NATOM)
      DIMENSION ILCATM(NATOM),ATMASS(NATOM), POLTMP1(3,3),
     &          POLTMP2(3,3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100),IFLAGS2(500)
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      DATA ONEM /-1.0/
C
      INQUIRE(FILE='POLDER',EXIST=YESNO)
      IF(.NOT.YESNO)RETURN
C
      IF(IFLAGS(1).GE.1)THEN
       WRITE(6,1000)
      ENDIF
      ININE=9
C
      OPEN(UNIT=61,FILE='POLDER',FORM='FORMATTED',STATUS='OLD')
      REWIND(61)
      CALL ZERO(PLD,9*3*NATOM)
      DO 10 IXYZ=1,9
       IOFF=1
       READ(61,*)
       DO 11 IATOM=1,NATOM
        IF(ATMASS(IATOM).NE.ZILCH)THEN
         READ(61,'(4F20.10)')ZJUNK,(PLD(J,IXYZ),J=IOFF,IOFF+2)
        ENDIF
        IOFF=IOFF+3
11     CONTINUE
10    CONTINUE
      CLOSE(UNIT=61,STATUS='KEEP')
C      Write(6,*) "Just after Readin"
C      Call output(PLD, 1,3*NATOM, 1,9, 3*NATOM, 9, 1)
C
C GET SOME INFORMATION FROM JOBARC
C
      CALL FILTER(PLD, 27*NATOM,1.0d-8)
      CALL DGETREC(20,'JOBARC','FULLSYOP',9*IORDER,SYOP)
      CALL IGETREC(20,'JOBARC','FULLPERM',NATOM*IORDER,IPTR)
      CALL DGETREC(20,'JOBARC','ORIENTMT',ININE,SCR)
      CALL TRNOPS(SYOP,SCR,IORDER)
C
C FILL BASIS VECTOR.  SKIP DUMMY ATOMS.
C
      DO 5 IATOM=1,NATOM
       NBFATM(IATOM)=3
       ILCATM(IATOM)=3*(IATOM-1)+1
5     CONTINUE
C
      ZNORM=ONE/DFLOAT(IORDER)
C
C HALF-PROJECT CARTESIAN POLARIZABILITY DERIVATIVES
C
C LOOP OVER CARTESIAN DIRECTIONS
C
      IOFFS=54*NATOM+1
      CALL ZERO(PLDSCR,27*NATOM)
      DO 30 IOP=1,IORDER
         ITMP=IPTR(IOP,1)
         CALL XDCOPY(27*NATOM,PLD,1,SCR,1)
         IOFFA=1
         IOFFB=27*NATOM+1
         DO 20 IBAS=1,9
            CALL IMAGE(NATOM,3*NATOM,1,IOP,IPTR,NBFATM,ILCATM,
     &                 SCR(IOFFA),SCR(IOFFB),SCR(IOFFS),1,3*NATOM,
     &                 SYOP,0)
            IOFFA=IOFFA+3*NATOM
            IOFFB=IOFFB+3*NATOM
 20      CONTINUE 
C         Write(6,*) "Half Projection"
C         Call output(SCR(27*NATOM+1), 1,9, 1,3*NATOM, 9, 3*NATOM, 1)
C
         CALL TRANSP(SCR(27*NATOM+1),SCR(54*NATOM+1),9,3*NATOM)
         CALL XDCOPY(27*NATOM,SCR(54*NATOM+1),1,SCR(27*NATOM+1),1)
      
         IOFFA=27*NATOM+1
         IOFFB=1
         IPTR(IOP,1)=1
C
         DO 22 IBAS=1,3*NATOM
            CALL XDCOPY(9, SCR(IOFFA), 1, POLTMP1, 1)
            DO 23 I = 1, 3
               CALL IMAGE(1,3,1,IOP,IPTR,NBFATM,ILCATM,
     &                    POLTMP1(1,I),POLTMP2(1,I),SCR(IOFFS),1,3,
     &                    SYOP,0)
 23         CONTINUE
            CALL TRANSP(POLTMP2, POLTMP1, 3, 3)
            DO 24 I = 1, 3
               CALL IMAGE(1,3,1,IOP,IPTR,NBFATM,ILCATM,
     &                    POLTMP1(1,I),POLTMP2(1,I),SCR(IOFFS),1,3,
     &                    SYOP,0)
 24         CONTINUE
            CALL XDCOPY(9, POLTMP2, 1, SCR(IOFFA), 1)
            IOFFA=IOFFA+9
 22   CONTINUE
C
C         Write(6,*) "Full Projection-1"
C         Call output(SCR, 1,9, 1,3*NATOM, 9, 3*NATOM, 1)

         IPTR(IOP,1)=ITMP
         CALL XDAXPY(27*NATOM,ONE,SCR,1,PLDSCR,1)
 30   CONTINUE
C
C      Write(6,*) "Full Projection-2"
C      Call output(PLDSCR, 1,9,1,3*NATOM,9,3*NATOM, 1)
C
      CALL XDSCAL(27*NATOM,ZNORM,PLDSCR,1)
      CALL TRANSP(PLDSCR,SCR,3*NATOM,9)
      CALL XDAXPY(27*NATOM,ONEM,SCR,1,PLD,1)
      ILOC=ISAMAX(27*NATOM,PLD,1)
      DIFMAX=PLD(ILOC,1)
      IF(IFLAGS(1).GE.1)WRITE(6,1001)DIFMAX
      IF(DIFMAX.GT.1.D-4)THEN
       WRITE(6,1002)
      ENDIF
C      CALL TRANSP(PLDSCR,PLD,3*NATOM,9)
C
C      Write(6,*) "Full Projection-3"
C      Call output(PLDSCR, 1,3*NATOM, 1,9, 3*NATOM,9,  1)

C
      OPEN(UNIT=61,FILE='POLDER',FORM='FORMATTED',STATUS='OLD')
      REWIND(61)
      ZJUNK=0.0D0
      DO 12 IXYZ=1,9
       IOFF=1
       WRITE(61,'(I5)')IXYZ
       DO 13 IATOM=1,NATOM
        IF(ATMASS(IATOM).NE.ZILCH)THEN
         WRITE(61,'(4F20.10)')ZJUNK,(PLDSCR(J,IXYZ),J=IOFF,IOFF+2)
        ENDIF
        IOFF=IOFF+3
13     CONTINUE
12    CONTINUE
      CLOSE(UNIT=61,STATUS='KEEP')
C
      RETURN
1000  FORMAT(T3,'@PROJDLD-I, Projecting polarizability derivatives ',
     &          'on to totally symmetric subspace.')
1001  FORMAT(T3,'Largest difference between matrix elements of ',
     &          'symmetrized',/,T3,'and unsymmetrized PLD : ',F15.10,
     &          '.')
1002  FORMAT(T3,'@PROJPLD-W, The input PLD was not totally symmetric.')
      END
