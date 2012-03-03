
C THIS ROUTINE PRINTS OUT A SUMMARY OF THE SYMMETRY ADAPTED COORDINATES.

c This routine was essentially the same as the old DUMPCORD routine but had
c a few extra print statements. The BFLAG argument controls that extra
c printing. Whether or not we need both options is for someone else to decide,
c but functionally, the program is the same as before minus the existence
c of DUMPCORD.

c INPUT
c integer NATOM
c integer NIRREP
c integer NIRREPC
c char*4  DOIT
c logical BFLAG

c OUTPUT
c double  DSCR(*)
c integer ISCR(*)
c char*8  LABEL(NIRREP+NIRREPC)

c RECORDS
c get DOIT//'SYMQ'
c get DOIT//'SYQT'
c get DOIT//'LABL'
c get 'NUMVIBRT'
c get 'LINEAR  '
c get 'SBGRPSYM'
c get 'COMPLABL'
c get DOIT//'PTGP'
c get 'COMPPTGP'

      SUBROUTINE PRTCOORD(NATOM,NIRREP,NIRREPC,DSCR,ISCR,
     &                    LABEL,DOIT,BFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION DSCR(*),ISCR(*)
      CHARACTER*8 LABEL(NIRREP+NIRREPC)
      CHARACTER*4 DOIT
      LOGICAL BFLAG

      CHARACTER*4 FPGRP,CPGRP

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "control.com"

      NSIZE=3*NATOM

      CALL GETREC(20,'JOBARC',DOIT//'SYMQ',NSIZE*NSIZE*IINTFP,DSCR)
      CALL GETREC(20,'JOBARC',DOIT//'SYQT',NSIZE,ISCR)
      CALL GETREC(20,'JOBARC',DOIT//'LABL',NIRREP*IINTFP,LABEL)
      CALL GETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      CALL GETREC(20,'JOBARC','LINEAR  ',1,ILINEAR)
      IF (BFLAG) THEN
      CALL GETREC(20,'JOBARC','SBGRPSYM',NSIZE,ISCR(NSIZE+1))
      CALL GETREC(20,'JOBARC','COMPLABL',NIRREPC*IINTFP,LABEL(NIRREP+1))
      CALL GETCREC(20,'JOBARC',DOIT//'PTGP',4,FPGRP)
      CALL GETCREC(20,'JOBARC','COMPPTGP',  4,CPGRP)
      END IF

      IOFF=1
      WRITE(6,1000)
1000  FORMAT(T3,' *** TABLE OF SYMMETRY COORDINATES *** ')
      IF (BFLAG) THEN
      WRITE(6,1010)
1010  FORMAT(T3,'                                       ')
      IF (ROTPROJ) THEN
         WRITE(6,1020)
1020     FORMAT(T3,' Eckart rotations and translations projected ',
     &             'from vibrational coordinates')
      ELSE
         WRITE(6,1021)
1021     FORMAT(T3,' Translations projected from vibrational ',
     &             'coordinates')
      END IF
      WRITE(6,1030)FPGRP,CPGRP
1030  FORMAT(T3,' Point group used : ',A,' ; Abelian subgroup : ',A)
      WRITE(6,1031)DOIT
1031  FORMAT(T3,'Coordinate type : ',A)
c     END IF (BFLAG)
      END IF
      DO IMODE=1,NMODE
         IF (BFLAG) THEN
         WRITE(6,1003)IMODE,LABEL(ISCR(IMODE)),
     &                LABEL(NIRREP+ISCR(NSIZE+IMODE))
1003     FORMAT(T3,' Mode : ',I3,' Symmetry : ',A,' Subgroup symm : ',A)
         ELSE
         WRITE(6,1001)IMODE,LABEL(ISCR(IMODE))
1001     FORMAT(T3,' Mode : ',I3,' Symmetry : ',A)
         END IF
         WRITE(6,1004)
1004     FORMAT(T11,'ATOM',T26,'X',T40,'Y',T54,'Z')
         DO IATOM=1,NATOM
            WRITE(6,1005)IATOM,(DSCR(J),J=IOFF,IOFF+2)
            IOFF=IOFF+3
1005        FORMAT(T11,I3,T19,F13.10,T33,F13.10,T47,F13.10)
         END DO
         WRITE(6,*)
      END DO

      IF (ROTPROJ) THEN
         WRITE(6,2000)
2000     FORMAT(T3,'@PRTCOORD: Rotational coordinates (x,y,z)')
         DO IXYZ=1,3-ILINEAR
            WRITE(6,1004)
            DO IATOM=1,NATOM
               WRITE(6,1005)IATOM,(DSCR(J),J=IOFF,IOFF+2)
               IOFF=IOFF+3
            END DO
         END DO
      END IF

      WRITE(6,3000)
3000  FORMAT(T3,'@PRTCOORD: Translational coordinates (x,y,z)')
      DO IXYZ=1,3
         WRITE(6,1004)
         DO IATOM=1,NATOM
            WRITE(6,1005)IATOM,(DSCR(J),J=IOFF,IOFF+2)
            IOFF=IOFF+3
         END DO
      END DO

      RETURN
      END

