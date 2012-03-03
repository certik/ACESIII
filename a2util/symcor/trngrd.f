
C THIS ROUTINE TRANSFORMS THE SYMMETRY COORDINATE GRADIENT VECTOR
C TO THE CARTESIAN REPRESENTATION.

c INPUT
c integer NATOM           : number of atoms
c double  SYMGRD(3*NATOM) : symmetry coordinate gradient
c integer NDSCR           : the amount of double scratch at DSCR
c char*4  TYPE            : (FULL|COMP) point group
c logical PRINTQ          : verbose printing flag

c OUTPUT
c double CARTGRD(3*NATOM) : cartesian gradient
c double DSCR(NDSCR)      : double scratch

c RECORDS
c get TYPE//'SYMQ'
c get 'ATOMMASS'

      SUBROUTINE TRNGRD(NATOM,SYMGRD,CARTGRD,DSCR,NDSCR,TYPE,PRINTQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION DSCR(NDSCR)
      DOUBLE PRECISION SYMGRD(3*NATOM),CARTGRD(3*NATOM)
      CHARACTER*4 TYPE
      LOGICAL PRINTQ

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      NSIZE=3*NATOM

      IF (NDSCR.LT.NSIZE*NSIZE) THEN
         print *, '@TRNGRD: Insufficient memory to load ',TYPE,'SYMQ'
         print *, '         have ',NDSCR,' doubles'
         print *, '         need ',NSIZE*NSIZE,' doubles'
         call aces_exit(1)
      END IF

c   o transform to mass-weighted cartesian coordinates
      CALL GETREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE*IINTFP,DSCR)
      CALL XGEMM('N','N',NSIZE,1,NSIZE,
     &           1.d0,DSCR,   NSIZE,
     &                SYMGRD, NSIZE,
     &           0.d0,CARTGRD,NSIZE)

c   o remove mass weighting
      CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,DSCR)
      IOFF=1
      DO IATOM=1,NATOM
         X=SQRT(DSCR(IATOM))
         CARTGRD(IOFF+0)=X*CARTGRD(IOFF+0)
         CARTGRD(IOFF+1)=X*CARTGRD(IOFF+1)
         CARTGRD(IOFF+2)=X*CARTGRD(IOFF+2)
         IOFF=IOFF+3
      END DO

      IF (PRINTQ) THEN
         write(6,*)' Gradient before transformation:'
         write(6,'((3f20.10))')(symgrd(i),i=1,nsize)
         write(6,*)' Numerical Cartesian gradient:'
         write(6,'((3f20.10))')(cartgrd(i),i=1,nsize)
      END IF

      RETURN
      END

