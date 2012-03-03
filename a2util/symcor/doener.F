
C THIS ROUTINE DETERMINES WHICH ENERGIES ARE CALCULATED TO
C EVALUATE THE FORCE CONSTANT MATRIX
C AND GENERATES THE APPROPRIATE *CARTESIAN* GEOMETRIES

c INPUT
c integer NATOM              : NUMBER OF ATOMS IN THE MOLECULE (WITH DUMMIES)
c integer NDIM               : NUMBER OF COORDINATES IN THIS IRREP
c double  SYMQ(3*NATOM,NDIM) : SYMMETRY ADAPTED COORDINATES FOR THIS IRREP
c double  COORD(3*NATOM)     : CARTESIAN COORDINATES OF THE REFERENCE STRUCTURE
c double  STPSIZ             : STEP SIZE FOR DISPLACEMENTS IN MASS-WEIGHTED
c                              CARTESIAN COORDINATES
c double  VMASS(NATOM)       : INVERSE SQUARE ROOTS OF THE ATOMIC MASSES
c integer INVOP(NDIM)
c logical PRINTQ             : VERBOSE PRINTING FLAG
c integer NDSCR              : AMOUNT OF DOUBLE SCRATCH AT DSCR

c OUTPUT
c double  POINTS(3*NATOM,*) : COORDINATES USED IN FD ENERGY CALCULATIONS
c integer NPOINT            : NUMBER OF ENERGY CALCS REQ'D FOR THIS IRREP
c double  DSCR(NDSCR)       : (scr) double scratch

#include "flags.h"

      SUBROUTINE DOENER(NATOM,NDIM,SYMQ,COORD,STPSIZ,VMASS,
     &                  POINTS,NPOINT,
     &                  INVOP,PRINTQ,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SYMQ(3*NATOM,NDIM),COORD(3*NATOM),VMASS(NATOM)
      DIMENSION POINTS(*),INVOP(NDIM),DSCR(NDSCR)
      LOGICAL PRINTQ

      LOGICAL ONEGRD
      CHARACTER*5 PHASE(2)

#include "control.com" /* for GMTRYOPT */

      DATA PHASE /'plus ','minus'/

      NSIZE=3*NATOM
      PRINTQ = .TRUE. 

      if (ndscr.lt.nsize*2) then
         print *, '@DOENER: Insufficient memory.'
         print *, '         have ',ndscr,' doubles'
         print *, '         need ',nsize*2,' doubles'
         call aces_exit(1)
      end if

C LOOP OVER DIMENSIONALITY OF SUBSPACE.  FOR ENERGIES, WE NEED
C PURE AND MIXED DISPLACEMENTS FOR ALL POSSIBLE MODES.  THERE
C WILL BE [NDIM*(NDIM+1)]/2 OF THESE.

      IOFFP=1
      NPOINT=0
      DO IDIM1=1,NDIM

C SCALE SYMMETRY COORDINATE VECTOR FOR DIMENSION IDIM1 AND TRANSFORM
C TO PURE CARTESIAN COORDINATES

         CALL DCOPY(NSIZE,SYMQ(1,IDIM1),1,DSCR,1)
         CALL DSCAL(NSIZE,STPSIZ,DSCR,1)
         NDX = 1
         DO IATOM=1,NATOM
            DSCR(NDX+0) = DSCR(NDX+0)*VMASS(IATOM)
            DSCR(NDX+1) = DSCR(NDX+1)*VMASS(IATOM)
            DSCR(NDX+2) = DSCR(NDX+2)*VMASS(IATOM)
            NDX = NDX+3
         END DO

C LOOP OVER SECOND DISPLACEMENT AND DO THE SAME STUFF

         DO IDIM2=1,IDIM1
            ONEGRD=(INVOP(IDIM1).GT.0).AND.(INVOP(IDIM2).GT.0)
         if (.not.gmtryopt.or.idim1.eq.idim2) then
            CALL DCOPY(NSIZE,SYMQ(1,IDIM2),1,DSCR(NSIZE+1),1)
            CALL DSCAL(NSIZE,STPSIZ,DSCR(NSIZE+1),1)
            NDX = 1+NSIZE
            DO IATOM=1,NATOM
               DSCR(NDX+0) = DSCR(NDX+0)*VMASS(IATOM)
               DSCR(NDX+1) = DSCR(NDX+1)*VMASS(IATOM)
               DSCR(NDX+2) = DSCR(NDX+2)*VMASS(IATOM)
               NDX = NDX+3
            END DO

            IF (IDIM1.EQ.IDIM2) THEN

c            o generate positive displacement
               CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,1.d0)
               IOFFP=IOFFP+NSIZE
               NPOINT=NPOINT+1

c            o generate negative displacement (totally symmetric irrep only)
               IF (.NOT.ONEGRD) THEN
                  CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,-1.d0)
                  IOFFP=IOFFP+NSIZE
                  NPOINT=NPOINT+1
               END IF

            ELSE

c            o generate +/+ displacements
               CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,1.d0)
               CALL DAXPY(NSIZE,1.d0,DSCR(NSIZE+1),1,POINTS(IOFFP),1)
               IOFFP=IOFFP+NSIZE
               NPOINT=NPOINT+1

c            o generate -/- displacements
               IF (.NOT.ONEGRD) THEN
                  CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,-1.d0)
                  CALL DAXPY(NSIZE,-1.d0,DSCR(NSIZE+1),1,
     &                                   POINTS(IOFFP),1)
                  IOFFP=IOFFP+NSIZE
                  NPOINT=NPOINT+1
               END IF

            END IF

         end if
         END DO

      END DO

      IF (PRINTQ) THEN
         WRITE(6,1000)
1000     FORMAT(T3,'Geometries used for energy calculations : ')
         IOFF=1
         IF (ONEGRD) THEN
            ITOP=1
         ELSE
            ITOP=2
         END IF
         DO IDIM1=1,NDIM
            DO IDIM2=1,IDIM1
            if (.not.gmtryopt.or.idim1.eq.idim2) then
               DO IPHASE=1,ITOP
                  IF (IDIM1.EQ.IDIM2) THEN
                     WRITE(6,1001)IDIM1,PHASE(IPHASE)
1001                 FORMAT(T3,'Symmetry coordinate : ',i3,' Phase : ',
     &                      a,' Type : Energy')
                  ELSE
                     WRITE(6,1002)IDIM1,IDIM2,PHASE(IPHASE)
1002                 FORMAT(T3,'Symmetry coordinates (mixed) : ',i3,i3,
     &                      ' Phase : ',a,' Type : Energy')
                  END IF
                  DO IATOM=1,NATOM
                     WRITE(6,1003)IATOM,(POINTS(IPOS),IPOS=IOFF,IOFF+2)
1003                 FORMAT(T3,I5,3F20.10)
                     IOFF=IOFF+3
                  END DO
               END DO
            end if
            END DO
         END DO
      END IF

      RETURN
      END

