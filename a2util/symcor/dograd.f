
C THIS ROUTINE DETERMINES WHICH GRADIENTS ARE CALCULATED TO
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
c double  POINTS(3*NATOM,*) : COORDINATES USED IN FD GRADIENT CALCULATIONS
c integer NPOINT            : NUMBER OF GRADIENT CALCS REQ'D FOR THIS IRREP
c double  DSCR(NDSCR)       : (scr) double scratch

      SUBROUTINE DOGRAD(NATOM,NDIM,SYMQ,COORD,STPSIZ,VMASS,
     &                  POINTS,NPOINT,
     &                  INVOP,PRINTQ,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SYMQ(3*NATOM,NDIM),COORD(3*NATOM),VMASS(NATOM)
      DIMENSION POINTS(*),INVOP(NDIM),DSCR(NDSCR)
      LOGICAL PRINTQ

      LOGICAL ONEGRD
      CHARACTER*5 PHASE

      NSIZE=3*NATOM

      if (ndscr.lt.nsize) then
         print *, '@DOGRAD: Insufficient memory.'
         print *, '         have ',ndscr,' doubles'
         print *, '         need ',nsize,' doubles'
         call aces_exit(1)
      end if

      IOFFP=1
      NPOINT=0

c   o loop over dimensionality of subspace
      DO IDIM=1,NDIM
         ONEGRD=(INVOP(IDIM).GT.0)

c      o scale symmetry coordinate vector
         CALL DCOPY(NSIZE,SYMQ(1,IDIM),1,DSCR,1)
         CALL DSCAL(NSIZE,STPSIZ,DSCR,1)

c      o transform from mass-weighted cartesians to pure cartesians
         NDX = 1
         DO IATOM=1,NATOM
            DSCR(NDX+0) = DSCR(NDX+0)*VMASS(IATOM)
            DSCR(NDX+1) = DSCR(NDX+1)*VMASS(IATOM)
            DSCR(NDX+2) = DSCR(NDX+2)*VMASS(IATOM)
            NDX = NDX+3
         END DO

c      o generate positive displacement
         CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,1.d0)
         IOFFP=IOFFP+NSIZE
         NPOINT=NPOINT+1

c      o generate negative displacement (totally symmetric irrep only)
         IF (.NOT.ONEGRD) THEN
            CALL VADD(POINTS(IOFFP),COORD,DSCR,NSIZE,-1.d0)
            IOFFP=IOFFP+NSIZE
            NPOINT=NPOINT+1
         END IF

      END DO

      IF (PRINTQ) THEN
         IOFF=1
         DO IPOINT=1,NPOINT
            IF (.NOT.ONEGRD.AND.MOD(IPOINT,2).EQ.0) THEN
               PHASE='minus'
               ISYCOR=1+(IPOINT-1)/2
            ELSE IF (.NOT.ONEGRD.AND.MOD(IPOINT,2).EQ.1) THEN
               ISYCOR=(IPOINT+1)/2
               PHASE='plus '
            ELSE IF (ONEGRD) THEN
               ISYCOR=IPOINT
               PHASE='plus'
            END IF
            WRITE(6,1001)ISYCOR,PHASE
1001        FORMAT(T3,'Symmetry coordinate : ',i3,' Phase : ',a,
     &             ' Type : Gradient')
            DO IATOM=1,NATOM
               WRITE(6,1002)IATOM,(POINTS(IPOS),IPOS=IOFF,IOFF+2)
1002           FORMAT(T3,I5,3F20.10)
               IOFF=IOFF+3
            END DO
         END DO
      END IF

      RETURN
      END

