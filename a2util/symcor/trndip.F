
C THIS ROUTINE TRANSFORMS THE SYMMETRY COORDINATE DIPOLE DERIVATIVE
C MATRIX TO THE CARTESIAN REPRESENTATION

c INPUT
c integer NMODE      :
c integer NATOM      :
c integer NIRREP     :
c integer ISYMIRR(*) :
c integer IDEGEN(*)  :
c char*4  TYPE       :
c integer NDSCR      :

c OUTPUT
c double DIPDFUL(3,3*NATOM) :
c double POLDFUL(9,3*NATOM) :
c double DSCR(NDSCR)        :

c INPUT/OUTPUT
c double DIPDER(3,3*NATOM) :
c double POLDER(9,3*NATOM) :

c RECORDS
c get TYPE//'SYMQ'
c get 'ATOMMASS'

      SUBROUTINE TRNDIP(NMODE,NATOM,NIRREP,
     &                  DIPDER,DIPDFUL,POLDER,POLDFUL,
     &                  ISYMIRR,IDEGEN,TYPE,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION DIPDER(3,3*NATOM),DIPDFUL(3,3*NATOM)
      DIMENSION POLDER(9,3*NATOM),POLDFUL(9,3*NATOM)
      DIMENSION ISYMIRR(*),IDEGEN(*)
      CHARACTER*4 TYPE
      double precision dscr(ndscr)

      LOGICAL bExist

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "control.com" /* for RAMAN */

      DATA ZILCH /0.0/

      if (eneronly) then
c      o there are no derivatives
         CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,DSCR)
         OPEN(UNIT=40,FILE='DIPDER',STATUS='UNKNOWN',FORM='FORMATTED')
         DO IXYZ=1,3
            WRITE(40,'(I5)')IXYZ
            DO IATOM=1,NATOM
               IF (DSCR(IATOM).NE.ZILCH) THEN
                  WRITE(40,'((4F20.10))')0.d0,0.d0,0.d0,0.d0,0.d0
               END IF
            END DO
         END DO
         CLOSE(UNIT=40,STATUS='KEEP')
         return
      end if

      NATOM3=3*NATOM
      NATOM9=9*NATOM

      IF (NDSCR.LT.(NATOM3+9)*NATOM3) THEN
         print *, '@TRNDIP: Insufficient memory.'
         print *, '         have ',NDSCR,' doubles'
         print *, '         need ',(NATOM3+9)*NATOM3,' doubles'
         call aces_exit(1)
      END IF

      CALL ZERO(DIPDFUL,3*NATOM3)
      CALL ZERO(POLDFUL,9*NATOM3)

      ICOUNT=1
      IOFF1 =1
      NLEFT =NMODE
      DO IRREP=1,NIRREP

c      o find first occurance of this irrep
         ILOC=ISRCHEQ(NMODE,ISYMIRR,1,IRREP)
         IF (ILOC.NE.NMODE+1) THEN
            ILAST=ISRCHNE(NLEFT,ISYMIRR(ILOC),1,IRREP)
            NVIBSYM=ILAST-1
            IDEG=IDEGEN(IRREP)
            NDIM=NVIBSYM/IDEG
c         o expand dipder and polder matrix into full matrices
            CALL DSCAL(3*NDIM,DFLOAT(IDEG),DIPDER(1,ioff1),1)
            CALL DCOPY(3*NDIM,DIPDER(1,IOFF1),1,DIPDFUL(1,ICOUNT),1)
            CALL DSCAL(9*NDIM,DFLOAT(IDEG),POLDER(1,ioff1),1)
            CALL DCOPY(9*NDIM,POLDER(1,IOFF1),1,POLDFUL(1,ICOUNT),1)
            ICOUNT=ICOUNT+NVIBSYM
            IOFF1 =IOFF1 +NDIM
            NLEFT =NLEFT -NVIBSYM
         END IF

      END DO

c   o transform to mass-weighted cartesian coordinates
      I = 1+3*NATOM9
      CALL GETREC(20,'JOBARC',TYPE//'SYMQ',NATOM3*NATOM3*IINTFP,DSCR(I))
      CALL XGEMM('N','T',3,NATOM3,NATOM3,
     &           1.d0,DIPDFUL,3,
     &                DSCR(I),NATOM3,
     &           0.d0,DSCR,   3)
      CALL DCOPY(NATOM9,DSCR,1,DIPDFUL,1)
      CALL XGEMM('N','T',9,NATOM3,NATOM3,
     &           1.d0,POLDFUL,9,
     &                DSCR(I),NATOM3,
     &           0.d0,DSCR,   9)
      CALL DCOPY(3*NATOM9,DSCR,1,POLDFUL,1)

c   o remove mass weighting
      CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,DSCR)
      IOFF1=1
      DO IATOM=1,NATOM
         X=SQRT(DSCR(IATOM))
         DO IXYZ=1,3
            DIPDFUL(1,IOFF1)=X*DIPDFUL(1,IOFF1)
            DIPDFUL(2,IOFF1)=X*DIPDFUL(2,IOFF1)
            DIPDFUL(3,IOFF1)=X*DIPDFUL(3,IOFF1)
            DO I = 1, 9
               POLDFUL(I,IOFF1)=X*POLDFUL(I,IOFF1)
            END DO
            IOFF1=IOFF1+1
         END DO
      END DO

      OPEN(UNIT=40,FILE='DIPDER',STATUS='UNKNOWN',FORM='FORMATTED')
      DO IXYZ=1,3
         WRITE(40,'(I5)')IXYZ
         IOFF=1
         DO IATOM=1,NATOM
            IF (DSCR(IATOM).NE.ZILCH) THEN
               WRITE(40,'((4F20.10))')0.d0,
     &                                (DIPDFUL(IXYZ,I),I=IOFF,IOFF+2)
            END IF
            IOFF=IOFF+3
         END DO
      END DO
      CLOSE(UNIT=40,STATUS='KEEP')

      INQUIRE(FILE='POLDER', EXIST=bExist)
      IF (.NOT.bExist.AND.RAMAN) THEN
         OPEN(UNIT=60,FILE='POLDER',STATUS='UNKNOWN',FORM='FORMATTED')
         DO IXYZ=1,9
            WRITE(60,'(I5)')IXYZ
            IOFF=1
            DO IATOM=1,NATOM
               IF (DSCR(IATOM).NE.ZILCH) THEN
                  WRITE(60,'((4F20.10))')0.d0,
     &                                   (POLDFUL(IXYZ,I),I=IOFF,IOFF+2)
               END IF
               IOFF=IOFF+3
            END DO
         END DO
         CLOSE(UNIT=60,STATUS='KEEP')
      END IF

      RETURN
      END

