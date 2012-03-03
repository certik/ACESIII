
C THIS ROUTINE BUILDS A CHARACTER TABLE FOR ANY POINT GROUP
C AND WRITES IT TO DISK ALONG WITH LABELS FOR EACH IRREDUCIBLE
C REPRESENTATION.  THE ALGORITHM IS AS FOLLOWS:
C
C   1. A SET OF POINTS IS GENERATED WHICH CONSTITUTES THE LARGEST
C      POSSIBLE ORBIT (POPULATION h) OF THE POINT GROUP.  SUCH A
C      SET CONTAINS THE REGULAR REPRESENTATION, MEANING THAT EACH
C      IRREP IS INCLUDED EXACTLY ONE TIME.
C
C   2. THE DISTANCE MATRIX FOR THESE POINTS IS CONSTRUCTED AND
C      THEN DIAGONALIZED.  THE EIGENVECTORS OF THIS MATRIX TRANSFORM
C      AS THE IRREPS.
C
C   3. EACH EIGENVECTOR IS THEN NORMALIZED, AND CHARACTERS ARE COMPUTED
C      FOR ALL OPERATIONS OF THE GROUP.

c INPUT
c integer IORDER
c char*4  TYPE

c OUTPUT
c double  SYOPS(9*IORDER)
c double  REGREP(9*IORDER*IORDER)
c double  DIST(IORDER,IORDER)
c double  CHAR(IORDER,IORDER)
c integer ISCR(3*IORDER)
c integer IPTR(IORDER,IORDER)
c integer IDEGEN(IORDER)
c char*8  IRRNM(IORDER)
c double  SCR(3)
c double  DSCR(3*IORDER)

c RECORDS
c get TYPE//'PTGP'
c get TYPE//'SYOP'
c put TYPE//'CHAR'
c put TYPE//'DEGN'
c put TYPE//'LABL'
c put TYPE//'NIRX'

      SUBROUTINE CHRTABLE(IORDER,SYOPS,REGREP,DIST,CHAR,ISCR,IPTR,
     &                    IDEGEN,IRRNM,SCR,TYPE,DSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SYOPS(9*IORDER),REGREP(9*IORDER*IORDER)
      DIMENSION DIST(IORDER,IORDER),CHAR(IORDER,IORDER),ISCR(3*IORDER)
      DIMENSION IPTR(IORDER,IORDER),IDEGEN(IORDER)
      CHARACTER*8 IRRNM(IORDER)
      DIMENSION SCR(3)
      CHARACTER*4 TYPE
      DIMENSION DSCR(3*IORDER)

      CHARACTER*4 PTGRP

      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWd

      DATA TOL  /1.D-8/
      DATA TOL2 /1.D-5/

      CALL GETCREC(20,'JOBARC',TYPE//'PTGP',4,PTGRP)

c   o pick an arbitrary point
      SCR(1)=SQRT(3.432)
      SCR(2)=SQRT(4.920)
      SCR(3)=SQRT(13.25626)

c   o read all symmetry operations
      CALL GETREC(20,'JOBARC',TYPE//'SYOP',IORDER*9*IINTFP,SYOPS)

c   o apply all operations of the group to the reference point
c     (operation N maps the reference atom into atom N)
      IOFFOPS=1
      IOFFREP=1
      DO I=1,IORDER
         CALL XGEMM('N','N',3,1,3,1.d0,SYOPS(IOFFOPS),3,SCR,3,
     &              0.d0,REGREP(IOFFREP),3)
         IOFFREP=IOFFREP+3
         IOFFOPS=IOFFOPS+9
      END DO

c   o sort coordinates to facilitate construction of permutation vector
      DO I=1,IORDER
         ISCR(I)=I
      END DO
      CALL SORT2(REGREP,SCR,ISCR,IORDER)

C NOW COMPUTE PERMUTATION VECTOR WHICH RELATES ALL POSITIONS
C TO THEIR IMAGE UNDER THE TRANSFORMATION

      IOFFOPS=1
      IOFFS=3*IORDER+1
      DO I=1,IORDER
         CALL XGEMM('N','N',3,IORDER,3,1.d0,SYOPS(IOFFOPS),3,REGREP,3,
     &              0.d0,SCR,3)
         CALL SORT2(SCR,SCR(IOFFS),ISCR(IORDER+1),IORDER)
         CALL STPTR(IORDER,ISCR,ISCR(IORDER+1),IPTR(1,I))
         IOFFOPS=IOFFOPS+9
      END DO

c   o construct distance matrix for points in regrep
      CALL ZERO(DIST,IORDER*IORDER)
      DO I=2,IORDER
         DO J=1,I-1
            IOFFI=3*(I-1)+1
            IOFFJ=3*(J-1)+1
            CALL VADD(SCR,REGREP(IOFFI),REGREP(IOFFJ),3,-1.d0)
            X=DNRM2(3,SCR,1)
            DIST(I,J)=X
            DIST(J,I)=X
         END DO
      END DO

c   o diagonalize distance matrix
      CALL DSYEV('V','L',IORDER,DIST,IORDER,SCR,DSCR,3*IORDER,INFO)
      INFO = 0
      IF (INFO.NE.0) THEN
         WRITE(*,*)
     &        '@CHRTABLE: There was a problem diagonalizing the matrix.'
         CALL ERREX
      END IF
      DO I = 1, IORDER-1
         DO J = I+1, IORDER
            IF ( SCR(I) .LT. SCR(J) ) then
               DTMP   = SCR(I)
               SCR(I) = SCR(J)
               SCR(J) = DTMP
               CALL DSWAP(IORDER,DIST(1,I),1,DIST(1,J),1)
            END IF
         END DO
      END DO

c   o collect sets of nondegenerate eigenvectors
      XTEST=-9999999.
      IOFF=1
      NIRREP=0
      IDEG=0
      DO I=1,IORDER
         Z=ABS(SCR(I+1)-XTEST)
         XTEST=SCR(I+1)
         IF (Z.GT.TOL.OR.I.EQ.IORDER) THEN
            ILOC1=I-IDEG
            LENGTH=(IDEG+1)*IORDER
            NIRREP=NIRREP+1
C NOW PERMUTE THE ATOMS ACCORDING TO THE PERMUTATION IMPLIED BY THE
C OPERATION AND COMPUTE THE OVERLAP.  THIS IS THE CHARACTER.
            DO IOP=1,IORDER
               IOFF=IORDER
               IPOS=ILOC1
               DO IVEC=1,IDEG+1
                  DO IPOSIN=1,IORDER
                     IPOSOUT=IPTR(IPOSIN,IOP)
                     SCR(IPOSOUT+IOFF)=DIST(IPOSIN,IPOS)
                  END DO
                  IOFF=IOFF+IORDER
                  IPOS=IPOS+1
               END DO
               CHAR(IOP,NIRREP)=DDOT(LENGTH,DIST(1,ILOC1),1,
     &                                      SCR(IORDER+1),1)
               IDEGEN(NIRREP)=IDEG+1
            END DO
            IDEG=0
         ELSE
            IDEG=IDEG+1
         END IF
      END DO

c   o squeeze out redundant irreps
      IRRUNQ=1
      DO I=2,NIRREP
         X=0.
         DO J=1,IRRUNQ
            X=X+DDOT(IORDER,CHAR(1,I),1,CHAR(1,J),1)
         END DO
         IF (ABS(X).LT.TOL2) THEN
            IRRUNQ=IRRUNQ+1
            CALL DCOPY(IORDER,CHAR(1,I),1,CHAR(1,IRRUNQ),1)
            IDEGEN(IRRUNQ)=IDEGEN(I)
         END IF
      END DO

      WRITE(6,1000)IRRUNQ
1000  FORMAT(T3,'@CHRTABLE-I, There are ',I3,' unique irreducible ',
     &          'representations.')
      IF (IFLAGS(1).GE.1) THEN
         WRITE(6,2000)
2000     FORMAT(T3,'             Expanded character table follows ')
      END IF
      DO I=1,IRRUNQ
         CALL ASSIRR(PTGRP,CHAR(1,I),IDEGEN(I),IRRNM(I)(1:4))
         IF (IFLAGS(1).GE.1) THEN
            WRITE(6,110)I,IRRNM(I)(1:4),IDEGEN(I)
110         FORMAT(T3,' Irrep : ',I3,/,
     &             T3,' Label : ',A ,/,
     &             T3,' Degen : ',I1)
            WRITE(6,'((12F6.3))')(CHAR(J,I),J=1,IORDER)
         END IF
      END DO

      CALL PUTREC(20,'JOBARC',TYPE//'CHAR',IRRUNQ*IORDER*IINTFP,CHAR)
      CALL PUTREC(20,'JOBARC',TYPE//'DEGN',IRRUNQ,IDEGEN)
      CALL PUTREC(20,'JOBARC',TYPE//'LABL',IRRUNQ*IINTFP,IRRNM)
      CALL PUTREC(20,'JOBARC',TYPE//'NIRX',1,IRRUNQ)

      RETURN
      END

