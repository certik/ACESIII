
C THIS PROJECTS THE ROTATIONAL MODES FROM THE SYMMETRY
C ADAPTED COORDINATES

c INPUT
c integer NATOM

c OUTPUT
c double TPROJ(*)
c double SCR(*)
c double SYMQ(*)
c double ATMASS(*)

c RECORDS
c get 'LINEAR  '
c get 'ATOMMASS'
c get ROTREC(IXYZ)
c put ROTREC(IXYZ)
c get 'COMPNSYQ'
c get 'COMPSYMQ'
c put 'COMPSYMQ'

      SUBROUTINE ROTPRJ(NATOM,TPROJ,SCR,SYMQ,ATMASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION SCR(*),SYMQ(*),TPROJ(*),ATMASS(*)

      CHARACTER*8 ROTREC(3)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA ONE /1.0/
      DATA TOL /1.D-8/
      DATA ROTREC /'ROTVECX ','ROTVECY ','ROTVECZ '/

      CALL GETREC(20,'JOBARC','LINEAR  ',1,ILINEAR)

C CONSTRUCT THE ROTATIONAL PROJECTOR = 1 - SUM |Rq><Rq|
C                                          xyz

      NSIZE=3*NATOM
      CALL ZERO(TPROJ,NSIZE*NSIZE)
      CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,ATMASS)
      IOFF=1
      DO I=1,NSIZE
         TPROJ(IOFF)=ONE
         IOFF=IOFF+NSIZE+1
      END DO
      DO IXYZ=1,3-ILINEAR
c      o mass-weight and write back to disk
         CALL GETREC(20,'JOBARC',ROTREC(IXYZ),NSIZE*IINTFP,SCR)
         IOFF=1
         DO I=1,NATOM
            Z=SQRT(ATMASS(I))
            CALL DSCAL(3,Z,SCR(IOFF),1)
            IOFF=IOFF+3
         END DO
         Z=DNRM2(NSIZE,SCR,1)
         FACT=ONE/Z
         CALL DSCAL(NSIZE,FACT,SCR,1)
         CALL PUTREC(20,'JOBARC',ROTREC(IXYZ),NSIZE*IINTFP,SCR)
         CALL XGEMM('N','N',NSIZE,NSIZE,1,
     &              -1.d0,SCR,  NSIZE,
     &                    SCR,  1,
     &              ONE,  TPROJ,NSIZE)
      END DO

c   o hit symmetry adapted coordinates with projector
      CALL GETREC(20,'JOBARC','COMPNSYQ',1,NSYMOLD)
      CALL GETREC(20,'JOBARC','COMPSYMQ',NSYMOLD*NSIZE*IINTFP,SCR)
      CALL XGEMM('N','N',NSIZE,NSYMOLD,NSIZE,ONE,TPROJ,NSIZE,
     &           SCR,NSIZE,0.d0,SYMQ,NSIZE)

c   o renormalize
      IOFF=1
      DO I=1,NSYMOLD
         X=DNRM2(NSIZE,SYMQ(IOFF),1)
         IF (ABS(X).GT.TOL) THEN
            Z=ONE/X
            CALL DSCAL(NSIZE,Z,SYMQ(IOFF),1)
         END IF
         IOFF=IOFF+NSIZE
      END DO

      CALL PUTREC(20,'JOBARC','COMPSYMQ',NSIZE*NSYMOLD*IINTFP,SYMQ)

      RETURN
      END

