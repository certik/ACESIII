
C THIS ROUTINE TRANSFORMS THE SYMMETRY COORDINATE HESSIAN MATRIX
C TO THE CARTESIAN REPRESENTATION.

c INPUT
c integer NMODE      :
c integer NATOM      :
c integer NIRREP     :
c integer ISYMIRR(*) :
c integer IDEGEN(*)  :
c char*4  TYPE       :
c integer NDSCR      :

c INPUT/OUTPUT
c double FCM(9*NATOM*NATOM) :

c OUTPUT
c double FCMFUL(9*NATOM*NATOM) :
c double DSCR(NDSCR)           :

c RECORDS
c get TYPE//'SYMQ'
c get 'ATOMMASS'

      SUBROUTINE TRNFCM(NMODE,NATOM,NIRREP,FCM,FCMFUL,
     &                  ISYMIRR,IDEGEN,TYPE,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION FCM(9*NATOM*NATOM),FCMFUL(9*NATOM*NATOM)
      DIMENSION ISYMIRR(*),IDEGEN(*)
      CHARACTER*4 TYPE
      double precision dscr(ndscr)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      NSIZE=3*NATOM

      IF (NDSCR.LT.2*NSIZE*NSIZE) THEN
         print *, '@TRNFCM: Insufficient memory.'
         print *, '         have ',NDSCR,' doubles'
         print *, '         need ',2*NSIZE*NSIZE,' doubles'
         call aces_exit(1)
      END IF

      CALL ZERO(FCMFUL,NSIZE*NSIZE)

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
c         o expand fcm matrix into full matrix
            call dscal(ndim*ndim,dfloat(ideg),fcm(ioff1),1)
            CALL BLKCPY(FCM(IOFF1),NDIM,NDIM,FCMFUL,NSIZE,NSIZE,
     &                  ICOUNT,ICOUNT)
            ICOUNT=ICOUNT+NVIBSYM
            IOFF1 = IOFF1+NDIM*NDIM
            NLEFT = NLEFT-NVIBSYM
         END IF

      END DO

c   o transform to mass-weighted cartesian coordinates
      CALL GETREC(20,'JOBARC',TYPE//'SYMQ',NSIZE*NSIZE*IINTFP,DSCR)
      I = 1+NSIZE*NSIZE
      CALL XGEMM('N','N',NSIZE,NSIZE,NSIZE,
     &           1.d0,DSCR,   NSIZE,
     &                FCMFUL, NSIZE,
     &           0.d0,DSCR(I),NSIZE)
      CALL XGEMM('N','T',NSIZE,NSIZE,NSIZE,
     &           1.d0,DSCR(I),NSIZE,
     &                DSCR,   NSIZE,
     &           0.d0,FCMFUL, NSIZE)

c   o remove mass weighting
      CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,DSCR)
      IOFF1=1
      IOFF2=1
      DO IATOM=1,NATOM
         X=SQRT(DSCR(IATOM))
         DO IXYZ=1,3
            CALL DSCAL(NSIZE,X,FCMFUL(IOFF1),1)
            CALL DSCAL(NSIZE,X,FCMFUL(IOFF2),NSIZE)
            IOFF1=IOFF1+NSIZE
            IOFF2=IOFF2+1
         END DO
      END DO

c   o create the FCM file
      OPEN(UNIT=40,FILE='FCM',STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(40,'(2I5)')NSIZE,NSIZE
      WRITE(40,'((3F20.10))')(FCMFUL(I),I=1,NSIZE*NSIZE)
      CLOSE(UNIT=40,STATUS='KEEP')

      RETURN
      END

