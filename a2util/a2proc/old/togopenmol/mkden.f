      SUBROUTINE MKDEN(COEF,DENS,NBFP,
     $     NOCC,IFLG,ICORR,DENRELN,DENSCF,COEF2,TEMPMAT,IUHF,
     $     NDFT)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION COEF(NBFP*NBFP),DENS(NBFP*NBFP),
     $     NOCC(16),DENRELN(NBFP*NBFP),DENSCF(NBFP*NBFP),
     $     COEF2(NBFP*NBFP),TEMPMAT(NBFP*NBFP)
C     
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /IPAR/ LUOUT
      COMMON /ISYMINF/ NIRREP,NSOIRP(8)
C     
      NBAS2=NBFP*NBFP
C     
      FACTOR=0.5
      IF (IUHF.NE.0) FACTOR=1.0
C     
      CALL ZERO(DENSCF,NBAS2)
C
      IF(IFLG.EQ.0) THEN
C
C Get the matrix of alpha MO's
C
C Get the scf mo in the scf order
         IF (NDFT.EQ.0) THEN
            CALL GETREC(20,'JOBARC','SCFEVCAS',IINTFP*NBAS2,COEF)
         ELSE
            CALL GETREC(20,'JOBARC','SCFEVCA0',IINTFP*NBAS2,COEF)
         END IF
         ISOFF=0
         WRITE(*,*) " scf mo to so transformation"
         CALL PRINTMAT(COEF,NBFP,NBFP)
C
      ELSE
C Get the matrix of beta MO's
         IF (NDFT.EQ.0) THEN
            if (iuhf.eq.0) then
            else
               CALL GETREC(20,'JOBARC','SCFEVCBS',IINTFP*NBAS2,COEF)
            end if
         ELSE
            CALL GETREC(20,'JOBARC','SCFEVCB0',IINTFP*NBAS2,COEF)
         END IF
         ISOFF=8
         WRITE(*,*) " scf mo to so transformation"
         CALL PRINTMAT(COEF,NBFP,NBFP)
C     
         ENDIF
C     
         IF (ICORR.EQ.0) THEN
C     
C     Calculate the SCF density
            IOFF=0
            IOFF2=0
            DO 10 ISYM=1,NIRREP
               NOCCT=NOCC(ISOFF+ISYM)
               DO 20 IOCC=1,NOCCT
                  DO 30 I=1,NSOIRP(ISYM)
                     DO 40 J=1,NSOIRP(ISYM)
C     scf density
                        DENSCF((IOFF+I-1)*NBFP+IOFF+J)=
     $                       DENSCF((IOFF+I-1)*NBFP+IOFF+J)+
     $                       COEF((IOFF+IOCC-1)*NBFP+IOFF+I)
     $                       *COEF((IOFF+IOCC-1)*NBFP+IOFF+J)
C     
 40                  CONTINUE
 30               CONTINUE
 20            CONTINUE
               IOFF=IOFF+NSOIRP(ISYM)
               IOFF2=IOFF2+NOCCT*NBAS2
 10         CONTINUE
C     
            DO 50 I=1,NBAS2
               DENS(I)=DENSCF(I)
 50         CONTINUE
         ELSE
C     
            IF (IFLG.EQ.0) THEN
               CALL GETREC(20,'JOBARC','RELDENSA',IINTFP*NBAS2,DENRELN)
               WRITE(*,*) "correlated density matrix in mo basis"
               CALL PRINTMAT(DENRELN,NBFP,NBFP)
               CALL GETREC(20,'JOBARC','SCFEVCA0',IINTFP*NBAS2,COEF2)
               WRITE(*,*) "mo->ao trasformation"
               CALL PRINTMAT(COEF2,NBFP,NBFP)
            ELSE
               CALL GETREC(20,'JOBARC','RELDENSB',IINTFP*NBAS2,DENRELN)
               WRITE(*,*) "correlated density matrix in mo basis"
               CALL PRINTMAT(DENRELN,NBFP,NBFP)
               CALL GETREC(20,'JOBARC','SCFEVCB0',IINTFP*NBAS2,COEF2)
               WRITE(*,*) "mo->ao trasformation"
               CALL PRINTMAT(COEF2,NBFP,NBFP)
            END IF
C     
C     Transform the correlated density from the scf-mo basis to the scf-so basis
C     
            CALL XGEMM('N','N',NBFP,NBFP,NBFP,1.0D+00,COEF2,
     $           NBFP,DENRELN,NBFP,0.0D+00,TEMPMAT,NBFP)
            CALL XGEMM('N','T',NBFP,NBFP,NBFP,10.D+00,TEMPMAT,
     $           NBFP,COEF2,NBFP,0.0D+00,DENS,NBFP)
C     
            DO 60 I=1,NBAS2
               DENS(I)=FACTOR*DENS(I)
               DENRELN(I)=FACTOR*DENRELN(I)
 60         CONTINUE
C     
         WRITE(*,*) "SO based correlated density"
         CALL PRINTMAT(DENS,NBFP,NBFP)
C     WRITE(*,*) "DENRELN"
C     CALL PRINTMAT(DENRELN,NBFP,NBFP)
      END IF
C     
      RETURN
      END
      
