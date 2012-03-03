      SUBROUTINE HEADER(NATMS,NUC2,COORD2,NPTSX,TVALX,BVALX,
     $     NPTSY,TVALY,BVALY,NPTSZ,TVALZ,BVALZ,
     $     NBASPEFF,IDORO,IAORBD,ICORR,INORBS,IORBVECA,IORBVECB,
     $     NORBOUTA,NORBOUTB,IAORBO,NBASP,IORBOUT,IDENOUT)
C     
C     Print out the headers of the two output files.
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NUC2(NATMS),COORD2(NATMS*3),IORBVECA(9),
     $     IORBVECB(9)
C     
      STEPX=(TVALX-BVALX)/REAL(NPTSX)
      STEPY=(TVALY-BVALY)/REAL(NPTSY)
      STEPZ=(TVALZ-BVALZ)/REAL(NPTSZ)
      MODVAL=MOD((NBASPEFF+1),10)
      NGRPS=(NBASPEFF+1)/10
C     
C     Density Output:
      IF (IDORO.NE.1) THEN
         WRITE(IDENOUT,*) "ACES II CUBE OUTPUT"
         IF (ICORR.EQ.0) THEN
            IF (IAORBD.EQ.0) WRITE(IDENOUT,*) "SCF Alpha Density"
            IF (IAORBD.EQ.1) WRITE(IDENOUT,*) "SCF Beta Density"
            IF (IAORBD.EQ.2) WRITE(IDENOUT,*) "SCF Total Density"
         ELSE
            IF (IAORBD.EQ.0) WRITE(IDENOUT,*) "Correlated Alpha Density"
            IF (IAORBD.EQ.1) WRITE(IDENOUT,*) "Correlated Beta Density"
            IF (IAORBD.EQ.2) WRITE(IDENOUT,*) "Correlated Total Density"
         END IF
         WRITE(IDENOUT,"(I5,3 F12.6)") NATMS,BVALX,BVALY,BVALZ
         WRITE(IDENOUT,"(I5,3 F12.6)") (NPTSX+1),STEPX,0.0,0.0
         WRITE(IDENOUT,"(I5,3 F12.6)") (NPTSY+1),0.0,STEPY,0.0
         WRITE(IDENOUT,"(I5,3 F12.6)") (NPTSZ+1),0.0,0.0,STEPZ
         DO 20 I=1,NATMS
            WRITE(IDENOUT,"(I5,4 F12.6)") NUC2(I),REAL(NUC2(I)),
     $           (COORD2(3*(I-1)+J),J=1,3)
 20      CONTINUE
      END IF
C
C     Orbital Output
      IF (IDORO.NE.0) THEN
         WRITE(IORBOUT,*) "ACES II CUBE OUTPUT"
         IF (INORBS.EQ.0) THEN
            WRITE(IORBOUT,*) "SCF Molecular Orbitals"
         ELSE
            WRITE(IORBOUT,*) "Correlated Natural Orbitals"
         END IF
         WRITE(IORBOUT,"(I5,3 F12.6)") - NATMS,BVALX,BVALY,BVALZ
         WRITE(IORBOUT,"(I5,3 F12.6)") (NPTSX+1),STEPX,0.0,0.0
         WRITE(IORBOUT,"(I5,3 F12.6)") (NPTSY+1),0.0,STEPY,0.0
         WRITE(IORBOUT,"(I5,3 F12.6)") (NPTSZ+1),0.0,0.0,STEPZ
C     
         DO 30 I=1,NATMS
            WRITE(IORBOUT,"(I5,4 F12.6)") NUC2(I),REAL(NUC2(I)),
     $           (COORD2(3*(I-1)+J),J=1,3)
 30      CONTINUE
C
         IF (IAORBO.EQ.3) THEN
            WRITE(IORBOUT,"(10 I5)") 
     $           NBASPEFF,(IORBVECA(J),J=1,NORBOUTA),
     $           (IORBVECB(J)+NBASP,J=1,NORBOUTB)
         ELSE
C     
            IF (NGRPS.NE.0) THEN
               WRITE(IORBOUT,"(10 I5)") NBASPEFF,(J,J=1,9)
            ELSE
               WRITE(IORBOUT,"(10 I5)") NBASPEFF,(J,J=1,MODVAL-1)
            END IF
C     
            IF (NGRPS.GT.1) THEN
               DO 40 I=1,NGRPS-1
                  WRITE(IORBOUT, "(10 I5)") (9+( (I-1)*10)+J,J=1,10)
 40            CONTINUE
            END IF
C     
            IF ((MODVAL.NE.0).AND.(NGRPS.NE.0)) THEN
               WRITE(IORBOUT,"(10 I5)") (9+((NGRPS-1)*10)+
     $              J,J=1,MODVAL)
            END IF
C     
         END IF
      END IF
C     
      RETURN
      END SUBROUTINE HEADER

