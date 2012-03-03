      SUBROUTINE ORBOUT(NBASP,NORBOUTA,NORBOUTB,PCOEFFA,
     $     PCOEFFB,IORBVECA,IORBVECB,ORDERMATA,ORDERMATB,
     $     ITFCT,NORBOUT,IORBVEC)
c     
C     This subroutine redefines the pcoeffa and pcoeffb
C     when only a few orbitals are going to be displayed.
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      DIMENSION PCOEFFA(ITFCT*NBASP),PCOEFFB(ITFCT*NBASP),
     $     IORBVECA(9),IORBVECB(9),ORDERMATA(9*ITFCT),
     $     ORDERMATB(9*ITFCT),IORBVEC(9)
C     
C     Divide the iorbvec into its alpha and beta components.
      IA=0
      IB=0
      DO 5 I=1,9
         IORBVECA(I)=0
         IORBVECB(I)=0
         IF (IORBVEC(I).GT.NBASP) THEN
            IORBVECB(IB+1)=IORBVEC(I)-NBASP
            IB=IB+1
         ELSEIF(IORBVEC(I).NE.O) THEN
            IORBVECA(IA+1)=IORBVEC(I)
            IA=IA+1
         END IF
 5    CONTINUE
C     
C     Fill in the temporary matrices with the choice orbitals
      DO 10 I=1,9
         IF (IORBVECA(I).NE.O) THEN
            DO 20 J=1,ITFCT
               ORDERMATA((I-1)*ITFCT+J)=PCOEFFA((IORBVECA(I)-1)*ITFCT+J)
 20         CONTINUE 
         END IF
C     
         IF (IORBVECB(I).NE.O) THEN
            DO 30 J=1,ITFCT
               ORDERMATB((I-1)*ITFCT+J)=PCOEFFB((IORBVECB(I)-1)*ITFCT+J)
 30         CONTINUE
         END IF
 10   CONTINUE
C     
      NORBOUTA=IA
      NORBOUTB=IB
      NORBOUT=NORBOUTA+NORBOUTB
C     
C     Now refill the hereafter useful part of pcoeffa and pcoeffb.
C     
      DO 40 I=1,ITFCT*NORBOUTA
         PCOEFFA(I)=ORDERMATA(I)
 40   CONTINUE
C     
      DO 50 I=1,ITFCT*NOR80UT8
         PCOEFFB(I)=ORDERMATB(I)
 50   CONTINUE
C     
      RETURN
      END SUBROUTINE ORBOUT
      
