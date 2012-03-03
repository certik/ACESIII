      SUBROUTINE A2GET_DEN(COEF, DENS, SCFDEN, RELDEN, SPH2CART, 
     &                     TMP, NBAS, NBASP, NOCC, NSOIRP, IFLG)
     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SPHERICAL
      CHARACTER*8 LABELPSCF, LABELDSCF, LABELDENS, LABELSDEN
      CHARACTER*1 SPIN(2)
C
      DIMENSION COEF(NBASP*NBASP), DENS(NBASP*NBASP),SCFDEN(NBAS*NBAS),
     &          NOCC(16), RELDEN(NBAS*NBAS), FOCKMAT(NBASP, NBASP),
     &          SPH2CART(NBAS*NBASP), TMP(NBAS*NBAS), NSOIRP(8)
C
#include "machsp.com"
#include "syminf.com"
#include "jodaflags.com"
C
      DATA SPIN /'A', 'B'/
C
      SPHERICAL = (IFLAGS(62) .EQ. 1)
      NBAS2=NBASP*NBASP
C
      LABELPSCF='HFDENSTY'
      LABELDSCF='HDDENSTY'
      LABELDENS='TDENSITY'
      LABELSDEN='DDENSITY'
C
      CALL ZERO(DENS,   NBASP*NBASP)
      CALL ZERO(SCFDEN, NBAS*NBAS)
      CALL ZERO(RELDEN, NBAS*NBAS)
C
      IF (IFLG.EQ.1) THEN
C
         CALL GETREC(20,'JOBARC','SCFEVCA0',IINTFP*NBAS2,COEF)
         ISOFF=0
      ELSE
         CALL GETREC(20,'JOBARC','SCFEVCB0',IINTFP*NBAS2,COEF)
         ISOFF=8
      ENDIF
C
      IOFF  = 0
      IOFF2 = 0
      DO 10 ISYM=1,NIRREP
         NOCCT=NOCC(ISOFF+ISYM)
         DO 20 IOCC=1,NOCCT
            DO 30 I=1,NSOIRP(ISYM)
               DO 40 J=1,NSOIRP(ISYM)
                  DENS((IOFF+I-1)*NBASP+IOFF+J)=DENS((IOFF+I-1)*
     &                                          NBASP+IOFF+J)+
     &                                          COEF((IOFF+IOCC-1)*
     &                                          NBASP+IOFF+I)*
     &                                          COEF((IOFF+IOCC-1)*
     &                                          NBASP+IOFF+J)
C
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
C
         IOFF  = IOFF+NSOIRP(ISYM)
         IOFF2 = IOFF2+NOCCT*NBAS2
C
   10 CONTINUE
  
     IF (ZEROTH_ORDER_DENS) THEN
C
C Get the SCF and correlated (relaxed) density matirices in AO basis. For 
C open shell calculations, we also read the density diferences. Note that
C the "relaxed" density read here only contain the correlated contributions.
C
      CALL GETREC(-20, 'JOBARC', LABELPSCF, NBAS2*IINTFP, SCFDEN)
      CALL GETREC(-20, 'JOBARC', LABELDENS, NBAS2*IINTFP, RELDEN)
C      
      IF (IFLG .EQ. 2) THEN
         CALL GETREC(-20, 'JOBARC', LABELDSCF, NBAS2*IINTFP, SCFDEN)
         CALL GETREC(-20, 'JOBARC', LABELSDEN, NBAS2*IINTFP, RELDEN)
      ENDIF
C      
#ifdef _DEBUG_LVL0
      Write(6,*) "The Density matrix in AO basis"
      CALL OUTPUT(RELDEN, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
      ELSE IF (FIRST_ORDER_DENS) THEN
          
           DO IPERT = 1, NATOMS
C
              WRITE(RECNAMEA, "(A,I2)") "PTDENA", IPERT
              CALL GETREC(20, "JOBARC", RECNAMEA, NBAS2*IINTFP,
C
              WRITE(RECNAMEB, "(A,I2)") "PTDENB", IPERT
              CALL GETREC(20, "JOBARC", RECNAMEB, NBAS2*IINTFP,
     &                    )
              CALL SAXPY(NBAS2, 


C
C At this point we need to transform the density matrices to Cartesian
C basis in order to directly contract with the Cartesian AO integrals.
C 
C$$$      IF (SPHERICAL) THEN
C$$$         CALL GETREC(20, 'JOBARC','CMP2CART', NBAS*NBASP*IINTFP, 
C$$$     &               SPH2CART)

C$$$      Write(6,*) "The spherical to Cartesian transformation matrix"
C$$$      CALL OUTPUT(SPH2CART, 1, NBAS, 1, NBASP, NBAS, NBASP, 1)
C
C$$$         CALL TRANS_SPH2CART(SCFDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$         CALL TRANS_SPH2CART(RELDEN, SPH2CART, TMP, NBAS, NBASP)
C$$$      ENDIF
C
C
      RETURN
      END


