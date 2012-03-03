      SUBROUTINE A2BULT_PRDUCT(INDX, JNDX, ITYPE, JTYPE, MAXPRM,  
     &                         EXP1, EXP2, CENTER, CNTMU, CNTNU, 
     &                         TMP2)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION TMP2(MAXPRM, MAXPRM), CENTER(3), CNTMU(3), CNTNU(3)
C
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
C

      LI = LMNVAL(1, ITYPE) 
      MI = LMNVAL(2, ITYPE) 
      NI = LMNVAL(3, ITYPE) 
C
      LJ = LMNVAL(1, JTYPE) 
      MJ = LMNVAL(2, JTYPE) 
      NJ = LMNVAL(3, JTYPE) 
C
      DISTNCEI = ((CENTER(1)-CNTMU(1))**2 + (CENTER(2)-CNTMU(2))**2 +
     &           (CENTER(3)-CNTMU(3))**2)
      DISTNCEJ = ((CENTER(1)-CNTNU(1))**2 + (CENTER(2)-CNTNU(2))**2 +
     &           (CENTER(3)-CNTNU(3))**2)
      
      PREFCTI = ((CENTER(1)-CNTMU(1))**LI)*((CENTER(2)-CNTMU(2))**MI)*
     &          ((CENTER(3)-CNTMU(3))**NI)
      PREFCTJ = ((CENTER(1)-CNTNU(1))**LJ)*((CENTER(2)-CNTNU(2))**MJ)*
     &          ((CENTER(3)-CNTNU(3))**NJ)
C
      EXPFCTI = DEXP(-1.0D0*EXP1*DISTNCEI)
      EXPFCTJ = DEXP(-1.0D0*EXP2*DISTNCEJ)

      TMP2(INDX, JNDX) = PREFCTI*EXPFCTI*PREFCTJ*EXPFCTJ
C
      RETURN
      END



























