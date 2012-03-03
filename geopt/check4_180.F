      SUBROUTINE CHECK4_180(QRIC, DELTAQ, TOTNOFBND, TOTNOFANG,
     &                      TOTREDUNCO, NO_ITER)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER TOTNOFBND, TOTNOFANG, TOTREDUNCO
      LOGICAL NO_ITER
      DIMENSION QRIC(TOTREDUNCO), DELTAQ(TOTREDUNCO)
C
       PI = (ATAN(DFLOAT(1))*DFLOAT(4))
      TWOPI = (ATAN(DFLOAT(1))*DFLOAT(4))
C
#ifdef _DEBUG_LVL0
      Write(6,*) ' In CHECK4_180...'
      Print*, "TOTNFBND, TOTNOFANG: ",
     &         TOTNOFBND, TOTNOFANG
      Write(6,*)
      Write(6,*) "Angles @Entry:"
      Write(6,10) ((QRIC(I) + DELTAQ(I))/DINVPI, I = (TOTNOFBND + 1),
     &             TOTNOFANG+TOTNOFBND)
   10 Format (5(1X,F10.6)) 
#endif
C  
      DO I = (TOTNOFBND + 1), TOTNOFANG+TOTNOFBND
         DVAL = (QRIC(I) + DELTAQ(I)) - PI
         IF ( .NOT. ((DVAL .LE. 0.0D0) .OR. (DVAL .GT. 1.0D-8)) ) THEN
             NO_ITER = .TRUE.
             RETURN
         ENDIF
      ENDDO

C
      RETURN
      END
