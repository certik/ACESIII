

C I don't know who wrote this (or why) but we use it and we'll have to suffer
C the consequences. The old incarnation of this routine crashed the optimizer
C in Sun's compiler(s), but the new one doesn't. This is really a textbook
C example of how to fool optimizers and how not to write code. - ADY

      SUBROUTINE MODMATMUL(A,B,C,NA,NB,NC,NTA,NTB,NTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(NA,NB),C(NB,NC),A(NA,NC)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
      DOUBLE PRECISION SCRATCH(3*MxAtms,3*MxAtms)
      if (nta.gt.3*MxAtms) then
         print *, '@MODMATMUL: Assertion failed.'
         print *, '   scratch dimension = ',3*MxAtms
         print *, '   leading dimension = ',nta
         call aces_exit(1)
      end if
      DO J=1,NTC
         DO I=1,NTA
            SCRATCH(I,J)=0.D0
            DO K=1,NTB
               SCRATCH(I,J)=SCRATCH(I,J)+B(I,K)*C(K,J)
            END DO
         END DO
      END DO
      DO J=1,NTC
         DO I=1,NTA
            A(I,J)=SCRATCH(I,J)
         END DO
      END DO
      RETURN
      END

C       SUBROUTINE MODMATMUL(A,B,C,NA,NB,NC,NTA,NTB,NTC)
C C
C CJDW 1/6/98. This used to be called MATMUL. Name was changed since
C C            apparently MATMUL is a reserved name in Fortran 90. Note
C C            that it is not trivial to replace MATMUL by XGEMM since in
C C            BUILDB MATMUL (MODMATMUL) is called with A and C sharing same
C C            memory location.
C C
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C C     Maximum number of atoms currently allowed
C #include "mxatms.par"
C       DIMENSION B(NA,NB),C(NB,NC),A(NA,NC),SCRATCH(3*MxAtms,3*MxAtms)
C       DO 10 I=1,NTA
C       DO 10 J=1,NTC
C       Z=0.D0
C       DO 20 K=1,NTB
C 20    Z=Z+B(I,K)*C(K,J)
C       SCRATCH(I,J)=Z
C 10    CONTINUE
C       DO 15 I=1,NTA
C       DO 15 J=1,NTC
C 15    A(I,J)=SCRATCH(I,J)
C       RETURN
C       END
