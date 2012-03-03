C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
C I don't know who wrote this (or why) but we use it and we'll have to suffer
C the consequences. The old incarnation of this routine crashed the optimizer
C in Sun's compiler(s), but the new one doesn't. This is really a textbook
C example of how to fool optimizers and how not to write code. - ADY

      SUBROUTINE MODMATMUL(A,B,C,NA,NB,NC,NTA,NTB,NTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(NA,NB),C(NB,NC),A(NA,NC)
      DOUBLE PRECISION SCRATCH(nta,ntc)

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
