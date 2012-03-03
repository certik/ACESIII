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
         SUBROUTINE  NLO__CHECK_LINEAR_DEPENDENCY
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWS,DDCOLS,
     +                      ROW,COL,
     +                      X,S,
     +
     +                              DEPEND )
     +
C-----------------------------------------------------------------------
C  OPERATION   : NLO__CHECK_LINEAR_DEPENDENCY
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine checks, if a system of COL columns
C                is linearly dependent or not. The routine simply
C                calculates the overlap matrix S between these
C                columns and attempts a Cholesky decomposition on
C                S. Any linear dependency will show up as a very
C                small (near zero) diagonal element during the
C                Cholesky decomposition. As soon as the routine
C                encounters such a situation we exit with the
C                linear dependency indicator DEPEND set true.
C
C                The Cholesky decomposition is performed in place,
C                that is the originally calculated overlap matrix S
C                will be destroyed. The column matrix X is conserved.
C
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     DEPEND

         INTEGER     DDROWX,DDCOLX,DDROWS,DDCOLS
         INTEGER     I,J,K
         INTEGER     ROW,COL

         DOUBLE PRECISION   ROOT
         DOUBLE PRECISION   SQRARG
         DOUBLE PRECISION   SUM
         DOUBLE PRECISION   VSMALL
         DOUBLE PRECISION   ZERO

         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)
         DOUBLE PRECISION  S (1:DDROWS,1:DDCOLS)

         DATA   ZERO      /0.D0/
         DATA   VSMALL    /1.D-8/
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of X and S matrices supplied.
C
C
         IF (ROW.GT.DDROWX .OR. COL.GT.DDCOLX) THEN
             WRITE (*,*) ' Dimensions of matrix X too small! '
             WRITE (*,*) ' nlo__check_linear_dependency '
             WRITE (*,*) ' DDROWX,DDCOLX,ROW,COL = ',
     +                     DDROWX,DDCOLX,ROW,COL
             WRITE (1,*) ' Dimensions of matrix X too small! '
             WRITE (1,*) ' nlo__check_linear_dependency '
             WRITE (1,*) ' DDROWX,DDCOLX,ROW,COL = ',
     +                     DDROWX,DDCOLX,ROW,COL
             STOP
         END IF

         IF (COL.GT.DDROWS .OR. COL.GT.DDCOLS) THEN
             WRITE (*,*) ' Dimensions of matrix S too small! '
             WRITE (*,*) ' nlo__check_linear_dependency '
             WRITE (*,*) ' DDROWS,DDCOLS,COL = ',DDROWS,DDCOLS,COL
             WRITE (1,*) ' Dimensions of matrix S too small! '
             WRITE (1,*) ' nlo__check_linear_dependency '
             WRITE (1,*) ' DDROWS,DDCOLS,COL = ',DDROWS,DDCOLS,COL
             STOP
         END IF
C
C
C             ...calculate lower triangle of overlap matrix.
C
C
         DO 100 J = 1,COL
         DO 100 I = J,COL
            SUM = ZERO
            DO 110 K = 1,ROW
               SUM = SUM + X (K,I) * X (K,J)
  110       CONTINUE
            S (I,J) = SUM
  100    CONTINUE
C
C
C             ...do the Cholesky decomposition on S.
C
C
         DEPEND = .FALSE.

         IF (S (1,1).LE.VSMALL) THEN
             DEPEND = .TRUE.
             RETURN
         END IF

         ROOT  =  DSQRT ( S (1,1) )
         S (1,1)  =  ROOT

         DO 200 I = 2,COL
            S (I,1)  =  S (I,1) / ROOT
  200    CONTINUE

         DO 210 J = 2,COL

            SQRARG  =  S (J,J)
            DO 220  K = 1,J-1
               SQRARG  =  SQRARG  -  S (J,K) ** 2
  220       CONTINUE

            IF (SQRARG.LE.VSMALL) THEN
                DEPEND = .TRUE.
                RETURN
            END IF

            ROOT  =  DSQRT (SQRARG)
            S (J,J)  =  ROOT

            DO 230 I = J+1,COL
               SUM  =  S (I,J)
               DO 240 K = 1,J-1
                  SUM  =  SUM  -  S (I,K) * S (J,K)
  240          CONTINUE
               S (I,J)  =  SUM / ROOT
  230       CONTINUE

  210    CONTINUE
C
C
C             ...finished!
C
C
         RETURN
         END
