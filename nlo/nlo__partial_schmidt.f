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
         SUBROUTINE  NLO__PARTIAL_SCHMIDT
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      DDROWS,DDCOLS,
     +                      ROW,COLC,COLD,
     +                      COLMAX,COLMIN,
     +                      S,
     +                      SAVES,SAVEC,SAVED,
     +                      C,
     +                      XVEC,
     +                      XMAT,
     +
     +                              D )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PARTIAL_SCHMIDT
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine orthogonalizes a set of MO vectors to
C                an already orthogonalized set of MO vectors using
C                the Schmidt orthogonalization method. The overlaps
C                needed between the two sets are evaluated from the
C                MO expansion coefficients in terms of AO's and the
C                AO overlap matrix.
C
C                The Schmidt orthogonalization leads to the following
C                expression of the new orthogonal D vectors:
C
C
C                                       M1
C                        D(i) = D(i) - sum S(k,i) C(k)
C                                       k
C
C                where i refers to the initially non-orthogonal D
C                functions and S (k,i) denotes the MO overlap matrix
C                element between the k-th orthogonal C function and
C                the i-th D function.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     SAVES,SAVEC,SAVED

         INTEGER     COLMAX,COLMIN
         INTEGER     DDROWC,DDCOLC,DDROWD,DDCOLD,DDROWS,DDCOLS
         INTEGER     I,J,K
         INTEGER     ROW,COLC,COLD

         DOUBLE PRECISION  X

         DOUBLE PRECISION  XVEC  (1:COLMAX)

         DOUBLE PRECISION  C     (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  D     (1:DDROWD,1:DDCOLD)
         DOUBLE PRECISION  S     (1:DDROWS,1:DDCOLS)
         DOUBLE PRECISION  XMAT  (1:ROW,   1:COLMIN)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of C,D and S matrices supplied.
C
C
         IF (ROW.GT.DDROWC .OR. COLC.GT.DDCOLC) THEN
             WRITE (*,*) ' Dimensions of matrix C too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' DDROWC,DDCOLC,ROW,COLC = ',
     +                     DDROWC,DDCOLC,ROW,COLC
             WRITE (1,*) ' Dimensions of matrix C too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' DDROWC,DDCOLC,ROW,COLC = ',
     +                     DDROWC,DDCOLC,ROW,COLC
             STOP
         END IF

         IF (ROW.GT.DDROWD .OR. COLD.GT.DDCOLD) THEN
             WRITE (*,*) ' Dimensions of matrix D too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' DDROWD,DDCOLD,ROW,COLD = ',
     +                     DDROWD,DDCOLD,ROW,COLD
             WRITE (1,*) ' Dimensions of matrix D too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' DDROWD,DDCOLD,ROW,COLD = ',
     +                     DDROWD,DDCOLD,ROW,COLD
             STOP
         END IF

         IF (ROW.GT.DDROWS .OR. ROW.GT.DDCOLS) THEN
             WRITE (*,*) ' Dimensions of matrix S too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' DDROWS,DDCOLS,ROW = ',DDROWS,DDCOLS,ROW
             WRITE (1,*) ' Dimensions of matrix S too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' DDROWS,DDCOLS,ROW = ',DDROWS,DDCOLS,ROW
             STOP
         END IF
C
C
C             ...check dimensions of working arrays XVEC and XMAT.
C
C
         IF (COLMAX.LT.MAX (COLC,COLD)) THEN
             WRITE (*,*) ' Dimensions of array XVEC too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' COLMAX,MAX (COLC,COLD) = ',
     +                     COLMAX,MAX (COLC,COLD)
             WRITE (1,*) ' Dimensions of array XVEC too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' COLMAX,MAX (COLC,COLD) = ',
     +                     COLMAX,MAX (COLC,COLD)
             STOP
         END IF

         IF (COLMIN.LT.MIN (COLC,COLD)) THEN
             WRITE (*,*) ' Dimensions of array XMAT too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' COLMIN,MIN (COLC,COLD) = ',
     +                     COLMIN,MIN (COLC,COLD)
             WRITE (1,*) ' Dimensions of array XMAT too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' COLMIN,MIN (COLC,COLD) = ',
     +                     COLMIN,MIN (COLC,COLD)
             STOP
         END IF
C
C
C             ...proceed according to sizes of C and D. If the size of
C                C is smaller than D, then we generate the COLD x COLC
C                molecular orbital overlap matrix, otherwise we generate
C                its transpose of COLC x COLD size.
C
C
         IF (COLC.LT.COLD) THEN

             CALL    MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                    ( DDROWD,DDCOLD,
     +                      DDROWS,DDCOLS,
     +                      DDROWC,DDCOLC,
     +                      ROW,COLC,
     +                      COLD,
     +                      COLD,ROW,COLC,
     +                      0,0,
     +                      SAVED,SAVES,SAVEC,
     +                      D,S,C,
     +                      XVEC,
     +
     +                             XMAT )
     +
     +
             DO I = 1,COLC
             DO J = 1,COLD
                X = XMAT (J,I)
                DO K = 1,ROW
                   D (K,J) = D (K,J) - X * C (K,I)
                END DO
             END DO
             END DO

         ELSE

             CALL    MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWS,DDCOLS,
     +                      DDROWD,DDCOLD,
     +                      ROW,COLD,
     +                      COLC,
     +                      COLC,ROW,COLD,
     +                      0,0,
     +                      SAVEC,SAVES,SAVED,
     +                      C,S,D,
     +                      XVEC,
     +
     +                             XMAT )
     +
     +
             DO J = 1,COLD
             DO I = 1,COLC
                X = XMAT (I,J)
                DO K = 1,ROW
                   D (K,J) = D (K,J) - X * C (K,I)
                END DO
             END DO
             END DO

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
