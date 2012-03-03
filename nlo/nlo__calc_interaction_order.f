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
         SUBROUTINE  NLO__CALC_INTERACTION_ORDER
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWP,DDCOLP,
     +                      DDVECD,
     +                      DDVECX,DDVECY,
     +                      N,
     +                      I,
     +                      C,P,D,
     +                      X,Y,
     +
     +                              ORDER )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__CALC_INTERACTION_ORDER
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine calculates the occupation interaction
C                order of the I-th natural orbital D. The occupation
C                interaction order is defined as follows:
C
C                The N x N occupation matrix P is idempotent:
C
C                                 P = P * P
C
C                For any diagonal element in P we have thus:
C
C                                                N
C                    P (I,I) - P (I,I)**2  =    sum   P (I,K)**2
C                                             i neq k
C
C                The square root of the right sum constitutes the
C                occupation interaction order of the I-th natural
C                orbital. If the individual P (I,K) are accurate
C                to +- E, where E is the accuracy threshold, then
C                we have:
C
C
C                     ---------------------------
C                    /   N
C                   /   sum   [P (I,K) +- E ]**2
C                 \/  i neq k
C
C
C                              -------------------
C                             /   N
C                      =     /   sum   P (I,K)**2
C                          \/  i neq k
C
C
C                                       N
C                                      sum   P (I,K)
C                                    i neq k                         2
C                      +- E   ------------------------------  + O [E]
C                                   -------------------
C                                  /   N
C                                 /   sum   P (I,K)**2
C                               \/  i neq k
C
C
C                showing that the occupation interaction order is
C                also roughly accurate to within +- E.
C
C
C
C                  Input:
C
C                    DDROWC,DDCOLC  =  declared dimensions of natural
C                                      orbital coefficient matrix C
C                    DDROWP,DDCOLP  =  declared dimensions of occupation
C                                      matrix P
C                    DDVECD         =  declared dimensions of natural
C                                      orbital D coefficient vector
C                    DDVECX,DDVECY  =  declared dimensions of scratch
C                                      vectors X and Y
C                    N              =  order of occupation matrix P
C                    I              =  index label of natural orbital D
C                                      (this is necessary to remove the
C                                      I-th contribution to the sum)
C                    C              =  N x N natural orbital coefficient
C                                      matrix in AO basis
C                    P              =  N x N occupation matrix in AO
C                                      basis
C                    D              =  natural orbital coefficient
C                                      vector in AO basis
C                    X,Y            =  flp scratch vectors
C
C
C                  Output:
C
C                    ORDER          =  occupation interaction order
C                                      for I-th natural orbital D
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

         LOGICAL     SAVEC,SAVED,SAVEP

         INTEGER     DDROWC,DDCOLC
         INTEGER     DDROWP,DDCOLP
         INTEGER     DDVECD,DDVECX,DDVECY
         INTEGER     I,J,N

         DOUBLE PRECISION  ORDER
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  D     (1:DDVECD)
         DOUBLE PRECISION  X     (1:DDVECX)
         DOUBLE PRECISION  Y     (1:DDVECY)

         DOUBLE PRECISION  C     (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  P     (1:DDROWP,1:DDCOLP)

         PARAMETER   (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF (N.GT.DDROWC .OR. N.GT.DDCOLC) THEN
             WRITE (*,*) ' Dimensions of matrix C too small! '
             WRITE (*,*) ' nlo__calc_interaction_order '
             WRITE (*,*) ' DDROWC,DDCOLC,N = ',DDROWC,DDCOLC,N
             WRITE (1,*) ' Dimensions of matrix C too small! '
             WRITE (1,*) ' nlo__calc_interaction_order '
             WRITE (1,*) ' DDROWC,DDCOLC,N = ',DDROWC,DDCOLC,N
             STOP
         END IF

         IF (N.GT.DDROWP .OR. N.GT.DDCOLP) THEN
             WRITE (*,*) ' Dimensions of matrix P too small! '
             WRITE (*,*) ' nlo__calc_interaction_order '
             WRITE (*,*) ' DDROWP,DDCOLP,N = ',DDROWP,DDCOLP,N
             WRITE (1,*) ' Dimensions of matrix P too small! '
             WRITE (1,*) ' nlo__calc_interaction_order '
             WRITE (1,*) ' DDROWP,DDCOLP,N = ',DDROWP,DDCOLP,N
             STOP
         END IF

         IF (N.GT.DDVECD) THEN
             WRITE (*,*) ' Dimension of vector D too small! '
             WRITE (*,*) ' nlo__calc_interaction_order '
             WRITE (*,*) ' DDVECD,N = ',DDVECD,N
             WRITE (1,*) ' Dimension of vector D too small! '
             WRITE (1,*) ' nlo__calc_interaction_order '
             WRITE (1,*) ' DDVECD,N = ',DDVECD,N
             STOP
         END IF

         IF (N.GT.DDVECX .OR. N.GT.DDVECY) THEN
             WRITE (*,*) ' Dimensions of vectors X,Y too small! '
             WRITE (*,*) ' nlo__calc_interaction_order '
             WRITE (*,*) ' DDVECX,DDVECY,N = ',DDVECX,DDVECY,N
             WRITE (1,*) ' Dimensions of vectors X,Y too small! '
             WRITE (1,*) ' nlo__calc_interaction_order '
             WRITE (1,*) ' DDVECX,DDVECY,N = ',DDVECX,DDVECY,N
             STOP
         END IF
C
C
C             ...form the I-th column of the occupation matrix
C                corresponding to the natural orbital D.
C
C
         SAVEC = .TRUE.
         SAVEP = .TRUE.
         SAVED = .TRUE.

         CALL  MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +              ( DDROWC,DDCOLC,
     +                DDROWP,DDCOLP,
     +                DDVECD,1,
     +                DDVECY,1,
     +                DDVECX,
     +                N,N,1,
     +                0,0,
     +                SAVEC,SAVEP,SAVED,
     +                C,P,D,
     +                X,
     +
     +                         Y )
     +
     +
C
C
C             ...calculate the interaction order.
C
C
         Y (I) = ZERO

         ORDER = ZERO
         DO J = 1,N
            ORDER = ORDER + Y (J) ** 2
         END DO

         ORDER = DSQRT (ORDER)
C
C
C             ...ready!
C
C
         RETURN
         END
