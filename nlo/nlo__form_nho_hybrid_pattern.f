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
         SUBROUTINE  NLO__FORM_NHO_HYBRID_PATTERN
     +
     +                    ( NSHELL,
     +                      W,
     +
     +                             LENGTH,
     +                             HYBRID )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_NHO_HYBRID_PATTERN
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms the NHO hybridization pattern
C                from an input weight 'array' containing the sum of
C                the squares of the NHO expansion coefficients in
C                terms of NAOs for each angular momentum quantum
C                number. The pattern is formed as a long character
C                word and returned.
C
C                Procedure:
C
C                Consider the weight array W, where the component W (0)
C                has the sum of the squares of the NHO expansion
C                coefficients corresponding to the s-shell, W (1) the
C                one for the p-shell, ... and W (7) the one for all
C                shells > i. Our goal is to create the hybridization
C                pattern of the form:
C
C                      ab(x.xx)c(x.xx)d(x.xx)....
C
C                where the x.xx are the ratios of occurance of
C                the higher shells b,c,d,... in the NHO as compared
C                to the a-shell. The a-shell is defined as the
C                shell with the lowest weight > small threshold
C                weight. Hence all x.xx are >= 1.00. When choosing
C                the a-shell, it might happen that some of the
C                ratios x.xx are >= 10.00, in which case the above
C                pattern cannot be created. In this case we take
C                the maximum of all x.xx which have >= 10.00 and
C                define the a-shell as the one corresponding to
C                that maximum. The x.xx for the other shells are
C                then rescaled.
C
C                  Input:
C
C                    W       =  shell weight array
C
C                  Output:
C
C                    LENGTH  =  character length of hybridization
C                               pattern.
C                    HYBRID  =  character word containing hybridization
C                               pattern.
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

         CHARACTER*50   HYBRID
         CHARACTER*1    LMAIN

         CHARACTER*7   LRATIO  (1:7)
         CHARACTER*1   LSYMB   (0:7)

         LOGICAL     ABSOLUT
         LOGICAL     INCRESE

         INTEGER     I
         INTEGER     IBEG,IEND
         INTEGER     LENGTH
         INTEGER     NR
         INTEGER     NSHELL,SHELL

         INTEGER     SORT (1:8)

         DOUBLE PRECISION  WEIGHT,WLARGE,WSMALL,WLIMIT
         DOUBLE PRECISION  WMIN,WMAX
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  W (0:NSHELL)

         DATA  LSYMB   /'s','p','d','f','g','h','i','x'/

         PARAMETER  (ZERO   = 0.D0)
         PARAMETER  (ONE    = 1.D0)
         PARAMETER  (WLARGE = 1.D+1)
         PARAMETER  (WSMALL = 1.D-2)
         PARAMETER  (WLIMIT = 1.D-8)
C
C
C------------------------------------------------------------------------
C
C
C             ...find minimum weight > small limit threshold.
C
C
         WMIN = ONE
         DO 100 I = 0,NSHELL
            WEIGHT = W (I)
            IF (WEIGHT.GT.WLIMIT) THEN
                WMIN = DMIN1 (WMIN,WEIGHT)
            END IF
  100    CONTINUE
C
C
C             ...rescale the weights and find the maximum of all
C                rescaled weights >= 10.00.
C
C
         WMAX = ZERO
         DO 200 I = 0,NSHELL
            WEIGHT = W (I) / WMIN
            IF (WEIGHT.GE.WLARGE) THEN
                WMAX = DMAX1 (WMAX,WEIGHT)
            END IF
            W (I) = WEIGHT
  200    CONTINUE

         INCRESE = .TRUE.
C
C
C             ...if any rescaled weight >= 10.00 was found, rescale
C                weights once again.
C
C
         IF (WMAX.NE.ZERO) THEN
             DO 300 I = 0,NSHELL
                W (I) = W (I) / WMAX
  300        CONTINUE
             INCRESE = .FALSE.
         END IF
C
C
C             ...sort the shell indices appropriately.
C
C
         ABSOLUT = .FALSE.

         CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +              ( NSHELL+1,NSHELL+1,
     +                1,NSHELL+1,
     +                ABSOLUT,INCRESE,
     +                0,
     +                W,
     +
     +                        SORT )
     +
     +
C
C
C             ...eliminate those weights which are < WSMALL for
C                printing. These 'zeros' will be found at the
C                beginning or at the end of the sorting array,
C                depending on the ordering (increasing or decreasing)
C                of the weights.
C
C
         IF (INCRESE) THEN
             IBEG = 1
             IEND = NSHELL+1
             DO 400 I = 1,NSHELL+1
                SHELL = SORT (I) - 1
                WEIGHT = W (SHELL)
                IF (WEIGHT.LT.WSMALL) THEN
                    IBEG = IBEG + 1
                END IF
  400        CONTINUE
         ELSE
             IBEG = 1
             IEND = NSHELL+1
             DO 410 I = NSHELL+1,1,-1
                SHELL = SORT (I) - 1
                WEIGHT = W (SHELL)
                IF (WEIGHT.LT.WSMALL) THEN
                    IEND = IEND - 1
                END IF
  410        CONTINUE
         END IF
C
C
C             ...form the hybridization pattern. Remember, that
C                the sorted indices will be in the range 1 to NSHELL+1,
C                rather than 0 to NSHELL.
C
C
         SHELL = SORT (IBEG) - 1
         LMAIN = LSYMB (SHELL)

         NR = 0
         DO 500 I = IBEG+1,IEND
            SHELL = SORT (I) - 1
            WEIGHT = W (SHELL)
            NR = NR + 1
            WRITE (LRATIO (NR),9000) LSYMB (SHELL),'(',WEIGHT,')'
  500    CONTINUE

         WRITE (HYBRID,9100) LMAIN,(LRATIO(I),I=1,NR)

         LENGTH = 1 + 7 * NR
C
C
C             ...writing formats used.
C
C
 9000    FORMAT (A1,A1,F4.2,A1)
 9100    FORMAT (A1,7(A7))
C
C
C             ...ready!
C
C
         RETURN
         END
