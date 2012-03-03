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
         INTEGER FUNCTION  NLO__ANGULAR_MOMENTUM_NBO
     +
     +                          ( NANG,
     +                            ANGMOM,
     +                            FIRST,LAST,
     +                            L,
     +                            NDEG,
     +                            QSYMACC,LSYMACC,
     +                            W,Q,
     +                            LOCAL )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ANGULAR_MOMENTUM_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This function returns the first index position of
C                an NBO set inside the ANGMOM array, all of which
C                posses lowest angular momentum component L. If
C                there is any negative value in ANGMOM, the NBO will
C                not participate in the search, i.e. will not be
C                recognized. Three criteria must be met by the set
C                of NBOs chosen:
C
C                     1) the lowest angular momentum must
C                        have a value of L
C
C                     2) there must be NDEG NBOs having the
C                        same lowest angular momentum value L,
C                        where the degeneracy criterion is based
C                        on equal interaction order and locality
C                        content within a threshold limit.
C
C                     3) off all possible sets found obeying 1)
C                        and 2) we pick the one that has the closest
C                        weight to the target weight corresponding
C                        to the weight average of the excluded
C                        NBOs between the ranges FIRST and LAST.
C
C                If no NBO set is found, a value of 0 is returned.
C                Since the weights of all NBOs transmitted here are
C                already ordered, the search will be performed starting
C                from both ends -> and <- and move inwards to the
C                positions FIRST and LAST inside the ANGMOM array:
C
C
C                      * * * * * * * * * * * * * * * * * * * * *
C                      ->        |       |                    <-
C                              FIRST    LAST
C
C                This procedure will ensure close proximity to the
C                target weight in FIRST and LAST.
C
C
C
C                  Input:
C
C                    NANG         =  # of lowest angular momenta
C                                    present in array ANGMOM over which
C                                    the search will be performed
C                    ANGMOM (I)   =  lowest angular momentum component
C                                    of I-th NBO.
C                    FIRST,LAST   =  the range of excluded NBOs from
C                                    the search.
C                    L            =  lowest angular momentum component
C                                    of NBO wanted.
C                    NDEG         =  # of consecutive NBOs wanted with
C                                    the same weight and locality
C                                    content within accuracy limits.
C                    QSYMACC      =  threshold value for NBO interaction
C                                    order equalities.
C                    LSYMACC      =  threshold value for NBO locality
C                                    equalities.
C                    W            =  the NANG NBO weights.
C                    Q            =  the NANG NBO interaction orders.
C                    LOCAL        =  the NANG NBO atomic localization
C                                    content.
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

         LOGICAL     ACCEPT

         INTEGER     K,L,M,N
         INTEGER     FIRST,LAST
         INTEGER     NANG
         INTEGER     NDEG
         INTEGER     NFIRST,NLAST
         INTEGER     NSTEPS

         INTEGER     ANGMOM  (1:NANG)

         DOUBLE PRECISION  QDIFF,LDIFF
         DOUBLE PRECISION  QSYMACC,LSYMACC
         DOUBLE PRECISION  QVAL,LVAL
         DOUBLE PRECISION  WLAST,WFIRST
         DOUBLE PRECISION  WTARGET
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  LOCAL (1:NANG)
         DOUBLE PRECISION  Q     (1:NANG)
         DOUBLE PRECISION  W     (1:NANG)

         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...form target (average) weight.
C
C
         WTARGET = ZERO
         DO N = FIRST,LAST
            WTARGET = WTARGET + W (N)
         END DO
         WTARGET = WTARGET / DFLOAT (N)
C
C
C             ...determine closest NBO set first index positions
C                NLAST and NFIRST from both respective LAST and
C                FIRST ends.
C
C
         NLAST = 0
         NFIRST = 0

         NSTEPS = MAX (FIRST-1,NANG-LAST)

         DO N = NSTEPS,1,-1
            M = LAST + N
            IF (M.LE.NANG) THEN
                IF (ANGMOM (M).EQ.L) THEN
                    QVAL = Q (M)
                    LVAL = LOCAL (M)
                    ACCEPT = .TRUE.
                    DO K = M-NDEG+1,M-1
                       QDIFF  = ABS (QVAL - Q (K))
                       LDIFF  = ABS (LVAL - LOCAL (K))
                       ACCEPT = ACCEPT .AND. (K.GT.LAST)
                       ACCEPT = ACCEPT .AND. (ANGMOM (K).EQ.L)
                       ACCEPT = ACCEPT .AND. (QDIFF.LT.QSYMACC)
                       ACCEPT = ACCEPT .AND. (LDIFF.LT.LSYMACC)
                    END DO
                    IF (ACCEPT) THEN
                        NLAST = M-NDEG+1
                    END IF
                END IF
            END IF
            M = FIRST - N
            IF (M.GE.1) THEN
                IF (ANGMOM (M).EQ.L) THEN
                    QVAL = Q (M)
                    LVAL = LOCAL (M)
                    ACCEPT = .TRUE.
                    DO K = M+1,M+NDEG-1
                       QDIFF  = ABS (QVAL - Q (K))
                       LDIFF  = ABS (LVAL - LOCAL (K))
                       ACCEPT = ACCEPT .AND. (K.LT.FIRST)
                       ACCEPT = ACCEPT .AND. (ANGMOM (K).EQ.L)
                       ACCEPT = ACCEPT .AND. (QDIFF.LT.QSYMACC)
                       ACCEPT = ACCEPT .AND. (LDIFF.LT.LSYMACC)
                    END DO
                    IF (ACCEPT) THEN
                        NFIRST = M
                    END IF
                END IF
            END IF
         END DO
C
C
C             ...decide which to take by comparing with target
C                weight.
C
C
         IF (NLAST.NE.0 .AND. NFIRST.NE.0) THEN

             WLAST  = ABS (WTARGET - W (NLAST))
             WFIRST = ABS (WTARGET - W (NFIRST))

             IF (WLAST.LT.WFIRST) THEN
                 NLO__ANGULAR_MOMENTUM_NBO = NLAST
             ELSE
                 NLO__ANGULAR_MOMENTUM_NBO = NFIRST
             END IF

         ELSE IF (NLAST.NE.0) THEN

             NLO__ANGULAR_MOMENTUM_NBO = NLAST

         ELSE IF (NFIRST.NE.0) THEN

             NLO__ANGULAR_MOMENTUM_NBO = NFIRST

         ELSE

             NLO__ANGULAR_MOMENTUM_NBO = 0

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
