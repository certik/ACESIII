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
         SUBROUTINE  NLO__FIND_ROTATION_NBO_INDICES
     +
     +                    ( NBAS,NATOM,NDEG,
     +                      NBOSIZE,
     +                      OFFLP,OFFBD,OFFEP,
     +                      OFFAB,OFFRY,OFFRR,
     +                      CENTER,
     +                      SYMSIZE,SYMSETS,SYMCEN,
     +                      SYMNBO,
     +                      QSYMACC,LSYMACC,
     +                      BDNCEN,BDCEN,
     +                      NBOBD,
     +                      SYMCRIT,
     +                      PSUB,
     +                      LOCAL,
     +                      Q,
     +
     +                              ZEROSYM,
     +                              NROT,
     +                              ROTIDX )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FIND_ROTATION_NBO_INDICES
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine determines those NBO indices which will
C                be used for reorientation of axial and central
C                degenerate atomic NBOs.
C
C                Procedure:
C
C                We will only use those NBOs for reorientation which
C                are actually marked as usable for that purpose in the
C                SYMNBO array (see below). The following algorithm
C                is applied:
C
C                1) Determination of NBO subset {NBO} for rotation
C
C                   Loop over all NBOs
C                     Check NBO => is SYMNBO (NBO) = 2 ?
C                       If yes, get NBO atomic indices NIDX
C                         If NIDX = 2 (2-cen bond or antibond) then 
C                           Is one atomic index = rotation center ?
C                             If yes, is second index member any
C                             of the symmetry center sets ?
C                                If yes, add present NBO to set {NBO}
C                         If NIDX = 1 (single atom NBO) then
C                           Is atomic index member any of the
C                           symmetry center sets ?
C                              If yes, add present NBO to set {NBO}
C                   continue
C
C                2) Determination of final NBO set [NBO] for rotation
C
C                      i) determine the largest absolute element NBOx
C                         of all NDEG degenerate PSUB vectors
C                         within the {NBO} subset determined
C                     ii) if this element < SYMCRIT, exit routine
C                         with ZEROSYM = .true.
C                    iii) determine the atomic NBO index IDX different
C                         from the rotation center to which that
C                         element NBOx belongs
C                     iv) determine what type of NBO it belongs to
C                         (Core, Lone-pair, 2-cen Bond, etc ...)
C                      v) determine all atomic index partners [IDX]
C                         corresponding to IDX and the symmetry
C                         center sets
C                     vi) determine all those NBO indices [NBO] which
C                         belong to [IDX], have the same type and are
C                         symmetry related to NBOx
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    NDEG         =  # of degenerate NBOs. Can only
C                                    be equal to 2,3,4 or 5.
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    OFFxy        =  absolute offset indices for x-type
C                                    NBOs. These indices indicate the
C                                    totallity of x-type NBOs before the
C                                    x-type NBO set (xy:LP=Lone-pair,
C                                    BD=Bond, EP=Empty-pair, AB=Anti-
C                                    bond, RY,RR=Rydberg types)
C                    CENTER       =  atomic index on which NBO rotation
C                                    will be performed
C                    SYMSIZE      =  the size of the symmetry related
C                                    atomic index sets (rings or high
C                                    symmetry sets)
C                    SYMSETS      =  the # of symmetry related atomic
C                                    atomic index sets (rings or high
C                                    symmetry sets)
C                    SYMCEN (A)   =  all SYMSIZE * SYMSETS symmetry
C                                    related atomic indices A grouped
C                                    into sets of size SYMSIZE
C                    SYMNBO (I)   =  integer vector indicating if the
C                                    I-th NBO is considered to be used
C                                    for symmetry adaptation (=2 if to
C                                    be used).
C                    QSYMACC      =  threshold value for NBO interaction
C                                    order equalities.
C                    LSYMACC      =  threshold value for NBO locality
C                                    equalities.
C                    BDNCEN (J)   =  current # of atomic centers for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDCEN (I,J)  =  current I-th atomic center index
C                                    for J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    for the I-th atomic x-type NBO.
C                                    This array is the handle for
C                                    accessing and modifying info of
C                                    the NHO bonds sitting in the
C                                    arrays BDNCEN,BDCEN,BDBAS and
C                                    BDOCC.
C                    SYMCRIT      =  symmetry interaction criterion.
C                                    Any absolute occupation matrix
C                                    interaction element between
C                                    symmetrized NBOs and a degenerate
C                                    set of NBOs is assumed to be zero
C                                    if it falls below this value.
C                    PSUB         =  NBAS x NDEG section of the
C                                    occupation matrix.
C                    LOCAL        =  atomic localization content for
C                                    all NBOs.
C                    Q            =  complete NBO interaction order
C                                    vector at present stage of
C                                    symmetrization.
C
C
C                  Output:
C
C                    ZEROSYM      =  is false, if the search has been
C                                    successful in the sense that at
C                                    least one occupation matrix
C                                    element was found to be larger or
C                                    equal than the symmetry interaction
C                                    criterion value SYMCRIT. If true,
C                                    then all occupation matrix elements
C                                    were < SYMCRIT, an indication
C                                    that the current set of NBOs passed
C                                    here have no symmetry interaction
C                                    with those NBOs currently available
C                                    for symmetry analysis.
C                    NROT         =  # of NBO indices that will be used
C                                    for rotation.
C                    ROTIDX       =  the rotation NBO indices.
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

         LOGICAL     ATOMIC
         LOGICAL     QSYM,LSYM
         LOGICAL     ZEROSYM

         INTEGER     ATX,ATX1,ATX2,ATXOLD,ATXMAX
         INTEGER     BOND
         INTEGER     CENTER
         INTEGER     J,N
         INTEGER     MROT,NROT
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NBO,NBO1ST,NBOLST,NBOBEG,NBOEND,NBOSIZE,NBOMAX
         INTEGER     NCEN
         INTEGER     NDEG
         INTEGER     NITER,MXITER
         INTEGER     NSYM,NSYMOLD
         INTEGER     NXBA
         INTEGER     OFF,OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,OFFRR
         INTEGER     ROT
         INTEGER     SET,SETMAX
         INTEGER     SYM,SYM1ST,SYMLST,SYMSIZE,SYMSETS

         INTEGER     BDNCEN  (1:NBAS )
         INTEGER     NBOBD   (1:NBAS )
         INTEGER     ROTIDX  (1:NBAS )
         INTEGER     SYMCEN  (1:NATOM)
         INTEGER     SYMNBO  (1:NBAS )

         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS )

         DOUBLE PRECISION  PMAX,PVAL
         DOUBLE PRECISION  QMAX,LMAX
         DOUBLE PRECISION  QNBO,LNBO
         DOUBLE PRECISION  QSYMACC,LSYMACC
         DOUBLE PRECISION  QTHRESH,LTHRESH
         DOUBLE PRECISION  SYMCRIT
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  LOCAL  (1:NBAS)
         DOUBLE PRECISION  Q      (1:NBAS)

         DOUBLE PRECISION  PSUB  (1:NBAS,1:NDEG)

         PARAMETER  (ZERO   = 0.D0)
         PARAMETER  (MXITER = 50  )
C
C
C------------------------------------------------------------------------
C
C
C             ...loop over all NBOs and determine the NBO, atomic and
C                symmetry set index leading to the maximum PSUB value.
C
C
         PMAX = ZERO
         ATXMAX = 0
         NBOMAX = 0

         DO NBO = 1,NBAS
            IF (SYMNBO (NBO).EQ.2) THEN

                BOND = NBOBD (NBO)
                NCEN = BDNCEN (BOND)

                IF (NCEN.EQ.1) THEN
                    ATX = BDCEN (1,BOND)
                ELSE IF (NCEN.EQ.2) THEN
                    ATX  = 0
                    ATX1 = BDCEN (1,BOND)
                    ATX2 = BDCEN (2,BOND)
                    IF (ATX1.EQ.CENTER) THEN
                        ATX = ATX2
                    ELSE IF (ATX2.EQ.CENTER) THEN
                        ATX = ATX1
                    END IF
                END IF

                IF (ATX.NE.0) THEN
                    N = 0
                    DO SET = 1,SYMSETS
                    DO SYM = 1,SYMSIZE
                       N = N + 1
                       IF (ATX.EQ.SYMCEN (N)) THEN
                           DO J = 1,NDEG
                              PVAL = ABS (PSUB (NBO,J))
                              IF (PVAL.GT.PMAX) THEN
                                  ATXMAX = ATX
                                  NBOMAX = NBO
                                  SETMAX = SET
                                  PMAX = PVAL
                              END IF
                           END DO
                       END IF
                    END DO
                    END DO
                END IF

            END IF
         END DO

         ZEROSYM = ABS (PMAX) .LT. SYMCRIT

         IF (ZEROSYM) THEN
             RETURN
         END IF
C
C
C             ...determine type of found NBO index by narrowing
C                search within the complete NBO set. Determine
C                also the symmetry set range.
C
C
         IF (NBOMAX.EQ.0) THEN
             WRITE (*,*) ' No NBO index! '
             WRITE (*,*) ' NBO = ',NBO
             WRITE (*,*) ' nlo__find_rotation_nbo_indices '
             WRITE (1,*) ' No NBO index! '
             WRITE (1,*) ' NBO = ',NBO
             WRITE (1,*) ' nlo__find_rotation_nbo_indices '
             STOP
         ELSE IF (NBOMAX.LE.OFFLP) THEN
             NBO1ST = 1
             NBOLST = OFFLP
             ATOMIC = .TRUE.
         ELSE IF (NBOMAX.LE.OFFBD) THEN
             NBO1ST = OFFLP + 1
             NBOLST = OFFBD
             ATOMIC = .TRUE.
         ELSE IF (NBOMAX.LE.OFFEP) THEN
             NBO1ST = OFFBD + 1
             NBOLST = OFFEP
             ATOMIC = .FALSE.
         ELSE IF (NBOMAX.LE.OFFAB) THEN
             NBO1ST = OFFEP + 1
             NBOLST = OFFAB
             ATOMIC = .TRUE.
         ELSE IF (NBOMAX.LE.OFFRY) THEN
             NBO1ST = OFFAB + 1
             NBOLST = OFFRY
             ATOMIC = .FALSE.
         ELSE IF (NBOMAX.LE.OFFRR) THEN
             NBO1ST = OFFRY + 1
             NBOLST = OFFRR
             ATOMIC = .TRUE.
         ELSE
             NBO1ST = OFFRR + 1
             NBOLST = NBAS
             ATOMIC = .TRUE.
         END IF

         SYMLST = SETMAX * SYMSIZE
         SYM1ST = SYMLST - SYMSIZE + 1
C
C
C             ...at this point we have narrowed down the NBO type
C                and we know to which symmetry set the NBO belongs.
C                Determine next the other corresponding NBOs from
C                the symmetry set. Two posibilities arise:
C
C                  i) Atomic (i.e. 1 center) NBOs
C                     ---------------------------
C                     Here we know that all NBOs are clustered
C                     according to atomic index and each type
C                     within each atomic set occurs at the same
C                     relative position. Hence, in this case we
C                     determine the relative index of NBOMAX and
C                     find the symmetry related partners using
C                     this index. After that we need to include
C                     those which have the same interaction order
C                     and locality in their immediate neighborhood.
C
C                 ii) Non-Atomic (i.e. 2 center) NBOs
C                     -------------------------------
C                     Here we must find all symmetry related
C                     bonds within the bond type found. These
C                     bonds must be present due to symmetry and
C                     they have to obey the following criteria:
C                     a) same # of centers (i.e. 2), b) same
C                     NBO interaction order, c) same localization
C                     content.
C
C
         QMAX = Q     (NBOMAX)
         LMAX = LOCAL (NBOMAX)

         WRITE (*,*) ' NBOMAX = ',NBOMAX
         WRITE (*,*) ' PMAX = ',PMAX
         WRITE (*,*) ' QMAX,LMAX = ',QMAX,LMAX

         IF (ATOMIC) THEN

             OFF = 0
             NXBA = 0
             NROT = 0
             ATXOLD = 0

             WRITE (*,*) ' NBO1ST,NBOLST = ',NBO1ST,NBOLST

             DO NBO = NBO1ST,NBOLST

                BOND = NBOBD (NBO)
                NCEN = BDNCEN (BOND)

                IF (NCEN.NE.1) THEN
                    WRITE (*,*) ' # of centers mismatch! '
                    WRITE (*,*) ' ATOMIC,NCEN = ',ATOMIC,NCEN
                    WRITE (*,*) ' nlo__find_rotation_nbo_indices '
                    WRITE (1,*) ' # of centers mismatch! '
                    WRITE (1,*) ' ATOMIC,NCEN = ',ATOMIC,NCEN
                    WRITE (1,*) ' nlo__find_rotation_nbo_indices '
                    STOP
                END IF

                ATX = BDCEN (1,BOND)

                IF (ATX.EQ.ATXMAX) THEN
                    NXBA = NXBA + 1
                    IF (NBO.LT.NBOMAX) THEN
                        OFF = OFF + 1
                    END IF
                END IF

                IF (ATX.NE.ATXOLD) THEN
                    DO SYM = SYM1ST,SYMLST
                       IF (ATX.EQ.SYMCEN (SYM)) THEN
                           NROT = NROT + 1
                           ROTIDX (NROT) = NBO
                       END IF
                    END DO
                    ATXOLD = ATX
                END IF

             END DO

             IF (NROT.NE.SYMSIZE) THEN
                 WRITE (*,*) ' Bad # of NBO rotation indices! '
                 WRITE (*,*) ' SYMSIZE,MROT = ',SYMSIZE,MROT
                 WRITE (*,*) ' nlo__find_rotation_nbo_indices '
                 WRITE (1,*) ' Bad # of NBO rotation indices! '
                 WRITE (1,*) ' SYMSIZE,MROT = ',SYMSIZE,MROT
                 WRITE (1,*) ' nlo__find_rotation_nbo_indices '
                 STOP
             END IF

             WRITE (*,*) ' OFF = ',OFF
             WRITE (*,*) ' NXBA = ',NXBA

             DO ROT = 1,NROT
                WRITE (*,*) ' Rot index (intermediate ) = ',ROTIDX (ROT)
             END DO

             CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                    ( 6,
     +                      ' Interaction order vector ',
     +                      NBAS,1,
     +                      NBAS,1,
     +                      Q )
     +
     +
             CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                    ( 6,
     +                      ' Locality vector ',
     +                      NBAS,1,
     +                      NBAS,1,
     +                      LOCAL )
     +
     +
C
C
C             ...determine the symmetry partners. A maximum of MXITER
C                attempts is beeing done, each time relaxing more and
C                more the symmetry accuracies involved. After each
C                try the QTHRESH and the LTHRESH settings are increased
C                by a factor QLSTEP with an appropriate message. If
C                all attempts do not produce the same amount of partners
C                for all symmetry related centers, the program stops
C                with a message.
C
C                Once a stable # of symmetry partners has been found,
C                the initial NROT rotation indices are expanded to
C                NSYM * NROT rotation indices during this process in
C                the following sequence, starting from the last of
C                NROT original rotation indices as shown in the
C                example below for NROT = 3 and NSYM = 4. The numbers
C                in the symmetry expanded column indicate the order
C                of evaluation, which is done such that the original
C                ROT indices are not overwritten until the very end.
C
C
C                      original ROT          symmetry expanded ROT
C
C
C                                                     1
C                                                     2
C                                                     3
C
C                                                     4
C                                                     5
C                                                     6
C
C                                                     7
C                                                     8
C                                                     9
C
C                          A                         10
C                          B                         11
C                          C                         12
C
C
C
             NITER = 0
             QTHRESH = QSYMACC
             LTHRESH = LSYMACC

 1000        NITER = NITER + 1

             IF (NITER.GT.MXITER) THEN
                 WRITE (*,*) ' Maximum # of iterations reached! '
                 WRITE (*,*) ' NITER = ',NITER
                 WRITE (*,*) ' Couldnt find set of NBO rot indices! '
                 WRITE (*,*) ' nlo__find_rotation_nbo_indices '
                 WRITE (1,*) ' Maximum # of iterations reached! '
                 WRITE (1,*) ' NITER = ',NITER
                 WRITE (1,*) ' Couldnt find set of NBO rot indices! '
                 WRITE (1,*) ' nlo__find_rotation_nbo_indices '
                 STOP
             END IF

             DO ROT = NROT,1,-1
                NBOBEG = ROTIDX (ROT)
                NBOEND = NBOBEG + NXBA - 1
                NBOMAX = NBOBEG + OFF

                NBO1ST = NBOMAX
                NBOLST = NBOMAX - 1

                DO NBO = NBOMAX,NBOBEG,-1
                   QNBO = Q     (NBO)
                   LNBO = LOCAL (NBO)
                   QSYM = DABS (QMAX - QNBO) .LT. QTHRESH
                   LSYM = DABS (LMAX - LNBO) .LT. LTHRESH
                   IF (QSYM.AND.LSYM) THEN
                       NBO1ST = NBO
                   END IF
                END DO

                DO NBO = NBOMAX,NBOEND
                   QNBO = Q     (NBO)
                   LNBO = LOCAL (NBO)
                   QSYM = DABS (QMAX - QNBO) .LT. QTHRESH
                   LSYM = DABS (LMAX - LNBO) .LT. LTHRESH
                   IF (QSYM.AND.LSYM) THEN
                       NBOLST = NBO
                   END IF
                END DO

                NSYM = NBOLST - NBO1ST + 1

                IF (NSYM.LT.1) THEN
                    WRITE (*,*) ' No Rot index defined! '
                    WRITE (*,*) ' NSYM = ',NSYM
                    WRITE (*,*) ' nlo__find_rotation_nbo_indices '
                    WRITE (1,*) ' No Rot index defined! '
                    WRITE (1,*) ' NSYM = ',NSYM
                    WRITE (1,*) ' nlo__find_rotation_nbo_indices '
                    STOP
                END IF

                IF (NSYM.GT.3) THEN
                    WRITE (*,*) ' Warning! # Rot indices > 3! '
                    WRITE (*,*) ' NSYM = ',NSYM
                    WRITE (*,*) ' nlo__find_rotation_nbo_indices '
                    WRITE (1,*) ' Warning! # Rot indices > 3! '
                    WRITE (1,*) ' NSYM = ',NSYM
                    WRITE (1,*) ' nlo__find_rotation_nbo_indices '
                END IF

                IF (ROT.EQ.NROT) THEN
                    NSYMOLD = NSYM
                ELSE
                    IF (NSYM.NE.NSYMOLD) THEN
                        QTHRESH = QTHRESH + QSYMACC
                        LTHRESH = LTHRESH + LSYMACC
                        WRITE (*,*) ' Bad # of Rot indices! '
                        WRITE (*,*) ' NSYM,NSYMOLD = ',NSYM,NSYMOLD
                        WRITE (*,*) ' Increasing thresholds ... '
                        WRITE (*,*) ' Q,L-thresh = ',QTHRESH,LTHRESH
                        WRITE (1,*) ' Bad # of Rot indices! '
                        WRITE (1,*) ' NSYM,NSYMOLD = ',NSYM,NSYMOLD
                        WRITE (1,*) ' Increasing thresholds ... '
                        WRITE (1,*) ' Q,L-thresh = ',QTHRESH,LTHRESH
                        GOTO 1000
                    END IF
                END IF
             END DO
C
C
C             ...the # of symmetry partners NSYM is stable. Run the
C                symmetry partner determination loop one more time
C                using the NSYM and final QTHRESH and LTHRESH values
C                and store all partner rotation indices.
C
C
             MROT = NSYM * NROT + 1

             DO ROT = NROT,1,-1
                NBOBEG = ROTIDX (ROT)
                NBOEND = NBOBEG + NXBA - 1
                NBOMAX = NBOBEG + OFF
                NBO1ST = NBOMAX
                NBOLST = NBOMAX - 1
                DO NBO = NBOMAX,NBOBEG,-1
                   QNBO = Q     (NBO)
                   LNBO = LOCAL (NBO)
                   QSYM = DABS (QMAX - QNBO) .LT. QTHRESH
                   LSYM = DABS (LMAX - LNBO) .LT. LTHRESH
                   IF (QSYM.AND.LSYM) THEN
                       NBO1ST = NBO
                   END IF
                END DO
                DO NBO = NBOMAX,NBOEND
                   QNBO = Q     (NBO)
                   LNBO = LOCAL (NBO)
                   QSYM = DABS (QMAX - QNBO) .LT. QTHRESH
                   LSYM = DABS (LMAX - LNBO) .LT. LTHRESH
                   IF (QSYM.AND.LSYM) THEN
                       NBOLST = NBO
                   END IF
                END DO
                DO NBO = NBOLST,NBO1ST,-1
                   MROT = MROT - 1
                   ROTIDX (MROT) = NBO
                END DO
             END DO

             NROT = NSYM * NROT

         ELSE
C
C
C             ...the non-atomic 2-center bonds and antibonds case.
C                Again the code provides for the option of performing
C                several iterations with respect to symmetry accuracy.
C
C
             NITER = 0
             QTHRESH = QSYMACC
             LTHRESH = LSYMACC

 2000        NITER = NITER + 1

             IF (NITER.GT.MXITER) THEN
                 WRITE (*,*) ' Maximum # of iterations reached! '
                 WRITE (*,*) ' NITER = ',NITER
                 WRITE (*,*) ' Couldnt find set of NBO rot indices! '
                 WRITE (*,*) ' nlo__find_rotation_nbo_indices '
                 WRITE (1,*) ' Maximum # of iterations reached! '
                 WRITE (1,*) ' NITER = ',NITER
                 WRITE (1,*) ' Couldnt find set of NBO rot indices! '
                 WRITE (1,*) ' nlo__find_rotation_nbo_indices '
                 STOP
             END IF

             NROT = 0

             DO NBO = NBO1ST,NBOLST

                BOND = NBOBD (NBO)
                NCEN = BDNCEN (BOND)

                IF (NCEN.EQ.2) THEN

                    ATX  = 0
                    ATX1 = BDCEN (1,BOND)
                    ATX2 = BDCEN (2,BOND)

                    IF (ATX1.EQ.CENTER) THEN
                        ATX = ATX2
                    ELSE IF (ATX2.EQ.CENTER) THEN
                        ATX = ATX1
                    END IF

                    IF (ATX.NE.0) THEN
                        DO SYM = SYM1ST,SYMLST
                           IF (ATX.EQ.SYMCEN (SYM)) THEN
                               QNBO = Q     (NBO)
                               LNBO = LOCAL (NBO)
                               QSYM = DABS (QMAX - QNBO).LT.QTHRESH
                               LSYM = DABS (LMAX - LNBO).LT.LTHRESH
                               IF (QSYM.AND.LSYM) THEN
                                   NROT = NROT + 1
                                   ROTIDX (NROT) = NBO
                               END IF
                           END IF
                        END DO
                    END IF

                END IF
             END DO

             IF (NROT.NE.SYMSIZE) THEN
                 QTHRESH = QTHRESH + QSYMACC
                 LTHRESH = LTHRESH + LSYMACC
                 WRITE (*,*) ' Bad # of Rot indices! '
                 WRITE (*,*) ' SYMSIZE,NROT = ',SYMSIZE,NROT
                 WRITE (*,*) ' Increasing thresholds ... '
                 WRITE (*,*) ' Q,L-thresh = ',QTHRESH,LTHRESH
                 WRITE (1,*) ' Bad # of Rot indices! '
                 WRITE (1,*) ' SYMSIZE,NROT = ',SYMSIZE,NROT
                 WRITE (1,*) ' Increasing thresholds ... '
                 WRITE (1,*) ' Q,L-thresh = ',QTHRESH,LTHRESH
                 GOTO 2000
             END IF

         END IF
C
C
C             ...printout of rotation indices for checking.
C
C
         DO ROT = 1,NROT
            WRITE (*,*) ' Rot index = ',ROTIDX (ROT)
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
