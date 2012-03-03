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
         SUBROUTINE  NLO__CHARACTERIZE_NHO_SPACE
     +
     +                    ( NBAS,NATOM,
     +                      NVB,NVA,
     +                      MXNVBA,
     +                      VSHELL,
     +                      ATNVB,ATVIDX,ATVOFF,
     +                      RYD2HYB,
     +                      LOCMAP,
     +                      VALIDX,RYDIDX,
     +                      IVEC,
     +
     +                              ATORD,
     +                              NHB,NCB,NRB,
     +                              COLMAP,
     +                              HSHELL,CSHELL,RSHELL,
     +                              NHA,
     +                              MXNHBA,
     +                              ATNHB,ATHVAL,ATHIDX,ATHOFF,
     +                              NCA,
     +                              MXNCBA,
     +                              ATNCB,ATCIDX,ATCOFF,
     +                              NRA,
     +                              MXNRBA,
     +                              ATNRB,ATRIDX,ATROFF,
     +                              MXNBA,MXCOL )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__CHARACTERIZE_NHO_SPACE
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine characterizes the NHO space from the
C                info created by the NAO routines. The NAO orbitals
C                are classified here as Hybrid, Core and Rydberg:
C
C                The Hybrid NAOs are defined to be those which will
C                be used for NHO construction on each atomic site.
C                There are two ways, depending on the keyword RYD2HYB:
C
C                    1) RYD2HYB = false. In this case the Hybrid
C                       NAOs consist only of the Valence NAOs.
C
C                    2) RYD2HYB = true. In this case both the Valence
C                       and the Rydberg NAOs are included in the
C                       Hybrid NAOs. Their order will be such that
C                       on each atomic site the Valence NAOs come
C                       first followed by the Rydberg NAOs.
C
C                When exiting this routine, the ordering array COLMAP
C                will thus place all NHB Hybrid NAOs together, one
C                hybrid atom after the other and according to the
C                two ways 1) and 2) described above we would have
C                then the following NHB Hybrid NAO order:
C
C                 1) | Val at#1 | Val at#2 | Val at#3 | ...
C
C                 2) | Val at#1 | Ryd at#1 | Val at#2 | Ryd at#2 | ...
C
C
C                The Core NAOs remain unchanged, that is they are the
C                same as they come out from the NAO procedure.
C
C                The Rydberg NAOs are those that remain unused for
C                classifying the Hybrid NAOs. Note, that an atom
C                might have no Hybrid NAOs but still has Rydberg NAOs,
C                if it has no Valence NAOs (a rare situation, which
C                indicates an isolated atom without bonding inside the 
C                molecular framework).
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    NVB          =  # of Valence NAOs.
C                    NVA          =  # of valence atoms.
C                    MXNVBA       =  maximum # of Valence NAOs per
C                                    atom.
C                    VSHELL (I)   =  l-shell type for the I-th Valence
C                                    NAO.
C                    ATNVB (A)    =  # of Valence NAOs on valence
C                                    atom A.
C                    ATVIDX (A)   =  atomic index for valence atom A.
C                    ATVOFF (A)   =  index offset for Valence NAOs
C                                    for valence atom A. This index is
C                                    equal to the total number of
C                                    valence NAO's on all valence atoms
C                                    preceeding valence atom A.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    should be included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                                    If false, only the Valence NAO
C                                    space will be used.
C                    LOCMAP (I)   =  will contain the local map in the
C                                    active sense, such that LOCMAP (I)
C                                    contains the position index of the
C                                    I-th NAO (based on NVB/NCB/NRB
C                                    ordering) in the new NHB/NCB/NRB
C                                    ordering.
C                    VALIDX (A)   =  will contain a 'yes' or 'no' index
C                                    indicator to signal that a Rydberg
C                                    atom A has or has not Valence NAOs.
C                    RYDIDX (A)   =  will contain the Rydberg atom index
C                                    that matches valence atom A. If
C                                    valence atom A has no Rydberg NAOs
C                                    then it will be set to 0.
C                    IVEC         =  int scratch array.
C                    ATORD        =  original optimum atomic index
C                                    ordering array for all atoms.
C                    NCB          =  original # of Core NAOs.
C                    NRB          =  original # of Rydberg NAOs.
C                    COLMAP (I)   =  original column map in the active
C                                    sense, such that COLMAP (I)
C                                    contains the position index of the
C                                    I-th NAO (based on atomic ordering)
C                                    in the original NVB/NCB/NRB
C                                    ordering.
C                    CSHELL (I)   =  original l-shell type for the I-th
C                                    Core NAO.
C                    RSHELL (I)   =  original l-shell type for the I-th
C                                    Rydberg NAO.
C                    NCA          =  original # of atoms on which Core
C                                    NAOs have been found.
C                    MXNCBA       =  original maximum # of Core NAOs
C                                    per atom.
C                    ATNCB (A)    =  original # of Core NAOs on core
C                                    atom A.
C                    ATCIDX (A)   =  original atomic index for core
C                                    atom A.
C                    ATCOFF (A)   =  original index offset for Core
C                                    NAOs for core atom A. This index
C                                    is equal to the total number of
C                                    Core NAOs on all core atoms
C                                    preceeding core atom A.
C                    NRA          =  original # of atoms on which
C                                    Rydberg NAOs have been found.
C                    MXNRBA       =  original maximum # of Rydberg
C                                    NAOs per atom.
C                    ATNRB (A)    =  original # of Rydberg NAOs on
C                                    Rydberg atom A.
C                    ATRIDX (A)   =  original atomic index for Rydberg
C                                    atom A.
C                    ATROFF (A)   =  original index offset for Rydberg
C                                    NAOs for Rydberg atom A. This index
C                                    is equal to the total number of
C                                    Rydberg NAOs on all Rydberg atoms
C                                    preceeding Rydberg atom A.
C
C
C                  Output:
C
C                    ATORD        =  new optimum atomic index ordering
C                                    array for all hybrid atoms only.
C                    NHB,NCB,NRB  =  new # of Hybrid, Core and Rydberg
C                                    NAOs found.
C                    COLMAP (I)   =  new column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of the I-th NAO
C                                    (based on atomic ordering) in the
C                                    new NHB/NCB/NRB ordering.
C                    HSHELL (I)   =  l-shell type for the I-th Hybrid
C                                    NAO.
C                    CSHELL (I)   =  l-shell type for the I-th new
C                                    Core NAO.
C                    RSHELL (I)   =  l-shell type for the I-th new
C                                    Rydberg NAO.
C                    NHA          =  # of atoms on which Hybrid NAOs
C                                    have been found.
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    ATHIDX (A)   =  atomic index for hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid NAOs
C                                    on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    NCA          =  # of atoms on which new Core NAOs
C                                    have been found.
C                    MXNCBA       =  maximum # of new Core NAOs per
C                                    atom.
C                    ATNCB (A)    =  # of new Core NAOs on core atom A.
C                    ATCIDX (A)   =  atomic index for new core atom A.
C                    ATCOFF (A)   =  index offset for new Core NAOs for
C                                    new core atom A. This index is
C                                    equal to the total number of new
C                                    Core NAOs on all new core atoms
C                                    preceeding new core atom A.
C                    NRA          =  # of atoms on which new Rydberg
C                                    NAOs have been found.
C                    MXNRBA       =  maximum # of new Rydberg NAOs per
C                                    atom.
C                    ATNRB (A)    =  # of new Rydberg NAOs on new
C                                    Rydberg atom A.
C                    ATRIDX (A)   =  atomic index for new Rydberg
C                                    atom A.
C                    ATROFF (A)   =  index offset for new Rydberg NAOs
C                                    for new Rydberg atom A. This index
C                                    is equal to the total number of
C                                    new Rydberg NAOs on all new
C                                    Rydberg atoms preceeding new
C                                    Rydberg atom A.
C                    MXNBA        =  overall maximum between MXNHBA,
C                                    MXNCBA and MXNRBA (see above).
C                    MXCOL        =  maximum between NHB and MXNBA
C                                    (see above). This will denote the
C                                    maximum # of NBAS-sized columns
C                                    that are needed as flp scratch
C                                    space later on for NHO generation.
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

         LOGICAL     RYD2HYB

         INTEGER     ATOM
         INTEGER     I,M,N
         INTEGER     MBEG,NBEG
         INTEGER     MXNCBA,MXNVBA,MXNRBA,MXNHBA,MXNBA,MXCOL
         INTEGER     NBAS,NATOM
         INTEGER     NCA,NVA,NRA,NHA
         INTEGER     NCB,NVB,NRB,NHB
         INTEGER     NRAOLD
         INTEGER     NVBA,NRBA,NHBA
         INTEGER     OFFVAL,OFFRYD
         INTEGER     VALATOM,RYDATOM

         INTEGER     ATORD   (1:NATOM)
         INTEGER     ATNCB   (1:NATOM)
         INTEGER     ATNVB   (1:NATOM)
         INTEGER     ATNRB   (1:NATOM)
         INTEGER     ATNHB   (1:NATOM)
         INTEGER     ATCIDX  (1:NATOM)
         INTEGER     ATVIDX  (1:NATOM)
         INTEGER     ATRIDX  (1:NATOM)
         INTEGER     ATHIDX  (1:NATOM)
         INTEGER     ATCOFF  (1:NATOM)
         INTEGER     ATVOFF  (1:NATOM)
         INTEGER     ATROFF  (1:NATOM)
         INTEGER     ATHOFF  (1:NATOM)
         INTEGER     ATHVAL  (1:NATOM)
         INTEGER     COLMAP  (1:NBAS)
         INTEGER     CSHELL  (1:NBAS)
         INTEGER     VSHELL  (1:NBAS)
         INTEGER     RSHELL  (1:NBAS)
         INTEGER     HSHELL  (1:NBAS)

         INTEGER     IVEC    (1:NBAS)
         INTEGER     LOCMAP  (1:NBAS)
         INTEGER     RYDIDX  (1:NATOM)
         INTEGER     VALIDX  (1:NATOM)
C
C
C------------------------------------------------------------------------
C
C
C             ...update optimum atomic index ordering array to
C                include hybrid atoms only.
C
C
         NHA = NVA

         N = 0
         DO 10 I = 1,NATOM
            ATOM = ATORD (I)
            DO 20 M = 1,NVA
               IF (ATOM.EQ.ATVIDX (M)) THEN
                   N = N + 1
                   ATORD (N) = ATOM
                   GOTO 10
               END IF
   20       CONTINUE
   10    CONTINUE
C
C
C             ...handle the trivial case first (Hybrid NAOs identical
C                to Valence NAOs). In this case the Rydberg NAOs
C                do not change and also the column mapping vector
C                stays the same. In effect this is just a renaming
C                of the Valence NAO space.
C
C
         IF (.NOT.RYD2HYB) THEN

             NHB = 0
             MXNHBA = MXNVBA

             DO 100 M = 1,NVA
                NVBA = ATNVB (M)
                ATNHB (M) = NVBA
                ATHVAL (M) = NVBA
                ATHIDX (M) = ATVIDX (M)
                ATHOFF (M) = ATVOFF (M)
                DO 110 I = 1,NVBA
                   NHB = NHB + 1
                   HSHELL (NHB) = VSHELL (NHB)
  110           CONTINUE
  100        CONTINUE

             MXNBA = MAX (MXNHBA,MXNCBA,MXNRBA)

             RETURN

         END IF
C
C
C             ...the nontrivial case. Here we have to add the Rydberg
C                space (if it exists) of each valence atom on top of
C                the Valence space. Determine first, which of the
C                valence atoms has a nonzero Rydberg space and which
C                of the Rydberg atoms has a zero valence space.
C                The first info, stored in array RYDIDX as the
C                Rydberg atom index, is needed to establish the Hybrid
C                NAO space. The second info, stored in array VALIDX
C                as a simple yes (0) or no (1), is used to update the
C                new Rydberg NAO space.
C
C
         NBEG = 1
         DO 200 M = 1,NVA
            VALATOM = ATVIDX (M)
            RYDIDX (M) = 0
            DO 210 N = NBEG,NRA
               RYDATOM = ATRIDX (N)
               IF (RYDATOM.EQ.VALATOM) THEN
                   NBEG = N + 1
                   RYDIDX (M) = N
                   GOTO 200
               ELSE IF (RYDATOM.GT.VALATOM) THEN
                   GOTO 200
               END IF
  210       CONTINUE
  200    CONTINUE

         MBEG = 1
         DO 220 N = 1,NRA
            RYDATOM = ATRIDX (N)
            VALIDX (N) = 0
            DO 230 M = MBEG,NVA
               VALATOM = ATVIDX (M)
               IF (VALATOM.EQ.RYDATOM) THEN
                   MBEG = M + 1
                   VALIDX (N) = M
                   GOTO 220
               ELSE IF (VALATOM.GT.RYDATOM) THEN
                   GOTO 220
               END IF
  230       CONTINUE
  220    CONTINUE
C
C
C             ...assemble the Hybrid NAO space data and generate the
C                Hybrid NAO portion of the local map.
C
C
         NHB = 0
         MXNHBA = 0

         DO 300 M = 1,NVA
            NVBA = ATNVB (M)
            NHBA = NVBA
            ATNHB (M) = NVBA
            ATHVAL (M) = NVBA
            ATHIDX (M) = ATVIDX (M)
            ATHOFF (M) = NHB
            OFFVAL = ATVOFF (M)

            DO 310 I = 1,NVBA
               NHB = NHB + 1
               HSHELL (NHB) = VSHELL (OFFVAL+I)
               LOCMAP (OFFVAL+I) = NHB
  310       CONTINUE

            N = RYDIDX (M)
            IF (N.GT.0) THEN
                NRBA = ATNRB (N)
                NHBA = NHBA + NRBA
                OFFRYD = ATROFF (N)
                ATNHB (M) = ATNHB (M) + NRBA
                DO 320 I = 1,NRBA
                   NHB = NHB + 1
                   HSHELL (NHB) = RSHELL (OFFRYD+I)
                   LOCMAP (NVB+NCB+OFFRYD+I) = NHB
  320           CONTINUE
            END IF

            MXNHBA = MAX0 (NHBA,MXNHBA)

  300    CONTINUE
C
C
C             ...generate the Core NAO portion of the local map.
C
C
         DO 400 I = 1,NCB
            LOCMAP (NVB+I) = NHB + I
  400    CONTINUE
C
C
C             ...update the new Rydberg NAO space data and generate the
C                new Rydberg NAO portion of the local map.
C
C
         NRAOLD = NRA
         NRB = 0
         NRA = 0
         MXNRBA = 0

         DO 500 N = 1,NRAOLD
            IF (VALIDX (N).EQ.0) THEN
                NRA = NRA + 1
                NRBA = ATNRB (N)
                OFFRYD = ATROFF (N)
                ATNRB (NRA) = NRBA
                ATRIDX (NRA) = ATRIDX (N)
                ATROFF (NRA) = NRB
                DO 510 I = 1,NRBA
                   NRB = NRB + 1
                   RSHELL (NRB) = RSHELL (OFFRYD+I)
                   LOCMAP (NVB+NCB+OFFRYD+I) = NHB + NCB + NRB
  510           CONTINUE
                MXNRBA = MAX0 (NRBA,MXNRBA)
            END IF
  500    CONTINUE
C
C
C             ...combine the local map with the old column map
C                to produce the new column map.
C
C
         DO 600 I = 1,NBAS
            IVEC (I) = COLMAP (I)
  600    CONTINUE

         DO 610 I = 1,NBAS
            COLMAP (I) = LOCMAP (IVEC (I))
  610    CONTINUE
C
C
C             ...form the extra needed maximums between the # of
C                Hybrid, Core and Rydberg NAOs finally established.
C
C
         MXNBA = MAX (MXNHBA,MXNCBA,MXNRBA)
         MXCOL = MAX (NHB,MXNBA)
C
C
C             ...ready!
C
C
         RETURN
         END
