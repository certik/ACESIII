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
         SUBROUTINE  NLO__ANALYZE_NAO_ORBITALS
     +
     +                    ( NBAS,NATOM,
     +                      MXSHELL,MXNAL,
     +                      MAXOCC,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      ZATOM,
     +                      NMB,
     +                      LSIZE,
     +                      BASBEG,BASEND,
     +                      PSYMACC,
     +                      WNAOVAL,WNAORYD,
     +                      WPRE,W,
     +                      XVEC,
     +
     +                              NCB,NVB,NRB,
     +                              WRYDAT,
     +                              POPCOR,POPVAL,POPRYD,POPSUM,
     +                              NAOTYP,
     +                              COLMAP,SYMMAP,
     +                              CSHELL,VSHELL,RSHELL,
     +                              NCA,
     +                              MXNCBA,
     +                              ATNCB,ATCIDX,ATCOFF,
     +                              NVA,
     +                              MXNVBA,
     +                              ATNVB,ATVIDX,ATVOFF,
     +                              NRA,
     +                              MXNRBA,
     +                              ATNRB,ATRIDX,ATROFF,
     +                              CYCLERYD,CYCLESYM )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ANALYZE_NAO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine analyzes the NAO orbitals obtained by
C                looking at their weights and decomposing the NAO set
C                into Core, Valence and Rydberg NAO orbitals. The
C                criteria for this decomposition is governed by two
C                NAO weight thresholds:
C
C                  WNAOVAL  =<  weight  =<   MAXOCC     Core
C                  WNAORYD  =<  weight   <   WNAOVAL    Valence
C                        0  =<  weight   <   WNAORYD    Rydberg
C
C                Once the decomposition is established, the total
C                # of Rydberg NAO's is compared with the total # of
C                Rydberg NAO's obtained previously. If there is any
C                difference in these two numbers, then we know that
C                the original pre-NAO decomposition into minimal and
C                Rydberg sets was probably based on a wrong Rydberg
C                weight threshold. In this case we exit the present
C                routine with the appropriate pre-NAO Rydberg weight
C                threshold WRYDBG, which is based on the old pre-NAO
C                weights, and set the Rydberg cycle indicator true.
C
C                This routine also determines, if atomic centers are
C                related by symmetry or not, which is based on equal
C                # of shells, shell types and equal weight summations
C                per shell type. Any change in the atomic center
C                symmetry pattern sets the symmetry cycle indicator
C                true. Thus the NAO generation has been successful
C                only if the Rydberg space and the symmetry pattern
C                are both stable.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MXNAL        =  maximum size of atomic l-shell
C                                    space. The atomic l-shell space
C                                    is the total # of contractions for
C                                    an atomic l-shell.
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A (in increasing order!).
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    NMB          =  # of minimal NAO's.
C                    LSIZE (I)    =  I-th l-shell size
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
C                    PSYMACC      =  symmetry accuracy for occupation
C                                    matrix values (including weights)
C                    WNAOVAL      =  the Valence NAO weight threshold.
C                    WNAORYD      =  the Rydberg NAO weight threshold.
C                    WPRE         =  pre-NAO weight vector.
C                    W            =  NAO weight vector.
C                    XVEC         =  flp scratch array of vector type.
C                    NRB          =  original # of Rydberg NAO's.
C                    SYMMAP (A,B) =  original symmetry map matrix,
C                                    containing a 1, if atoms A and B
C                                    are symmetry related and a value
C                                    of 0 otherwise, at the present
C                                    stage of calculation. If this
C                                    matrix changes here, we have to
C                                    rerun the pre-NAO/NAO generation
C                                    sequence, hence CYCLE will be set
C                                    true.
C
C
C                  Output:
C
C                    NCB,NVB,NRB  =  # of Core, Valence and Rydberg
C                                    NAOs found.
C                    WRYDAT (A)   =  the new pre-NAO weight threshold
C                                    for atom A (if necessary).
C                    POPxxx (A)   =  # of electrons populating the
C                                    Core, Valence, Rydberg and all
C                                    NAOs on atom A (xxx = COR,VAL,
C                                    RYD and SUM, respectively).
C                    NAOTYP (I)   =  integer 2,1 or 0, identifying
C                                    the type of NAO encountered for
C                                    the I-th NAO (= 2 for Core, = 1
C                                    for Valence and = 0 for Rydberg).
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of the I-th NAO
C                                    (based on atomic ordering) in an
C                                    NVB/NCB/NRB ordering, in which the
C                                    valence NAOs are bundled together
C                                    for further processing.
C                    SYMMAP (A,B) =  is set equal to 1 if atoms A and
C                                    B are symmetry related as judged
C                                    by comparison of the NAO weight
C                                    sums of the individual shells.
C                                    If no symmetry related the value
C                                    given is 0. Note, that if the
C                                    atomic Z values differ or the
C                                    overall NAO dimension is different
C                                    for both atoms then they cannot be
C                                    symmetry related.
C                    CSHELL (I)   =  l-shell type for the I-th Core NAO.
C                    VSHELL (I)   =  l-shell type for the I-th Valence
C                                    NAO.
C                    RSHELL (I)   =  l-shell type for the I-th Rydberg
C                                    NAO.
C                    NCA          =  # of atoms on which Core NAOs have
C                                    been found.
C                    MXNCBA       =  maximum # of Core NAOs per atom.
C                    ATNCB (A)    =  # of Core NAOs on core atom A.
C                    ATCIDX (A)   =  atomic index for core atom A.
C                    ATCOFF (A)   =  index offset for Core NAOs for
C                                    core atom A. This index is equal
C                                    to the total number of Core NAOs
C                                    on all core atoms preceeding core
C                                    atom A.
C                    NVA          =  # of atoms on which Valence NAOs
C                                    have been found.
C                    MXNVBA       =  maximum # of Valence NAOs per
C                                    atom.
C                    ATNVB (A)    =  # of Valence NAOs on valence
C                                    atom A
C                    ATVIDX (A)   =  atomic index for valence atom A.
C                    ATVOFF (A)   =  index offset for Valence NAOs
C                                    for valence atom A. This index is
C                                    equal to the total number of
C                                    valence NAO's on all valence atoms
C                                    preceeding valence atom A.
C                    NRA          =  # of atoms on which Rydberg NAOs
C                                    have been found.
C                    MXNRBA       =  maximum # of Rydberg NAOs per atom.
C                    ATNRB (A)    =  # of Rydberg NAOs on core atom A.
C                    ATRIDX (A)   =  atomic index for Rydberg atom A.
C                    ATROFF (A)   =  index offset for Rydberg NAOs for
C                                    Rydberg atom A. This index is equal
C                                    to the total number of Rydberg NAOs
C                                    on all Rydberg atoms preceeding
C                                    Rydberg atom A.
C                    CYCLERYD     =  cycle indicator for rerunning the
C                                    pre-NAO/NAO generation sequence
C                                    in case a Rydberg space dimension
C                                    mismatch was found. If true, the
C                                    sequence will be rerun with new
C                                    pre-NAO Rydberg atomic limits.
C                    CYCLESYM     =  cycle indicator for rerunning the
C                                    pre-NAO/NAO generation sequence
C                                    in case higher symmetry was found.
C                                    If true, the sequence will be
C                                    rerun.
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

         LOGICAL     CYCLERYD,CYCLESYM
         LOGICAL     EQUALAB
         LOGICAL     ZVALEQ,NBASEQ,NLTYPEQ,NALEQ,LTYPEQ,LDIMEQ

         INTEGER     A,B,I,L,M,N
         INTEGER     ATOM
         INTEGER     BASOFFA,BASOFFB
         INTEGER     BASNR,BASNRA,BASNRB
         INTEGER     INDEX
         INTEGER     LDIM,LTOT,LTYPE
         INTEGER     MXSHELL,MXNAL,MXNCBA,MXNVBA,MXNRBA
         INTEGER     NBAS,NBAS2,NBASA,NBASB
         INTEGER     NAL
         INTEGER     NATOM
         INTEGER     NLTYPE
         INTEGER     NCA,NVA,NRA
         INTEGER     NCB,NVB,NRB,NMB
         INTEGER     NCBA,NVBA,NRBA
         INTEGER     NCBAOLD,NVBAOLD,NRBAOLD
         INTEGER     NRBPREV
         INTEGER     ZATOMA,ZATOMB

         INTEGER     ATNCB   (1:NATOM)
         INTEGER     ATNVB   (1:NATOM)
         INTEGER     ATNRB   (1:NATOM)
         INTEGER     ATCIDX  (1:NATOM)
         INTEGER     ATVIDX  (1:NATOM)
         INTEGER     ATRIDX  (1:NATOM)
         INTEGER     ATCOFF  (1:NATOM)
         INTEGER     ATVOFF  (1:NATOM)
         INTEGER     ATROFF  (1:NATOM)
         INTEGER     BASBEG  (1:NATOM)
         INTEGER     BASEND  (1:NATOM)
         INTEGER     COLMAP  (1:NBAS)
         INTEGER     CSHELL  (1:NBAS)
         INTEGER     VSHELL  (1:NBAS)
         INTEGER     RSHELL  (1:NBAS)
         INTEGER     LSIZE   (0:MXSHELL)
         INTEGER     NAOTYP  (1:NBAS)
         INTEGER     NSHELLS (1:NATOM)
         INTEGER     ZATOM   (1:NATOM)

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM    ,1:NATOM)

         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  PSYMACC
         DOUBLE PRECISION  THRESH
         DOUBLE PRECISION  WEIGHT,WNAORYD,WNAOVAL,WSUMA,WSUMB,WDIFF
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  POPCOR  (1:NATOM)
         DOUBLE PRECISION  POPVAL  (1:NATOM)
         DOUBLE PRECISION  POPRYD  (1:NATOM)
         DOUBLE PRECISION  POPSUM  (1:NATOM)
         DOUBLE PRECISION  W       (1:NBAS )
         DOUBLE PRECISION  WPRE    (1:NBAS )
         DOUBLE PRECISION  WRYDAT  (1:NATOM)
         DOUBLE PRECISION  XVEC    (1:NATOM)

         PARAMETER    (ZERO   =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...set initial values. Zero the atomic populations and
C                the temporary array which will hold the new atomic
C                Rydberg pre-NAO weight bounds.
C
C
         CYCLERYD = .FALSE.
         CYCLESYM = .FALSE.

         DO ATOM = 1,NATOM
            POPCOR (ATOM) = ZERO
            POPVAL (ATOM) = ZERO
            POPRYD (ATOM) = ZERO
            POPSUM (ATOM) = ZERO
            XVEC   (ATOM) = ZERO
         END DO
C
C
C             ...analyze the NAO weights and decompose the NAO set
C                into Core, Valence and Rydberg. Find here also
C                the new atomic Rydberg pre-NAO weight bounds, even
C                if they will not be used later, and store them in
C                a temporary array.
C
C
         NBAS2 = NBAS + NBAS
         BASNR = 0
         NRBPREV = NRB

         NCB = 0
         NVB = 0
         NRB = 0
         NCA = 0
         NVA = 0
         NRA = 0
         MXNCBA = 0
         MXNVBA = 0
         MXNRBA = 0

         DO 1000 ATOM = 1,NATOM

            NCBA = 0
            NVBA = 0
            NRBA = 0
            NCBAOLD = NCB
            NVBAOLD = NVB
            NRBAOLD = NRB

            NLTYPE = NSHELLS (ATOM)

            DO 2000 N = 1,NLTYPE

               NAL = NBASAL (N,ATOM)
               LTYPE = SHELLS (N,ATOM)
               LDIM = LSIZE (LTYPE)
               LTOT = NAL * LDIM

               DO L = 1,NAL

                  WEIGHT = ZERO
                  DO M = 1,LDIM
                     WEIGHT = WEIGHT + W (BASNR+M)
                  END DO
                  WEIGHT = WEIGHT / DFLOAT (LDIM)

                  IF (WEIGHT.GT.(MAXOCC+PSYMACC)) THEN
                      WRITE (*,*) ' shell NAO weight > ',MAXOCC+PSYMACC
                      WRITE (*,*) ' WEIGHT = ',WEIGHT
                      WRITE (*,*) ' nlo__analyze_nao_orbitals '
                      WRITE (1,*) ' shell NAO weight > ',MAXOCC+PSYMACC
                      WRITE (1,*) ' WEIGHT = ',WEIGHT
                      WRITE (1,*) ' nlo__analyze_nao_orbitals '
                      STOP
                  ELSE IF (WEIGHT.GE.WNAOVAL) THEN
                      DO M = 1,LDIM
                         NCB = NCB + 1
                         NCBA = NCBA + 1
                         COLMAP (BASNR+M) = NCB + NBAS
                         NAOTYP (BASNR+M) = 2
                         CSHELL (NCB) = LTYPE
                         POPCOR (ATOM) = POPCOR (ATOM) + W (BASNR+M)
                         POPSUM (ATOM) = POPSUM (ATOM) + W (BASNR+M)
                      END DO
                  ELSE IF (WEIGHT.GE.WNAORYD) THEN
                      DO M = 1,LDIM
                         NVB = NVB + 1
                         NVBA = NVBA + 1
                         COLMAP (BASNR+M) = NVB
                         NAOTYP (BASNR+M) = 1
                         VSHELL (NVB) = LTYPE
                         POPVAL (ATOM) = POPVAL (ATOM) + W (BASNR+M)
                         POPSUM (ATOM) = POPSUM (ATOM) + W (BASNR+M)
                      END DO
                  ELSE IF (WEIGHT.GE.ZERO) THEN
                      DO M = 1,LDIM
                         NRB = NRB + 1
                         NRBA = NRBA + 1
                         COLMAP (BASNR+M) = NRB + NBAS2
                         NAOTYP (BASNR+M) = 0
                         RSHELL (NRB) = LTYPE
                         POPRYD (ATOM) = POPRYD (ATOM) + W (BASNR+M)
                         POPSUM (ATOM) = POPSUM (ATOM) + W (BASNR+M)
                      END DO
                      M = BASNR + 1
                      XVEC (ATOM) = MAX (XVEC (ATOM),WPRE (M))
                  ELSE IF (WEIGHT.GE.-PSYMACC) THEN

                      WRITE (*,*) ' WARNING! Negative NAO weight! '
                      WRITE (*,*) -PSYMACC,' < shell NAO weight < 0 '
                      WRITE (*,*) ' WEIGHT = ',WEIGHT
                      WRITE (*,*) ' nlo__analyze_nao_orbitals '
                      WRITE (*,*) ' Proceeding ... '
                      WRITE (1,*) ' WARNING! Negative NAO weight! '
                      WRITE (1,*) -PSYMACC,' < shell NAO weight < 0 '
                      WRITE (1,*) ' WEIGHT = ',WEIGHT
                      WRITE (1,*) ' nlo__analyze_nao_orbitals '
                      WRITE (1,*) ' Proceeding ... '

                      DO M = 1,LDIM
                         NRB = NRB + 1
                         NRBA = NRBA + 1
                         COLMAP (BASNR+M) = NRB + NBAS2
                         NAOTYP (BASNR+M) = 0
                         RSHELL (NRB) = LTYPE
                         POPRYD (ATOM) = POPRYD (ATOM) + W (BASNR+M)
                         POPSUM (ATOM) = POPSUM (ATOM) + W (BASNR+M)
                      END DO
                      M = BASNR + 1
                      XVEC (ATOM) = MAX (XVEC (ATOM),WPRE (M))
                  ELSE
                      WRITE (*,*) ' shell NAO weight < ',-PSYMACC
                      WRITE (*,*) ' WEIGHT = ',WEIGHT
                      WRITE (*,*) ' nlo__analyze_nao_orbitals '
                      WRITE (1,*) ' max shell NAO weight < ',-PSYMACC
                      WRITE (1,*) ' WEIGHT = ',WEIGHT
                      WRITE (1,*) ' nlo__analyze_nao_orbitals '
                      STOP
                  END IF

                  BASNR = BASNR + LDIM

               END DO
 2000       CONTINUE

            IF (NCBA.GT.0) THEN
                NCA = NCA + 1
                ATNCB (NCA) = NCBA
                ATCIDX (NCA) = ATOM
                ATCOFF (NCA) = NCBAOLD
                MXNCBA = MAX0 (NCBA,MXNCBA)
            END IF

            IF (NVBA.GT.0) THEN
                NVA = NVA + 1
                ATNVB (NVA) = NVBA
                ATVIDX (NVA) = ATOM
                ATVOFF (NVA) = NVBAOLD
                MXNVBA = MAX0 (NVBA,MXNVBA)
            END IF

            IF (NRBA.GT.0) THEN
                NRA = NRA + 1
                ATNRB (NRA) = NRBA
                ATRIDX (NRA) = ATOM
                ATROFF (NRA) = NRBAOLD
                MXNRBA = MAX0 (NRBA,MXNRBA)
            END IF

 1000    CONTINUE
C
C
C             ...check, if the # of Rydberg NAO's found conforms
C                with the # of Rydberg pre-NAO's found previously.
C                If not, transfer the new atomic Rydberg pre-NAO
C                weight bounds from the temporary array and exit
C                routine with Rydberg cycle indicator set true.
C
C
         CYCLERYD = NRB .NE. NRBPREV

         IF (CYCLERYD) THEN
             WRITE (*,*) ' Rydberg CYCLE! '
             DO ATOM = 1,NATOM
                WRYDAT (ATOM) = XVEC (ATOM)
             END DO
             RETURN
         END IF
C
C
C             ...check, if the # of Core and Valence NAOs add up ok
C                to the # of minimal pre-NAOs.
C
C
         IF ((NCB+NVB).NE.NMB) THEN
             WRITE (*,*) ' Problems in finding correct # of NAOs! '
             WRITE (*,*) ' NCB,NVB,NMB = ',NCB,NVB,NMB
             WRITE (*,*) ' nlo__analyze_nao_orbitals '
             WRITE (1,*) ' Problems in finding correct # of NAOs! '
             WRITE (1,*) ' NCB,NVB,NMB = ',NCB,NVB,NMB
             WRITE (1,*) ' nlo__analyze_nao_orbitals '
             STOP
         END IF
C
C
C             ...no more NAO determination cycles necessary. Update
C                the NCB and NRB sections of the column map.
C
C
         DO I = 1,NBAS
            INDEX = COLMAP (I)
            IF (INDEX.GT.NBAS2) THEN
                INDEX = INDEX - NBAS2 + NVB + NCB
                COLMAP (I) = INDEX
            ELSE IF (INDEX.GT.NBAS) THEN
                INDEX = INDEX - NBAS + NVB
                COLMAP (I) = INDEX
            END IF
         END DO
C
C
C             ...determine the new symmetry map by looping over all
C                unique atomic pairs and compare with the old symmetry
C                map values. Set the CYCLESYM indicator = .true., if
C                any change occurs.
C
C
         DO A = 1,NATOM
            ZATOMA  = ZATOM (A)
            BASOFFA = BASBEG (A) - 1
            NBASA   = BASEND (A) - BASOFFA
            NLTYPE  = NSHELLS (A)

            DO B = A+1,NATOM
               ZATOMB  = ZATOM (B)
               BASOFFB = BASBEG (B) - 1
               NBASB   = BASEND (B) - BASOFFB

               ZVALEQ  = ZATOMA .EQ. ZATOMB
               NBASEQ  = NBASA  .EQ. NBASB
               NLTYPEQ = NLTYPE .EQ. NSHELLS (B)

               EQUALAB = ZVALEQ.AND.NBASEQ.AND.NLTYPEQ

               IF (EQUALAB) THEN

                   BASNRA = BASOFFA
                   BASNRB = BASOFFB

                   DO N = 1,NLTYPE
                      NAL   = NBASAL (N,A)
                      LTYPE = SHELLS (N,A)
                      LDIM  = LSIZE (LTYPE)

                      NALEQ  = NAL   .EQ. NBASAL (N,B)
                      LTYPEQ = LTYPE .EQ. SHELLS (N,B)
                      LDIMEQ = LDIM  .EQ. LSIZE (SHELLS (N,B))

                      EQUALAB = EQUALAB.AND.NALEQ.AND.LTYPEQ.AND.LDIMEQ

                      IF (EQUALAB) THEN
                          DO L = 1,NAL
                             WSUMA = ZERO
                             WSUMB = ZERO
                             DO M = 1,LDIM
                                WSUMA = WSUMA + W (BASNRA+M)
                                WSUMB = WSUMB + W (BASNRB+M)
                             END DO
                             WDIFF = DABS (WSUMA - WSUMB)
                             THRESH = PSYMACC * DFLOAT (LDIM)
                             EQUALAB = EQUALAB .AND. (WDIFF.LT.THRESH)
                             BASNRA = BASNRA + LDIM
                             BASNRB = BASNRB + LDIM
                          END DO
                      END IF

                   END DO
               END IF

               IF (EQUALAB) THEN
                   CYCLESYM = SYMMAP (A,B) .NE. 1
                   SYMMAP (A,B) = 1
                   SYMMAP (B,A) = 1
               END IF

            END DO
         END DO

         IF (CYCLESYM) THEN
             WRITE (*,*) ' Symmetry CYCLE! '
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
