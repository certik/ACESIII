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
         SUBROUTINE  NLO__PRINT_NHO_RESULTS
     +
     +                    ( UNITID,
     +                      NBAS,NATOM,NBOND,
     +                      NBOSIZE,
     +                      MXNBA,MXSHELL,
     +                      NHB,NCB,NRB,
     +                      NHA,NCA,NRA,
     +                      ZATOM,
     +                      RYD2HYB,
     +                      ATNHB,ATHVAL,ATHIDX,ATHOFF,
     +                      ATNCB,ATCIDX,ATCOFF,
     +                      ATNRB,ATRIDX,ATROFF,
     +                      HSHELL,CSHELL,RSHELL,
     +                      BDNCEN,BDCEN,BDNBAS,BDBAS,BDOCC,
     +                      COLMAP,NAOIDX,
     +                      ANGNHO,
     +                      H )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NHO_RESULTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints details of the final set of
C                natural hybrid orbitals (NHOs) to be used to establish
C                the natural bond orbitals (NBOs) to a file specified
C                by its unit identification number.
C
C                For each atomic center the following is printed:
C
C                   1) Expansion of NHOs in terms of the valence NAOs.
C
C                   2) NHO hybridization pattern.
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    NBOND        =  total # of bonds
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MXNBA        =  overall maximum between the
C                                    maximum # of Hybrid, Core and
C                                    Rydberg NAOs per atom.
C                    MXSHELL      =  largest l-shell value
C                    NxB          =  total # of Hybrid, Core and
C                                    Rydberg NAOs (x=H,C,R)
C                    NxA          =  total # of hybrid, core and
C                                    Rydberg atoms (x=H,C,R)
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    has been included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    ATHIDX (A)   =  atomic index for hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid
C                                    NAOs on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    ATNCB (A)    =  # of Core NAOs on core atom A.
C                    ATCIDX (A)   =  atomic index for core atom A.
C                    ATCOFF (A)   =  index offset for Core NAOs for
C                                    core atom A. This index is equal
C                                    to the total number of Core NAOs
C                                    on all core atoms preceeding
C                                    core atom A.
C                    ATNRB (A)    =  # of Rydberg NAOs on Rydberg atom A.
C                    ATRIDX (A)   =  atomic index for Rydberg atom A.
C                    ATROFF (A)   =  index offset for Rydberg NAOs
C                                    for Rydberg atom A. This index
C                                    is equal to the total number of
C                                    Rydberg NAOs on all Rydberg atoms
C                                    preceeding Rydberg atom A.
C                    HSHELL (I)   =  l-shell type for I-th Hybrid NAO.
C                    CSHELL (I)   =  l-shell type for I-th Core NAO.
C                    RSHELL (I)   =  l-shell type for I-th Rydberg NAO.
C                    BDNCEN (J)   =  # of atomic centers for J-th bond.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond.
C                    BDNBAS (J)   =  # of basis functions (NHOs) for
C                                    J-th bond.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond.
C                    BDOCC (J)    =  # of occupied levels for J-th bond.
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of I-th NAO, based
C                                    on atomic order, in the NHB/NCB/NRB
C                                    order. This array is needed for
C                                    relating the printing of the NHOs
C                                    in terms of the original NAO
C                                    indices as they came out of the
C                                    NAO generation routine (i.e. before
C                                    the atomic -> NHB/NCB/NRB ordering.
C                    NAOIDX (I)   =  will contain the inverse of the
C                                    COLMAP array. Hence NAOIDX (I) will
C                                    contain the original NAO index
C                                    in the atomic order of the I-th
C                                    NHB/NCB/NRB reordered NAO.
C                    ANGNHO (J,L) =  NBAS x (MXSHELL+1) matrix
C                                    containing the sum of the square
C                                    of the NAO coefficients per L-type
C                                    angular momentum for the J-th NHO.
C                                    The NHO order is NHB/NCB/NRB.
C                    H (I,J)      =  MXNBA x NBAS matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB/NCB/NRB
C                                    NHOs in that order, with all NHOs
C                                    belonging to a specific atomic
C                                    center being grouped together.
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

         CHARACTER*2   ATCHAR
         CHARACTER*1   BAR,DASH

         CHARACTER*2   ATSYMB  (1:104)
         CHARACTER*1   HYBRID  (1:50)
         CHARACTER*3   LSYMB   (0:7)
         CHARACTER*7   NAOSYMB (0:2)

         LOGICAL     RYD2HYB

         INTEGER     ATIDX
         INTEGER     ATOM,BOND
         INTEGER     I,J,L
         INTEGER     LENGTH
         INTEGER     LTYPE
         INTEGER     MXNBA,MXSHELL
         INTEGER     NAO,NHO
         INTEGER     NATOM,NBOND
         INTEGER     NBAS
         INTEGER     NBOSIZE
         INTEGER     NOCC,NCEN,NNHO
         INTEGER     NHB,NCB,NRB
         INTEGER     NHA,NCA,NRA
         INTEGER     NHBA,NVBA,NCBA,NRBA
         INTEGER     NSHELL
         INTEGER     OFF
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     ATNHB  (1:NHA  )
         INTEGER     ATNCB  (1:NCA  )
         INTEGER     ATNRB  (1:NRA  )
         INTEGER     ATHIDX (1:NHA  )
         INTEGER     ATCIDX (1:NCA  )
         INTEGER     ATRIDX (1:NRA  )
         INTEGER     ATHOFF (1:NHA  )
         INTEGER     ATCOFF (1:NCA  )
         INTEGER     ATROFF (1:NRA  )
         INTEGER     ATHVAL (1:NHA  )
         INTEGER     BDNBAS (1:NHB  )
         INTEGER     BDNCEN (1:NHB  )
         INTEGER     BDOCC  (1:NHB  )
         INTEGER     COLMAP (1:NBAS )
         INTEGER     NAOIDX (1:NBAS )
         INTEGER     HSHELL (1:NHB  )
         INTEGER     CSHELL (1:NCB  )
         INTEGER     RSHELL (1:NRB  )
         INTEGER     ZATOM  (1:NATOM)

         INTEGER     BDBAS   (1:NBOSIZE,1:NHB)
         INTEGER     BDCEN   (1:NBOSIZE,1:NHB)

         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  W (0:7)

         DOUBLE PRECISION  ANGNHO (1:NBAS,0:MXSHELL)

         DOUBLE PRECISION  H (1:MXNBA,1:NBAS)

         DATA ATSYMB /' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     +                'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     +                'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     +                'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     +                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     +                'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     +                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +                'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     +                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     +                'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     +                'Md','No','Lr','xx'/
         DATA  BAR     /'|'/
         DATA  DASH    /'-'/
         DATA  LSYMB   /' s ',' p ',' d ',' f ',' g ',' h ',' i ','> i'/
         DATA  NAOSYMB /'Rydberg','Valence','   Core'/

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...before starting the print we generate the inverse
C                NAOIDX of the COLMAP array. Also we determine the
C                number of shells that will be considered for printing.
C
C
         DO I = 1,NBAS
            NAOIDX (COLMAP (I)) = I
         END DO

         NSHELL = MIN (MXSHELL,7)
C
C
C             ...print out title for valence part.
C
C
         WRITE (UNITID,8100) 'Valence Natural Hybrid Orbital Analysis'
 8100    FORMAT (//,10X,A39,//)
C
C
C             ...print out header line.
C
C
         WRITE (UNITID,9000) 'NHO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'NAO #',BAR,'LM-state',BAR,'NAO-type',
     +                        BAR,'NAO coeffs',BAR,'Hybridization'
         WRITE (UNITID,9010) (DASH,I=1,83)
C
C
C             ...print out atomic valence NHO hybrid information.
C
C
         DO 1000 ATOM = 1,NHA
            OFF = ATHOFF (ATOM)
            NHBA = ATNHB (ATOM)
            NVBA = ATHVAL (ATOM)
            ATIDX = ATHIDX (ATOM)

            ZVAL = ZATOM (ATIDX)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            DO J = 1,NVBA
               NHO = OFF + J

               DO L = 0,NSHELL
                  W (L) = ANGNHO (NHO,L)
               END DO
               IF (MXSHELL.GT.7) THEN
                   DO L = 8,MXSHELL
                      W (7) = W (7) + ANGNHO (NHO,L)
                   END DO
               END IF

               CALL  NLO__FORM_NHO_HYBRID_PATTERN
     +
     +                    ( NSHELL,
     +                      W,
     +
     +                              LENGTH,
     +                              HYBRID )
     +
     +
               NAO = OFF + 1
               LTYPE = HSHELL (NAO)
               WRITE (UNITID,9020) NHO,BAR,ATCHAR,BAR,ATIDX,BAR,
     +                             NAOIDX (NAO),BAR,LSYMB (LTYPE),
     +                             BAR,NAOSYMB (1),BAR,H (1,NHO),
     +                             BAR,(HYBRID(I),I=1,LENGTH)
               DO I = 2,NVBA
                  NAO = OFF + I
                  LTYPE = HSHELL (NAO)
                  WRITE (UNITID,9030) BAR,BAR,BAR,NAOIDX (NAO),BAR,
     +                                LSYMB (LTYPE),BAR,NAOSYMB (1),
     +                                BAR,H (I,NHO),BAR
               END DO

               DO I = NVBA+1,NHBA
                  NAO = OFF + I
                  LTYPE = HSHELL (NAO)
                  WRITE (UNITID,9030) BAR,BAR,BAR,NAOIDX (NAO),BAR,
     +                                LSYMB (LTYPE),BAR,NAOSYMB (0),
     +                                BAR,H (I,NHO),BAR
               END DO

               WRITE (UNITID,9010) (DASH,I=1,83)

            END DO

 1000    CONTINUE
C
C
C             ...print out title for core part.
C
C
         WRITE (UNITID,8200) 'Core Natural Hybrid Orbital Analysis'
 8200    FORMAT (//,15X,A36,//)
C
C
C             ...print out header line.
C
C
         WRITE (UNITID,9000) 'NHO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'NAO #',BAR,'LM-state',BAR,'NAO-type',
     +                        BAR,'NAO coeffs',BAR,'Hybridization'
         WRITE (UNITID,9010) (DASH,I=1,83)
C
C
C             ...print out atomic core NHO hybrid information.
C
C
         DO 2000 ATOM = 1,NCA
            OFF = ATCOFF (ATOM)
            NCBA = ATNCB (ATOM)
            ATIDX = ATCIDX (ATOM)

            ZVAL = ZATOM (ATIDX)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            DO J = 1,NCBA
               NHO = NHB + OFF + J

               DO L = 0,NSHELL
                  W (L) = ANGNHO (NHO,L)
               END DO
               IF (MXSHELL.GT.7) THEN
                   DO L = 8,MXSHELL
                      W (7) = W (7) + ANGNHO (NHO,L)
                   END DO
               END IF

               CALL  NLO__FORM_NHO_HYBRID_PATTERN
     +
     +                    ( NSHELL,
     +                      W,
     +
     +                              LENGTH,
     +                              HYBRID )
     +
     +
               NAO = OFF + 1
               LTYPE = CSHELL (NAO)
               WRITE (UNITID,9020) NHO,BAR,ATCHAR,BAR,ATIDX,BAR,
     +                             NAOIDX (NHB+NAO),BAR,LSYMB (LTYPE),
     +                             BAR,NAOSYMB (2),BAR,H (1,NHO),
     +                             BAR,(HYBRID(I),I=1,LENGTH)
               DO I = 2,NCBA
                  NAO = OFF + I
                  LTYPE = CSHELL (NAO)
                  WRITE (UNITID,9030) BAR,BAR,BAR,NAOIDX (NHB+NAO),BAR,
     +                                LSYMB (LTYPE),BAR,NAOSYMB (2),
     +                                BAR,H (I,NHO),BAR
               END DO

               WRITE (UNITID,9010) (DASH,I=1,83)

            END DO

 2000    CONTINUE
C
C
C             ...print out title for Rydberg part.
C
C
         WRITE (UNITID,8300) 'Rydberg Natural Hybrid Orbital Analysis'
 8300    FORMAT (//,10X,A39,//)
C
C
C             ...print out header line.
C
C
         WRITE (UNITID,9000) 'NHO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'NAO #',BAR,'LM-state',BAR,'NAO-type',
     +                        BAR,'NAO coeffs',BAR,'Hybridization'
         WRITE (UNITID,9010) (DASH,I=1,83)
C
C
C             ...print out atomic Rydberg NHO hybrid information
C                sitting in the hybrid section (if any).
C
C
         IF (RYD2HYB) THEN

             DO 3000 ATOM = 1,NHA
                OFF = ATHOFF (ATOM)
                NHBA = ATNHB (ATOM)
                NVBA = ATHVAL (ATOM)
                NRBA = NHBA - NVBA
                ATIDX = ATHIDX (ATOM)

                ZVAL = ZATOM (ATIDX)

                IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                    ATCHAR = ATSYMB (104)
                ELSE
                    ATCHAR = ATSYMB (ZVAL)
                END IF

                DO J = 1,NRBA
                   NHO = OFF + NVBA + J

                   DO L = 0,NSHELL
                      W (L) = ANGNHO (NHO,L)
                   END DO
                   IF (MXSHELL.GT.7) THEN
                       DO L = 8,MXSHELL
                          W (7) = W (7) + ANGNHO (NHO,L)
                       END DO
                   END IF

                   CALL  NLO__FORM_NHO_HYBRID_PATTERN
     +
     +                        ( NSHELL,
     +                          W,
     +
     +                                  LENGTH,
     +                                  HYBRID )
     +
     +
                   NAO = OFF + 1
                   LTYPE = HSHELL (NAO)
                   WRITE (UNITID,9020) NHO,BAR,ATCHAR,BAR,ATIDX,BAR,
     +                                 NAOIDX (NAO),BAR,LSYMB (LTYPE),
     +                                 BAR,NAOSYMB (1),BAR,H (1,NHO),
     +                                 BAR,(HYBRID(I),I=1,LENGTH)
                   DO I = 2,NVBA
                      NAO = OFF + I
                      LTYPE = HSHELL (NAO)
                      WRITE (UNITID,9030) BAR,BAR,BAR,NAOIDX (NAO),
     +                                    BAR,LSYMB (LTYPE),BAR,
     +                                    NAOSYMB (1),BAR,H (I,NHO),BAR
                   END DO

                   DO I = NVBA+1,NHBA
                      NAO = OFF + I
                      LTYPE = HSHELL (NAO)
                      WRITE (UNITID,9030) BAR,BAR,BAR,NAOIDX (NAO),BAR,
     +                                    LSYMB (LTYPE),BAR,NAOSYMB (0),
     +                                    BAR,H (I,NHO),BAR
                   END DO

                   WRITE (UNITID,9010) (DASH,I=1,83)

                END DO

 3000        CONTINUE

         END IF
C
C
C             ...print out atomic Rydberg NHO hybrid information.
C
C
         DO 4000 ATOM = 1,NRA
            OFF = ATROFF (ATOM)
            NRBA = ATNRB (ATOM)
            ATIDX = ATRIDX (ATOM)

            ZVAL = ZATOM (ATIDX)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            DO J = 1,NRBA
               NHO = NHB + NCB + OFF + J

               DO L = 0,NSHELL
                  W (L) = ANGNHO (NHO,L)
               END DO
               IF (MXSHELL.GT.7) THEN
                   DO L = 8,MXSHELL
                      W (7) = W (7) + ANGNHO (NHO,L)
                   END DO
               END IF

               CALL  NLO__FORM_NHO_HYBRID_PATTERN
     +
     +                    ( NSHELL,
     +                      W,
     +
     +                              LENGTH,
     +                              HYBRID )
     +
     +
               NAO = OFF + 1
               LTYPE = RSHELL (NAO)
               WRITE (UNITID,9020) NHO,BAR,ATCHAR,BAR,ATIDX,BAR,
     +                             NAOIDX (NHB+NCB+NAO),BAR,
     +                             LSYMB (LTYPE),BAR,NAOSYMB (0),
     +                             BAR,H (1,NHO),BAR,
     +                             (HYBRID(I),I=1,LENGTH)
               DO I = 2,NRBA
                  NAO = OFF + I
                  LTYPE = RSHELL (NAO)
                  WRITE (UNITID,9030) BAR,BAR,BAR,NAOIDX (NHB+NCB+NAO),
     +                                BAR,LSYMB (LTYPE),BAR,NAOSYMB (0),
     +                                BAR,H (I,NHO),BAR
               END DO

               WRITE (UNITID,9010) (DASH,I=1,83)

            END DO

 4000    CONTINUE
C
C
C             ...formats for printing the atomic NHO information.
C
C
 9000    FORMAT (1X,A5,1X,  A1, 1X,A4,1X,  A1, 1X,A6,1X, A1,
     +           1X,A5,1X,  A1, 1X,A8,1X,  A1, 1X,A8,1X, A1,
     +           1X,A10,1X, A1, 2X,A13)
 9010    FORMAT (1X,83A1)
 9020    FORMAT (2X,I3,2X,    A1, 2X,A2,2X,  A1, 2X,I3,3X, A1,
     +           1X,I4,2X,    A1, 4X,A3,3X,  A1, 2X,A7,1X, A1,
     +           1X,F10.6,1X, A1, 2X,50(A1))
 9030    FORMAT (2X,3X,2X,    A1, 2X,2X,2X,  A1, 2X,3X,3X, A1,
     +           1X,I4,2X,    A1, 4X,A3,3X,  A1, 2X,A7,1X, A1,
     +           1X,F10.6,1X, A1)
C
C
C             ...print out the NHO angular momentum weight table.
C
C
         WRITE (UNITID,8400) 'NHO Angular Momentum Weight Table'
 8400    FORMAT (//,20X,A33)

         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows = NHOs ; Cols = L-value (1=s,2=p,...) ',
     +                  NBAS,MXSHELL+1,
     +                  NBAS,MXSHELL+1,
     +                  ANGNHO )
     +
     +
C
C
C             ...print out title.
C
C
         WRITE (UNITID,8000) 'Bond formation analysis'
 8000    FORMAT (//,20X,A23,//)
C
C
C             ...print out header line and bond formation
C                information # 1.
C
C
         WRITE (UNITID,9040) 'BOND #',BAR,'Occ #',BAR,
     +                       '# of Centers',BAR,'Atomic Centers  '
         WRITE (UNITID,9050) (DASH,I=1,62)

         DO BOND = 1,NBOND
            NOCC = BDOCC (BOND)
            NCEN = BDNCEN (BOND)
            WRITE (UNITID,9060) BOND,BAR,NOCC,BAR,NCEN,BAR,
     +                          (BDCEN (I,BOND),I=1,NCEN)
         END DO
         WRITE (UNITID,9050) (DASH,I=1,62)
C
C
C             ...print out header line and bond formation
C                information # 2.
C
C
         WRITE (UNITID,9040) 'BOND #',BAR,'Occ #',BAR,
     +                       '  # of NHOs ',BAR,'Bond NHO indices'
         WRITE (UNITID,9050) (DASH,I=1,62)

         DO BOND = 1,NBOND
            NOCC = BDOCC (BOND)
            NNHO = BDNBAS (BOND)
            WRITE (UNITID,9060) BOND,BAR,NOCC,BAR,NNHO,BAR,
     +                          (BDBAS (I,BOND),I=1,NNHO)
         END DO
         WRITE (UNITID,9050) (DASH,I=1,62)
C
C
C             ...formats for printing the bond formation information.
C
C
 9040    FORMAT (1X,A6,1X, A1, 1X,A5,1X,  A1, 1X,A12,1X, A1, 2X,A16)
 9050    FORMAT (1X,62A1)
 9060    FORMAT (2X,I4,2X, A1, 2X,I3,2X,  A1, 5X,I4,5X, A1, 2X,80(I4))
C
C
C             ...ready!
C
C
         RETURN
         END
