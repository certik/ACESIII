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
         SUBROUTINE  NLO__PRINT_NBO_RESULTS
     +
     +                    ( UNITID,
     +                      NBAS,NATOM,
     +                      NBOSIZE,
     +                      MXSHELL,
     +                      NCA,NLA,NEA,NYA,NRA,
     +                      NCB,NLB,NBB,NEB,NAB,NYB,NRB,
     +                      ZATOM,
     +                      ATNCB,ATCIDX,ATCOFF,
     +                      ATNLB,ATLIDX,ATLOFF,
     +                      ATNEB,ATEIDX,ATEOFF,
     +                      ATNYB,ATYIDX,ATYOFF,
     +                      ATNRB,ATRIDX,ATROFF,
     +                      BDNCEN,BDCEN,BDNBAS,BDBAS,
     +                      NBOBD,
     +                      LOCAL,
     +                      ANGNBO,
     +                      W,Q,B )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NBO_RESULTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints details of the set of obtained
C                natural bond orbitals (NBOs) to a file specified
C                by its unit identification number. When entering this
C                routine, the NBO expansion coefficients in terms of
C                the atomic NHOs have been ordered such that all
C                Lone-pair, Bond, Empty-pair, Antibond and Rydberg NBOs
C                are clustered together and each set is in the order
C                mentioned. Also the overall weight vector has been
C                ordered to match the order of the NBO expansion
C                coefficients in the NBO part, and additionally the
C                remaining Core and Rydberg NHOs not used for NBO
C                construction were placed in front and at the end,
C                respectively. A picture of what is present at this
C                stage helps:
C
C                Order of NBO expansion coefficients (B matrix):
C
C                        |  LP  |   B  |  EP  |  AB  |RY(nbo)|
C                          (NLB)  (NBB)  (NEB)  (NAB)  (NYB)
C
C                Order of overall weights (W array):
C
C                 | C(nho) |  LP  |   B  |  EP  |  AB  |RY(nbo)|RY(nho)|
C                    (NCB)  (NLB)  (NBB)  (NEB)  (NAB)   (NYB)   (NRB)
C
C                where the variables below in parenthesis indicate the
C                number of elements present in each set.
C
C                The following is printed:
C
C                   1) The occupation numbers (weights) and
C                      occupation interaction orders of the
C                      NBOs and their characterization as Core,
C                      Lone-pair (1-center), Bond and Antibond
C                      (x-center, x=2,3,...), and Rydberg type
C                      NBOs.
C
C                   2) A population chart indicating the NBO
C                      populations broken down in Lewis, Delocalized
C                      and Rest populations:
C
C                       Lewis: Core + Lone-pairs + 1- and 2-cen Bonds
C                       Delocalized: >2-cen Bonds
C                       Rest: Antibonds + Rydbergs
C
C                   3) The NBO atomic localization map.
C
C                   4) The NBO angular momentum weight table.
C
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    MXSHELL      =  largest l-shell value
C                    NxA          =  total # of Core, Lone-pair,
C                                    Empty-pair and Rydberg atoms
C                                    found (x=C,L,E,Y,R).
C                    NxB          =  total # of NHO Core, NBO Lone-
C                                    pairs, NBO Bonds, NBO Empty-pairs,
C                                    NBO Antibonds, NBO Rydbergs and NHO
C                                    Rydbergs found(x=C,L,B,E,A,Y,R).
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    ATNCB (A)    =  # of Core NHOs on core atom A.
C                    ATCIDX (A)   =  atomic index for core atom A.
C                    ATCOFF (A)   =  index offset for Core NHOs for
C                                    core atom A. This index is equal
C                                    to the total number of Core NHOs
C                                    on all core atoms preceeding core
C                                    atom A.
C                    ATNLB (A)    =  # of Lone-pair NBOs on Lone-pair
C                                    atom A.
C                    ATLIDX (A)   =  atomic index for Lone-pair atom A.
C                    ATLOFF (A)   =  index offset for Lone-pair NBOs for
C                                    Lone-pair atom A. This index is
C                                    equal to the total number of
C                                    Lone-pair NBOs on all Lone-pair
C                                    atoms preceeding Lone-pair atom A.
C                    ATNEB (A)    =  # of Empty-pair NBOs on Empty-pair
C                                    atom A.
C                    ATEIDX (A)   =  atomic index for Empty-pair atom A.
C                    ATEOFF (A)   =  index offset for Empty-pair NBOs for
C                                    Empty-pair atom A. This index is
C                                    equal to the total number of
C                                    Empty-pair NBOs on all Empty-pair
C                                    atoms preceeding Empty-pair atom A.
C                    ATNYB (A)    =  # of Rydberg NBOs on Rydberg atom A.
C                    ATYIDX (A)   =  atomic index for Rydberg atom A.
C                    ATYOFF (A)   =  index offset for Rydberg NBOs for
C                                    Rydberg atom A. This index is
C                                    equal to the total number of
C                                    Rydberg NBOs on all Rydberg atoms
C                                    preceeding Rydberg atom A.
C                    ATNRB (A)    =  # of Rydberg NHOs on core atom A.
C                    ATRIDX (A)   =  atomic index for Rydberg atom A.
C                    ATROFF (A)   =  index offset for Rydberg NHOs for
C                                    Rydberg atom A. This index is equal
C                                    to the total number of Rydberg NHOs
C                                    on all Rydberg atoms preceeding
C                                    Rydberg atom A.
C                    BDNCEN (J)   =  # of atomic centers for J-th bond.
C                                    The bond order is NHB-bonds/NCB/NRB.
C                                    Note that the # of NHB-bonds might
C                                    be < NHB size.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDNBAS (J)   =  # of basis functions (NHOs) for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    for the I-th NBO. This array is
C                                    the handle for accessing info of
C                                    the NHO bonds sitting in the arrays
C                                    BDNCEN,BDCEN,BDBAS and BDOCC. The
C                                    NBOBD elements are in NCB/NLB/NBB/
C                                    NEB/NAB/NYB/NRB order.
C                    LOCAL (A,I)  =  NBO atom locality map for each
C                                    atom A and with NBOs in NCB/NLB/
C                                    NBB/NEB/NAB/NYB/NRB order.
C                    ANGNBO (J,L) =  NBAS x (MXSHELL+1) matrix containing
C                                    the sum of the square of the NAO
C                                    coefficients per L-type angular
C                                    momentum for the J-th NBO. The NBO
C                                    order is NCB/NLB/NBB/NEB/NAB/NYB/
C                                    NRB.
C                    W            =  NBO weight vector with elements
C                                    in NCB/NLB/NBB/NEB/NAB/NYB/NRB
C                                    order.
C                    Q            =  NBO occupation interaction order
C                                    vector with elements in NCB/NLB/
C                                    NBB/NEB/NAB/NYB/NRB order.
C                    B (I,J)      =  NBOSIZE x NBAS matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO. The NBO
C                                    column index is in NCB/NLB/NBB/NEB/
C                                    NAB/NYB/NRB order.
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
         CHARACTER*12  NBOCHAR

         CHARACTER*2   ATSYMB  (1:104)
         CHARACTER*12  NBOSYMB (1:6)

         INTEGER     ATOM
         INTEGER     BOND
         INTEGER     I,J,N
         INTEGER     MXSHELL
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NBO
         INTEGER     NBOSIZE
         INTEGER     NCA,NLA,NEA,NYA,NRA
         INTEGER     NCB,NLB,NBB,NEB,NAB,NYB,NRB
         INTEGER     NCBA,NLBA,NEBA,NYBA,NRBA
         INTEGER     NCEN
         INTEGER     NHO,NNHO,NNHOOLD
         INTEGER     OFF
         INTEGER     OFFLP,OFFBB,OFFEP,OFFAB,OFFRY,OFFRR
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     ATNCB  (1:NCA  )
         INTEGER     ATNLB  (1:NLA  )
         INTEGER     ATNEB  (1:NEA  )
         INTEGER     ATNYB  (1:NYA  )
         INTEGER     ATNRB  (1:NRA  )
         INTEGER     ATCIDX (1:NCA  )
         INTEGER     ATLIDX (1:NLA  )
         INTEGER     ATEIDX (1:NEA  )
         INTEGER     ATYIDX (1:NYA  )
         INTEGER     ATRIDX (1:NRA  )
         INTEGER     ATCOFF (1:NCA  )
         INTEGER     ATLOFF (1:NLA  )
         INTEGER     ATEOFF (1:NEA  )
         INTEGER     ATYOFF (1:NYA  )
         INTEGER     ATROFF (1:NRA  )
         INTEGER     BDNBAS (1:NBAS )
         INTEGER     BDNCEN (1:NBAS )
         INTEGER     NBOBD  (1:NBAS )
         INTEGER     ZATOM  (1:NATOM)

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS)
         INTEGER     BDCEN   (1:NBOSIZE,1:NBAS)

         DOUBLE PRECISION  COEFF
         DOUBLE PRECISION  LEWIS,DELOC,REST
         DOUBLE PRECISION  LOCVAL
         DOUBLE PRECISION  ORDER
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  Q (1:NBAS)
         DOUBLE PRECISION  W (1:NBAS)

         DOUBLE PRECISION  B      (1:NBOSIZE,1:NBAS   )
         DOUBLE PRECISION  LOCAL  (1:NATOM  ,1:NBAS   )
         DOUBLE PRECISION  ANGNBO (1:NBAS   ,0:MXSHELL)

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
         DATA  NBOSYMB /'    Core    ',
     +                  '  Lone-pair ',
     +                  '    Bond    ',
     +                  ' Empty-pair ',
     +                  '  Anti-bond ',
     +                  '   Rydberg  '/

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...print out title.
C
C
         WRITE (UNITID,8000) 'Natural bond orbital analysis'
 8000    FORMAT (//,30X,A29,//)
C
C
C             ...print out header line.
C
C
         WRITE (UNITID,9000) 'NBO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'NHO #',BAR,'NHO coeffs',BAR,
     +                       'NBO-type',BAR,'Occupancy',BAR,
     +                       'Int. Orders',BAR,'Locality'
         WRITE (UNITID,9010) (DASH,I=1,104)
C
C
C             ...print out Core NAO information.
C
C
         LEWIS = ZERO
         DELOC = ZERO
         REST  = ZERO

         IF (NCA.GT.0) THEN

             NBOCHAR = NBOSYMB (1)

             DO 1000 N = 1,NCA
                OFF = ATCOFF (N)
                NCBA = ATNCB (N)
                ATOM = ATCIDX (N)
                ZVAL = MIN (ZATOM (ATOM),104)
                ATCHAR = ATSYMB (ZVAL)
                NNHOOLD = 0

                DO J = 1,NCBA
                   NBO = OFF + J
                   BOND = NBOBD (NBO)
                   NNHO = BDNBAS (BOND)
                   WEIGHT = W (NBO)
                   ORDER  = Q (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   LEWIS = LEWIS + WEIGHT

                   IF (NNHO.EQ.1) THEN
                       NHO = BDBAS (1,BOND)
                       COEFF = B (1,NBO)
                       IF (J.EQ.1) THEN
                           WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,ATOM,
     +                                         BAR,NHO,BAR,COEFF,BAR,
     +                                         NBOCHAR,BAR,WEIGHT,BAR,
     +                                         ORDER,BAR,LOCVAL
                       ELSE
                           WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,BAR,
     +                                         COEFF,BAR,NBOCHAR,BAR,
     +                                         WEIGHT,BAR,ORDER,BAR,
     +                                         LOCVAL
                       END IF
                   ELSE
                       IF (NNHOOLD.EQ.1) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                       DO I = 1,NNHO
                          NHO = BDBAS (I,BOND)
                          COEFF = B (I,NBO)

                          IF (J.EQ.1 .AND. I.EQ.1) THEN
                              WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,
     +                                            ATOM,BAR,NHO,BAR,
     +                                            COEFF,BAR,NBOCHAR,BAR,
     +                                            WEIGHT,BAR,ORDER,BAR,
     +                                            LOCVAL
                          ELSE IF (I.EQ.1) THEN
                              WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,NBOCHAR,
     +                                            BAR,WEIGHT,BAR,ORDER,
     +                                            BAR,LOCVAL
                          ELSE
                              WRITE (UNITID,9920) BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,BAR,
     +                                            BAR,BAR
                          END IF
                       END DO
                       IF (J.LT.NCBA) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                   END IF

                   NNHOOLD = NNHO

                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 1000        CONTINUE
         END IF
C
C
C             ...print out Lone-pair NBO information.
C
C
         IF (NLB.GT.0) THEN

             OFFLP = NCB
             NBOCHAR = NBOSYMB (2)

             DO 2000 N = 1,NLA
                OFF = ATLOFF (N)
                NLBA = ATNLB (N)
                ATOM = ATLIDX (N)
                ZVAL = MIN (ZATOM (ATOM),104)
                ATCHAR = ATSYMB (ZVAL)
                NNHOOLD = 0

                DO J = 1,NLBA
                   NBO = OFFLP + OFF + J
                   BOND = NBOBD (NBO)
                   NNHO = BDNBAS (BOND)
                   WEIGHT = W (NBO)
                   ORDER  = Q (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   LEWIS = LEWIS + WEIGHT

                   IF (NNHO.EQ.1) THEN
                       NHO = BDBAS (1,BOND)
                       COEFF = B (1,NBO)
                       IF (J.EQ.1) THEN
                           WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,ATOM,
     +                                         BAR,NHO,BAR,COEFF,BAR,
     +                                         NBOCHAR,BAR,WEIGHT,BAR,
     +                                         ORDER,BAR,LOCVAL
                       ELSE
                           WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,BAR,
     +                                         COEFF,BAR,NBOCHAR,BAR,
     +                                         WEIGHT,BAR,ORDER,BAR,
     +                                         LOCVAL
                       END IF
                   ELSE
                       IF (NNHOOLD.EQ.1) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                       DO I = 1,NNHO
                          NHO = BDBAS (I,BOND)
                          COEFF = B (I,NBO)

                          IF (J.EQ.1 .AND. I.EQ.1) THEN
                              WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,
     +                                            ATOM,BAR,NHO,BAR,
     +                                            COEFF,BAR,NBOCHAR,BAR,
     +                                            WEIGHT,BAR,ORDER,BAR,
     +                                            LOCVAL
                          ELSE IF (I.EQ.1) THEN
                              WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,NBOCHAR,
     +                                            BAR,WEIGHT,BAR,ORDER,
     +                                            BAR,LOCVAL
                          ELSE
                              WRITE (UNITID,9920) BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,BAR,
     +                                            BAR,BAR
                          END IF
                       END DO
                       IF (J.LT.NLBA) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                   END IF

                   NNHOOLD = NNHO

                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 2000        CONTINUE
         END IF
C
C
C             ...print out Bond NBO information.
C
C
         IF (NBB.GT.0) THEN

             OFFBB = NCB + NLB
             NBOCHAR = NBOSYMB (3)

             DO 3000 N = 1,NBB
                NBO = OFFBB + N
                BOND = NBOBD (NBO)
                NCEN = BDNCEN (BOND)
                WEIGHT = W (NBO)
                ORDER  = Q (NBO)

                IF (NCEN.EQ.2) THEN
                    LEWIS = LEWIS + WEIGHT
                ELSE
                    DELOC = DELOC + WEIGHT
                END IF

                LOCVAL = ZERO
                DO I = 1,NCEN
                   ATOM = BDCEN (I,BOND)
                   LOCVAL = LOCVAL + LOCAL (ATOM,NBO)
                END DO

                DO I = 1,NCEN
                   NHO = BDBAS (I,BOND)
                   ATOM = BDCEN (I,BOND)
                   COEFF = B (I,NBO)
                   ZVAL = MIN (ZATOM (ATOM),104)
                   ATCHAR = ATSYMB (ZVAL)

                   IF (I.EQ.1) THEN
                       WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,NBOCHAR,
     +                                     BAR,WEIGHT,BAR,ORDER,BAR,
     +                                     LOCVAL
                   ELSE
                       WRITE (UNITID,9050) BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,BAR,
     +                                     BAR,BAR
                   END IF
                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 3000        CONTINUE
         END IF
C
C
C             ...print out Empty-pair NBO information.
C
C
         IF (NEB.GT.0) THEN

             OFFEP = NCB + NLB + NBB
             NBOCHAR = NBOSYMB (4)

             DO 4000 N = 1,NEA
                OFF = ATEOFF (N)
                NEBA = ATNEB (N)
                ATOM = ATEIDX (N)
                ZVAL = MIN (ZATOM (ATOM),104)
                ATCHAR = ATSYMB (ZVAL)
                NNHOOLD = 0

                DO J = 1,NEBA
                   NBO = OFFEP + OFF + J
                   BOND = NBOBD (NBO)
                   NNHO = BDNBAS (BOND)
                   WEIGHT = W (NBO)
                   ORDER  = Q (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   REST = REST + WEIGHT

                   IF (NNHO.EQ.1) THEN
                       NHO = BDBAS (1,BOND)
                       COEFF = B (1,NBO)
                       IF (J.EQ.1) THEN
                           WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,ATOM,
     +                                         BAR,NHO,BAR,COEFF,BAR,
     +                                         NBOCHAR,BAR,WEIGHT,BAR,
     +                                         ORDER,BAR,LOCVAL
                       ELSE
                           WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,BAR,
     +                                         COEFF,BAR,NBOCHAR,BAR,
     +                                         WEIGHT,BAR,ORDER,BAR,
     +                                         LOCVAL
                       END IF
                   ELSE
                       IF (NNHOOLD.EQ.1) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                       DO I = 1,NNHO
                          NHO = BDBAS (I,BOND)
                          COEFF = B (I,NBO)

                          IF (J.EQ.1 .AND. I.EQ.1) THEN
                              WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,
     +                                            ATOM,BAR,NHO,BAR,
     +                                            COEFF,BAR,NBOCHAR,BAR,
     +                                            WEIGHT,BAR,ORDER,BAR,
     +                                            LOCVAL
                          ELSE IF (I.EQ.1) THEN
                              WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,NBOCHAR,
     +                                            BAR,WEIGHT,BAR,ORDER,
     +                                            BAR,LOCVAL
                          ELSE
                              WRITE (UNITID,9920) BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,BAR,
     +                                            BAR,BAR
                          END IF
                       END DO
                       IF (J.LT.NEBA) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                   END IF

                   NNHOOLD = NNHO

                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 4000        CONTINUE
         END IF
C
C
C             ...print out Antibond NBO information.
C
C
         IF (NAB.GT.0) THEN

             OFFAB = NCB + NLB + NBB + NEB
             NBOCHAR = NBOSYMB (5)

             DO 5000 N = 1,NAB
                NBO = OFFAB + N
                BOND = NBOBD (NBO)
                NCEN = BDNCEN (BOND)
                WEIGHT = W (NBO)
                ORDER  = Q (NBO)
                REST = REST + WEIGHT

                LOCVAL = ZERO
                DO I = 1,NCEN
                   ATOM = BDCEN (I,BOND)
                   LOCVAL = LOCVAL + LOCAL (ATOM,NBO)
                END DO

                DO I = 1,NCEN
                   NHO = BDBAS (I,BOND)
                   ATOM = BDCEN (I,BOND)
                   COEFF = B (I,NBO)
                   ZVAL = MIN (ZATOM (ATOM),104)
                   ATCHAR = ATSYMB (ZVAL)

                   IF (I.EQ.1) THEN
                       WRITE (UNITID,9040) NBO,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,NBOCHAR,
     +                                     BAR,WEIGHT,BAR,ORDER,BAR,
     +                                     LOCVAL
                   ELSE
                       WRITE (UNITID,9050) BAR,ATCHAR,BAR,ATOM,BAR,
     +                                     NHO,BAR,COEFF,BAR,BAR,
     +                                     BAR,BAR
                   END IF
                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 5000        CONTINUE
         END IF
C
C
C             ...print out Rydberg NBO information.
C
C
         IF (NYB.GT.0) THEN

             OFFRY = NCB + NLB + NBB + NEB + NAB
             NBOCHAR = NBOSYMB (6)

             DO 6000 N = 1,NYA
                OFF = ATYOFF (N)
                NYBA = ATNYB (N)
                ATOM = ATYIDX (N)
                ZVAL = MIN (ZATOM (ATOM),104)
                ATCHAR = ATSYMB (ZVAL)
                NNHOOLD = 0

                DO J = 1,NYBA
                   NBO = OFFRY + OFF + J
                   BOND = NBOBD (NBO)
                   NNHO = BDNBAS (BOND)
                   WEIGHT = W (NBO)
                   ORDER  = Q (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   REST = REST + WEIGHT

                   IF (NNHO.EQ.1) THEN
                       NHO = BDBAS (1,BOND)
                       COEFF = B (1,NBO)
                       IF (J.EQ.1) THEN
                           WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,ATOM,
     +                                         BAR,NHO,BAR,COEFF,BAR,
     +                                         NBOCHAR,BAR,WEIGHT,BAR,
     +                                         ORDER,BAR,LOCVAL
                       ELSE
                           WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,BAR,
     +                                         COEFF,BAR,NBOCHAR,BAR,
     +                                         WEIGHT,BAR,ORDER,BAR,
     +                                         LOCVAL
                       END IF
                   ELSE
                       IF (NNHOOLD.EQ.1) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                       DO I = 1,NNHO
                          NHO = BDBAS (I,BOND)
                          COEFF = B (I,NBO)

                          IF (J.EQ.1 .AND. I.EQ.1) THEN
                              WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,
     +                                            ATOM,BAR,NHO,BAR,
     +                                            COEFF,BAR,NBOCHAR,BAR,
     +                                            WEIGHT,BAR,ORDER,BAR,
     +                                            LOCVAL
                          ELSE IF (I.EQ.1) THEN
                              WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,NBOCHAR,
     +                                            BAR,WEIGHT,BAR,ORDER,
     +                                            BAR,LOCVAL
                          ELSE
                              WRITE (UNITID,9920) BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,BAR,
     +                                            BAR,BAR
                          END IF
                       END DO
                       IF (J.LT.NYBA) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                   END IF

                   NNHOOLD = NNHO

                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 6000        CONTINUE
         END IF
C
C
C             ...print out Rydberg NBO information.
C
C
         IF (NRA.GT.0) THEN

             OFFRR = NCB + NLB + NBB + NEB + NAB + NYB
             NBOCHAR = NBOSYMB (6)

             DO 7000 N = 1,NRA

                OFF = ATROFF (N)
                NRBA = ATNRB (N)
                ATOM = ATRIDX (N)
                ZVAL = MIN (ZATOM (ATOM),104)
                ATCHAR = ATSYMB (ZVAL)
                NNHOOLD = 0

                DO J = 1,NRBA
                   NBO = OFFRR + OFF + J
                   BOND = NBOBD (NBO)
                   NNHO = BDNBAS (BOND)
                   WEIGHT = W (NBO)
                   ORDER  = Q (NBO)
                   LOCVAL = LOCAL (ATOM,NBO)
                   REST = REST + WEIGHT

                   IF (NNHO.EQ.1) THEN
                       NHO = BDBAS (1,BOND)
                       COEFF = B (1,NBO)
                       IF (J.EQ.1) THEN
                           WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,ATOM,
     +                                         BAR,NHO,BAR,COEFF,BAR,
     +                                         NBOCHAR,BAR,WEIGHT,BAR,
     +                                         ORDER,BAR,LOCVAL
                       ELSE
                           WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,BAR,
     +                                         COEFF,BAR,NBOCHAR,BAR,
     +                                         WEIGHT,BAR,ORDER,BAR,
     +                                         LOCVAL
                       END IF
                   ELSE
                       IF (NNHOOLD.EQ.1) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                       DO I = 1,NNHO
                          NHO = BDBAS (I,BOND)
                          COEFF = B (I,NBO)

                          IF (J.EQ.1 .AND. I.EQ.1) THEN
                              WRITE (UNITID,9900) NBO,BAR,ATCHAR,BAR,
     +                                            ATOM,BAR,NHO,BAR,
     +                                            COEFF,BAR,NBOCHAR,BAR,
     +                                            WEIGHT,BAR,ORDER,BAR,
     +                                            LOCVAL
                          ELSE IF (I.EQ.1) THEN
                              WRITE (UNITID,9910) NBO,BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,NBOCHAR,
     +                                            BAR,WEIGHT,BAR,ORDER,
     +                                            BAR,LOCVAL
                          ELSE
                              WRITE (UNITID,9920) BAR,BAR,BAR,NHO,
     +                                            BAR,COEFF,BAR,BAR,
     +                                            BAR,BAR
                          END IF
                       END DO
                       IF (J.LT.NRBA) THEN
                           WRITE (UNITID,9930) BAR,BAR,BAR,BAR,
     +                                         BAR,BAR,BAR,BAR
                       END IF
                   END IF

                   NNHOOLD = NNHO

                END DO

                WRITE (UNITID,9010) (DASH,I=1,104)
 7000        CONTINUE
         END IF
C
C
C             ...formats for printing the NBO information.
C
C
 9000    FORMAT (1X,A5,1X,     A1, 1X,A4,1X,    A1, 1X,A6,1X,  A1,
     +           1X,A5,1X,     A1, 1X,A10,1X,   A1, 3X,A8,3X,  A1,
     +           3X,A9,3X,     A1, 2X,A11,2X,   A1, 2X,A8,2X)
 9010    FORMAT (1X,104A1)
 9040    FORMAT (2X,I3,2X,     A1, 2X,A2,2X,    A1, 2X,I3,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,A12,1X, A1,
     +           2X,F12.10,1X, A1, 2X,F12.10,1X,A1, 2X,F8.4,2X)
 9050    FORMAT (2X,3X,2X,     A1, 2X,A2,2X,    A1, 2X,I3,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,12X,1X, A1,
     +           2X,12X,1X,    A1, 2X,12X,1X,   A1)
 9900    FORMAT (2X,I3,2X,     A1, 2X,A2,2X,    A1, 2X,I3,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,A12,1X, A1,
     +           2X,F12.10,1X, A1, 2X,F12.10,1X,A1, 2X,F8.4,2X)
 9910    FORMAT (2X,I3,2X,     A1, 2X,2X,2X,    A1, 2X,3X,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,A12,1X, A1,
     +           2X,F12.10,1X, A1, 2X,F12.10,1X,A1, 2X,F8.4,2X)
 9920    FORMAT (2X,3X,2X,     A1, 2X,2X,2X,    A1, 2X,3X,3X,  A1,
     +           2X,I3,2X,     A1, 1X,F10.6,1X, A1, 1X,12X,1X, A1,
     +           2X,12X,1X,    A1, 2X,12X,1X,   A1, 2X,8X,2X)
 9930    FORMAT (2X,3X,2X,     A1, 2X,2X,2X,    A1, 2X,3X,3X,  A1,
     +           2X,3X,2X,     A1, 1X,10X,1X,   A1, 1X,12X,1X, A1,
     +           2X,12X,1X,    A1, 2X,12X,1X,   A1, 2X,8X,2X)
C
C
C             ...print out population chart.
C
C
         WRITE (UNITID,9060) 'Population Chart'
 9060    FORMAT (//,20X,A16,//)

         WRITE (UNITID,9070) 'NBO-type',BAR,'Occupancy'
         WRITE (UNITID,9080) (DASH,I=1,28)
         WRITE (UNITID,9090) 'Lewis',BAR,LEWIS
         WRITE (UNITID,9100) 'Delocalized',BAR,DELOC
         WRITE (UNITID,9110) 'Rest',BAR,REST
         WRITE (UNITID,9080) (DASH,I=1,28)
         WRITE (UNITID,9120) 'Total',BAR,LEWIS+DELOC+REST

 9070    FORMAT (3X,A8,2X,  A1, 4X,A9)
 9080    FORMAT (2X,28A1)
 9090    FORMAT (7X,A5,1X,  A1, 1X,F12.8)
 9100    FORMAT (1X,A11,1X, A1, 1X,F12.8)
 9110    FORMAT (8X,A4,1X,  A1, 1X,F12.8)
 9120    FORMAT (7X,A5,1X,  A1, 1X,F12.8)
C
C
C             ...print out the NBO angular momentum weight table.
C
C
         WRITE (UNITID,9130) 'NBO Angular Momentum Weight Table'
 9130    FORMAT (//,20X,A33)

         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows = NBOs ; Cols = L-value (1=s,2=p,...) ',
     +                  NBAS,MXSHELL+1,
     +                  NBAS,MXSHELL+1,
     +                  ANGNBO )
     +
     +
C
C
C             ...print out the NBO atomic localization map.
C
C
         WRITE (UNITID,9140) 'NBO Atomic Localization Map'
 9140    FORMAT (//,20X,A27)

         CALL    MAT__PRINT_A_FLOAT_3_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows = Atoms ; Columns = NBOs ',
     +                  NATOM,NBAS,
     +                  NATOM,NBAS,
     +                  LOCAL )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
