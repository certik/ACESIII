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
         SUBROUTINE  NLO__ROTATE_TGROUP_NBO
     +
     +                    ( NBAS,
     +                      NBOSIZE,
     +                      NXB,NXA,NXBA,
     +                      ATXOFF,
     +                      LSYMACC,QSYMACC,
     +                      NPAIR,
     +                      IDXSYM,
     +                      NBOBD,
     +                      P,
     +                      LOCAL,
     +                      W,Q,
     +                      XVEC,
     +                      XMAT,
     +
     +                              SYMNBO,
     +                              SLTYPE,
     +                              BDNBAS,
     +                              BDBAS,
     +                              B,
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ROTATE_TGROUP_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine rotates a tetrahedral set (four) of
C                two-fold interaction order degenerate atomic NBOs,
C                such that they become equivalent under operations of
C                the tetrahedral groups T,Td and Th.
C
C                The first step is to divide the set into two pairs
C                and rotating each pair of degenerate atomic NBOs such
C                that their lobes match and form rectangular angles
C                between them (for the mathematics on how this is done,
C                please see the corresponding NBO pair rotation
C                routine):
C
C
C
C                                                  x
C                       x  *                       *
C                        x *    x                  x
C                         x* x                     *
C                    ******x******   ---->   x**x**x**x**x
C                        x *x                      *
C                     x    * x                     x
C                          *  x                    *
C                                                  x
C
C
C
C                In the tetrahedral case such an arrangement of the
C                lobes is always well defined, because of the existence
C                of a nonzero tilt angle between them. If the lobes
C                would be opposite then there is obviously an infinite
C                amount of such arrangements in space.
C
C                The second step tries to construct sp2 hybrids that
C                will transform properly under the four C3 axis of the
C                tetrahedral groups.
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    NXB          =  total # of atomic x-type NBOs.
C                    NXA          =  total # of atomic x-type NBO atoms.
C                    NXBA         =  # of atomic x-type NBOs on one of
C                                    the tetrahedrally related atoms.
C                    ATXOFF (A)   =  index offset for atomic NBOs for
C                                    x-type NBO atom A. This index is
C                                    equal to the total number of atomic
C                                    x-type NBOs on all x-type NBO atoms
C                                    preceeding x-type NBO atom A.
C                    LSYMACC      =  symmetry accuracy for orbital
C                                    localization contents
C                    QSYMACC      =  symmetry accuracy for orbital
C                                    interaction order values (sum of
C                                    squares of offdiagonal occupation
C                                    matrix elements for one orbital)
C                    NPAIR        =  local atomic index of the two-fold
C                                    weight degenerate atomic NBOs.
C                                    For each atom the local atomic
C                                    indices for the degenerate pair
C                                    are NPAIR and NPAIR+1.
C                    IDXSYM       =  array that contains the four
C                                    tetrahedrally related x-type NBO
C                                    atomic indices.
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    for the I-th atomic x-type NBO.
C                                    This array is the handle for
C                                    accessing and modifying info of
C                                    the NHO bonds sitting in the
C                                    arrays BDNCEN,BDCEN,BDBAS and
C                                    BDOCC.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    LOCAL        =  atomic localization content for
C                                    all NBOs on one of the 4 atoms.
C                    W            =  NBO weight vector at present stage
C                                    of symmetrization for all NBOs
C                                    on one of the 4 atoms.
C                    Q            =  NBO interaction order vector at
C                                    present stage of symmetrization
C                                    for all NBOs on one of the 4 atoms.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    SYMNBO       =  initial integer vector indicating
C                                    the different symmetry stages of
C                                    each NBO (=0 means no symmetry
C                                    adaptation, =1 means hybridization
C                                    was performed, =2 means symmetry
C                                    adapted in final form). This array
C                                    is extremely useful for deciding
C                                    which NBOs can be used for the
C                                    symmetry adaptation at different
C                                    stages. Note that only NBOs with
C                                    SYMNBO = 2 can be considered for
C                                    that purpose.
C                    SLTYPE       =  initial integer vector indicating
C                                    the smallest angular momentum
C                                    component present in the NBOs.
C                                    Useful for deciding which atomic
C                                    NBOs to use for hybrid formation.
C                                    Once an atomic NBO has been used
C                                    for hybrid formation, a -1 will
C                                    be placed at the corresponding
C                                    vector position and the NBO cannot
C                                    be used again for hybrid formation.
C                    BDNBAS (J)   =  initial # of basis functions (NHOs)
C                                    for J-th bond. The bond order is
C                                    NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    BDBAS (I,J)  =  initial I-th global basis (NHO)
C                                    index for J-th bond. The bond order
C                                    is NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    B (I,J)      =  initial NBOSIZE x NXB matrix
C                                    containing the J-th symmetry
C                                    adapted x-type NBO expansion
C                                    coefficients in terms of the I-th
C                                    atomic NHOs forming the J-th
C                                    symmetry adapted x-type NBO.
C                    C            =  initial NBAS x NXB section of the
C                                    NBO coefficient matrix in AO basis
C                                    before symmetrization corresponding
C                                    to the x-type NBOs treated here.
C
C
C                  Output:
C
C                    SYMNBO       =  updated integer vector indicating
C                                    the different symmetry stages of
C                                    each NBO (=0 means no symmetry
C                                    adaptation, =1 means hybridization
C                                    was performed, =2 means symmetry
C                                    adapted in final form). This array
C                                    is extremely useful for deciding
C                                    which NBOs can be used for the
C                                    symmetry adaptation at different
C                                    stages. Note that only NBOs with
C                                    SYMNBO = 2 can be considered for
C                                    that purpose.
C                    SLTYPE       =  updated integer vector indicating
C                                    the smallest angular momentum
C                                    component present in the NBOs.
C                                    Useful for deciding which atomic
C                                    NBOs to use for hybrid formation.
C                                    Once an atomic NBO has been used
C                                    for hybrid formation, a -1 will
C                                    be placed at the corresponding
C                                    vector position and the NBO cannot
C                                    be used again for hybrid formation.
C                    BDNBAS (J)   =  modified # of basis functions
C                                    (NHOs) for J-th bond. The bond
C                                    order is NHB-bonds/NCB/NRB. Note
C                                    that the # of NHB-bonds might
C                                    be < NHB size.
C                    BDBAS (I,J)  =  modified I-th global basis (NHO)
C                                    index for J-th bond. The bond order
C                                    is NHB-bonds/NCB/NRB. Note that the
C                                    # of NHB-bonds might be < NHB size.
C                    B (I,J)      =  modified NBOSIZE x NXB matrix
C                                    containing the J-th symmetry
C                                    adapted x-type NBO expansion
C                                    coefficients in terms of the I-th
C                                    atomic NHOs forming the J-th
C                                    symmetry adapted x-type NBO.
C                    C            =  modified NBAS x NXB section of the
C                                    NBO coefficient matrix in AO basis
C                                    corresponding to the symmetrized
C                                    x-type NBOs treated here.
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

         LOGICAL     HYBRID

         INTEGER     BDS1,BDS2,BDS3,BDS4
         INTEGER     BDP1X,BDP2X,BDP3X,BDP4X
         INTEGER     BDP1Y,BDP2Y,BDP3Y,BDP4Y
         INTEGER     NBAS
         INTEGER     NBOSIZE
         INTEGER     NHOS1,NHOS2,NHOS3,NHOS4
         INTEGER     NHOP1X,NHOP2X,NHOP3X,NHOP4X
         INTEGER     NHOP1Y,NHOP2Y,NHOP3Y,NHOP4Y
         INTEGER     NLO__ANGULAR_MOMENTUM_NBO
         INTEGER     NPAIR
         INTEGER     NXB,NXA,NXBA
         INTEGER     OFF1,OFF2,OFF3,OFF4
         INTEGER     P1,P2,P3,P4
         INTEGER     S1,S2,S3,S4
         INTEGER     TOTSYM

         INTEGER     ATXOFF  (1:NXA )
         INTEGER     BDNBAS  (1:NBAS)
         INTEGER     IDXSYM  (1:4   )
         INTEGER     NBOBD   (1:NXB )
         INTEGER     SLTYPE  (1:NXB )
         INTEGER     SYMNBO  (1:NXB )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS)

         DOUBLE PRECISION  COS1,COS2,COS3,COS4
         DOUBLE PRECISION  LSYMACC,QSYMACC
         DOUBLE PRECISION  SIN1,SIN2,SIN3,SIN4

         DOUBLE PRECISION  LOCAL (1:NXBA)
         DOUBLE PRECISION  Q     (1:NXBA)
         DOUBLE PRECISION  W     (1:NXBA)
         DOUBLE PRECISION  XVEC  (1:NBAS)

         DOUBLE PRECISION  HYB1  (1:3,1:3)
         DOUBLE PRECISION  HYB2  (1:3,1:3)
         DOUBLE PRECISION  HYB3  (1:3,1:3)
         DOUBLE PRECISION  HYB4  (1:3,1:3)

         DOUBLE PRECISION  B     (1:NBOSIZE,1:NXB )
         DOUBLE PRECISION  C     (1:NBAS   ,1:NXB )
         DOUBLE PRECISION  P     (1:NBAS   ,1:NBAS)
         DOUBLE PRECISION  XMAT  (1:NBAS   ,1:2   )
C
C
C------------------------------------------------------------------------
C
C
C             ...rotate first the pairs into appropriate positions
C                by calling twice the pair rotation routine.
C
C
         OFF1 = ATXOFF (IDXSYM (1))
         OFF2 = ATXOFF (IDXSYM (2))
         OFF3 = ATXOFF (IDXSYM (3))
         OFF4 = ATXOFF (IDXSYM (4))

         P1 = OFF1 + NPAIR
         P2 = OFF2 + NPAIR
         P3 = OFF3 + NPAIR
         P4 = OFF4 + NPAIR

         BDP1X = NBOBD (P1  )
         BDP1Y = NBOBD (P1+1)
         BDP2X = NBOBD (P2  )
         BDP2Y = NBOBD (P2+1)
         BDP3X = NBOBD (P3  )
         BDP3Y = NBOBD (P3+1)
         BDP4X = NBOBD (P4  )
         BDP4Y = NBOBD (P4+1)

         NHOP1X = BDBAS (1,BDP1X)
         NHOP1Y = BDBAS (1,BDP1Y)
         NHOP2X = BDBAS (1,BDP2X)
         NHOP2Y = BDBAS (1,BDP2Y)
         NHOP3X = BDBAS (1,BDP3X)
         NHOP3Y = BDBAS (1,BDP3Y)
         NHOP4X = BDBAS (1,BDP4X)
         NHOP4Y = BDBAS (1,BDP4Y)

         CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +              ( NBAS,
     +                P,
     +                XVEC,
     +                XMAT,
     +
     +                       COS1,SIN1,
     +                       COS2,SIN2,
     +                       C (1,P1),
     +                       C (1,P2) )
     +
     +
         CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +              ( NBAS,
     +                P,
     +                XVEC,
     +                XMAT,
     +
     +                       COS3,SIN3,
     +                       COS4,SIN4,
     +                       C (1,P3),
     +                       C (1,P4) )
     +
     +
         B (1,P1  ) =   COS1
         B (2,P1  ) = - SIN1
         B (1,P1+1) =   SIN1
         B (2,P1+1) =   COS1
         B (1,P2  ) =   COS2
         B (2,P2  ) = - SIN2
         B (1,P2+1) =   SIN2
         B (2,P2+1) =   COS2
         B (1,P3  ) =   COS3
         B (2,P3  ) = - SIN3
         B (1,P3+1) =   SIN3
         B (2,P3+1) =   COS3
         B (1,P4  ) =   COS4
         B (2,P4  ) = - SIN4
         B (1,P4+1) =   SIN4
         B (2,P4+1) =   COS4

         SYMNBO (P1  ) = 2
         SYMNBO (P1+1) = 2
         SYMNBO (P2  ) = 2
         SYMNBO (P2+1) = 2
         SYMNBO (P3  ) = 2
         SYMNBO (P3+1) = 2
         SYMNBO (P4  ) = 2
         SYMNBO (P4+1) = 2

         BDNBAS (BDP1X) = 2
         BDNBAS (BDP1Y) = 2
         BDNBAS (BDP2X) = 2
         BDNBAS (BDP2Y) = 2
         BDNBAS (BDP3X) = 2
         BDNBAS (BDP3Y) = 2
         BDNBAS (BDP4X) = 2
         BDNBAS (BDP4Y) = 2

         BDBAS (1,BDP1X) = NHOP1X
         BDBAS (2,BDP1X) = NHOP1Y
         BDBAS (1,BDP1Y) = NHOP1X
         BDBAS (2,BDP1Y) = NHOP1Y
         BDBAS (1,BDP2X) = NHOP2X
         BDBAS (2,BDP2X) = NHOP2Y
         BDBAS (1,BDP2Y) = NHOP2X
         BDBAS (2,BDP2Y) = NHOP2Y
         BDBAS (1,BDP3X) = NHOP3X
         BDBAS (2,BDP3X) = NHOP3Y
         BDBAS (1,BDP3Y) = NHOP3X
         BDBAS (2,BDP3Y) = NHOP3Y
         BDBAS (1,BDP4X) = NHOP4X
         BDBAS (2,BDP4X) = NHOP4Y
         BDBAS (1,BDP4Y) = NHOP4X
         BDBAS (2,BDP4Y) = NHOP4Y
C
C
C             ...attempt sp2 hybridization. At the moment we look
C                for totally symmetric NBOs in the vicinity of the
C                degenerate pair. If no such totally symmetric NBO
C                can be found, the hybridization is abandoned.
C
C
         TOTSYM = NLO__ANGULAR_MOMENTUM_NBO
     +
     +                 ( NXBA,
     +                   SLTYPE (OFF1+1),
     +                   NPAIR,NPAIR+1,
     +                   0,
     +                   1,
     +                   QSYMACC,LSYMACC,
     +                   W,Q,
     +                   LOCAL )
     +
     +
         IF (TOTSYM.EQ.0) THEN

             TOTSYM = NLO__ANGULAR_MOMENTUM_NBO
     +
     +                     ( NXBA,
     +                       SLTYPE (OFF1+1),
     +                       NPAIR,NPAIR+1,
     +                       1,
     +                       1,
     +                       QSYMACC,LSYMACC,
     +                       W,Q,
     +                       LOCAL )
     +
     +
         END IF

         IF (TOTSYM.EQ.0) THEN

             TOTSYM = NLO__ANGULAR_MOMENTUM_NBO
     +
     +                     ( NXBA,
     +                       SLTYPE (OFF1+1),
     +                       NPAIR,NPAIR+1,
     +                       2,
     +                       1,
     +                       QSYMACC,LSYMACC,
     +                       W,Q,
     +                       LOCAL )
     +
     +
         END IF

         HYBRID = TOTSYM .NE. 0

         IF (HYBRID) THEN

             IF (NBOSIZE.LT.3) THEN
                 WRITE (*,*) ' Cannot form sp2 hybrids! '
                 WRITE (*,*) ' NBOSIZE = ',NBOSIZE,' (> 2 required!) '
                 WRITE (*,*) ' nlo__rotate_tgroup_nbo '
                 WRITE (1,*) ' Cannot form sp2 hybrids! '
                 WRITE (1,*) ' NBOSIZE = ',NBOSIZE,' (> 2 required!) '
                 WRITE (1,*) ' nlo__rotate_tgroup_nbo '
                 STOP
             END IF

             S1 = OFF1 + TOTSYM
             S2 = OFF2 + TOTSYM
             S3 = OFF3 + TOTSYM
             S4 = OFF4 + TOTSYM

             BDS1 = NBOBD (S1)
             BDS2 = NBOBD (S2)
             BDS3 = NBOBD (S3)
             BDS4 = NBOBD (S4)

             NHOS1 = BDBAS (1,BDS1)
             NHOS2 = BDBAS (1,BDS2)
             NHOS3 = BDBAS (1,BDS3)
             NHOS4 = BDBAS (1,BDS4)

             CALL  NLO__FORM_SP2_PAIR_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +
     +                           HYB1,HYB2,
     +                           C (1,S1),
     +                           C (1,P1),
     +                           C (1,S2),
     +                           C (1,P2) )
     +
     +
             CALL  NLO__FORM_SP2_PAIR_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +
     +                           HYB3,HYB4,
     +                           C (1,S3),
     +                           C (1,P3),
     +                           C (1,S4),
     +                           C (1,P4) )
     +
     +
             B (1,S1  ) =          HYB1 (1,1)
             B (2,S1  ) =   COS1 * HYB1 (2,1) + SIN1 * HYB1 (3,1)
             B (3,S1  ) = - SIN1 * HYB1 (2,1) + COS1 * HYB1 (3,1)
             B (1,P1  ) =          HYB1 (1,2)
             B (2,P1  ) =   COS1 * HYB1 (2,2) + SIN1 * HYB1 (3,2)
             B (3,P1  ) = - SIN1 * HYB1 (2,2) + COS1 * HYB1 (3,2)
             B (1,P1+1) =          HYB1 (1,3)
             B (2,P1+1) =   COS1 * HYB1 (2,3) + SIN1 * HYB1 (3,3)
             B (3,P1+1) = - SIN1 * HYB1 (2,3) + COS1 * HYB1 (3,3)

             B (1,S2  ) =          HYB2 (1,1)
             B (2,S2  ) =   COS2 * HYB2 (2,1) + SIN2 * HYB2 (3,1)
             B (3,S2  ) = - SIN2 * HYB2 (2,1) + COS2 * HYB2 (3,1)
             B (1,P2  ) =          HYB2 (1,2)
             B (2,P2  ) =   COS2 * HYB2 (2,2) + SIN2 * HYB2 (3,2)
             B (3,P2  ) = - SIN2 * HYB2 (2,2) + COS2 * HYB2 (3,2)
             B (1,P2+1) =          HYB2 (1,3)
             B (2,P2+1) =   COS2 * HYB2 (2,3) + SIN2 * HYB2 (3,3)
             B (3,P2+1) = - SIN2 * HYB2 (2,3) + COS2 * HYB2 (3,3)

             B (1,S3  ) =          HYB3 (1,1)
             B (2,S3  ) =   COS3 * HYB3 (2,1) + SIN3 * HYB3 (3,1)
             B (3,S3  ) = - SIN3 * HYB3 (2,1) + COS3 * HYB3 (3,1)
             B (1,P3  ) =          HYB3 (1,2)
             B (2,P3  ) =   COS3 * HYB3 (2,2) + SIN3 * HYB3 (3,2)
             B (3,P3  ) = - SIN3 * HYB3 (2,2) + COS3 * HYB3 (3,2)
             B (1,P3+1) =          HYB3 (1,3)
             B (2,P3+1) =   COS3 * HYB3 (2,3) + SIN3 * HYB3 (3,3)
             B (3,P3+1) = - SIN3 * HYB3 (2,3) + COS3 * HYB3 (3,3)

             B (1,S4  ) =          HYB4 (1,1)
             B (2,S4  ) =   COS4 * HYB4 (2,1) + SIN4 * HYB4 (3,1)
             B (3,S4  ) = - SIN4 * HYB4 (2,1) + COS4 * HYB4 (3,1)
             B (1,P4  ) =          HYB4 (1,2)
             B (2,P4  ) =   COS4 * HYB4 (2,2) + SIN4 * HYB4 (3,2)
             B (3,P4  ) = - SIN4 * HYB4 (2,2) + COS4 * HYB4 (3,2)
             B (1,P4+1) =          HYB4 (1,3)
             B (2,P4+1) =   COS4 * HYB4 (2,3) + SIN4 * HYB4 (3,3)
             B (3,P4+1) = - SIN4 * HYB4 (2,3) + COS4 * HYB4 (3,3)

             SYMNBO (S1) = 2
             SYMNBO (S2) = 2
             SYMNBO (S3) = 2
             SYMNBO (S4) = 2

             SLTYPE (S1  ) = -1
             SLTYPE (P1  ) = -1
             SLTYPE (P1+1) = -1
             SLTYPE (S2  ) = -1
             SLTYPE (P2  ) = -1
             SLTYPE (P2+1) = -1
             SLTYPE (S3  ) = -1
             SLTYPE (P3  ) = -1
             SLTYPE (P3+1) = -1
             SLTYPE (S4  ) = -1
             SLTYPE (P4  ) = -1
             SLTYPE (P4+1) = -1

             BDNBAS (BDS1 ) = 3
             BDNBAS (BDP1X) = 3
             BDNBAS (BDP1Y) = 3
             BDNBAS (BDS2 ) = 3
             BDNBAS (BDP2X) = 3
             BDNBAS (BDP2Y) = 3
             BDNBAS (BDS3 ) = 3
             BDNBAS (BDP3X) = 3
             BDNBAS (BDP3Y) = 3
             BDNBAS (BDS4 ) = 3
             BDNBAS (BDP4X) = 3
             BDNBAS (BDP4Y) = 3

             BDBAS (1,BDS1 ) = NHOS1
             BDBAS (2,BDS1 ) = NHOP1X
             BDBAS (3,BDS1 ) = NHOP1Y
             BDBAS (1,BDP1X) = NHOS1
             BDBAS (2,BDP1X) = NHOP1X
             BDBAS (3,BDP1X) = NHOP1Y
             BDBAS (1,BDP1Y) = NHOS1
             BDBAS (2,BDP1Y) = NHOP1X
             BDBAS (3,BDP1Y) = NHOP1Y
             BDBAS (1,BDS2 ) = NHOS2
             BDBAS (2,BDS2 ) = NHOP2X
             BDBAS (3,BDS2 ) = NHOP2Y
             BDBAS (1,BDP2X) = NHOS2
             BDBAS (2,BDP2X) = NHOP2X
             BDBAS (3,BDP2X) = NHOP2Y
             BDBAS (1,BDP2Y) = NHOS2
             BDBAS (2,BDP2Y) = NHOP2X
             BDBAS (3,BDP2Y) = NHOP2Y
             BDBAS (1,BDS3 ) = NHOS3
             BDBAS (2,BDS3 ) = NHOP3X
             BDBAS (3,BDS3 ) = NHOP3Y
             BDBAS (1,BDP3X) = NHOS3
             BDBAS (2,BDP3X) = NHOP3X
             BDBAS (3,BDP3X) = NHOP3Y
             BDBAS (1,BDP3Y) = NHOS3
             BDBAS (2,BDP3Y) = NHOP3X
             BDBAS (3,BDP3Y) = NHOP3Y
             BDBAS (1,BDS4 ) = NHOS4
             BDBAS (2,BDS4 ) = NHOP4X
             BDBAS (3,BDS4 ) = NHOP4Y
             BDBAS (1,BDP4X) = NHOS4
             BDBAS (2,BDP4X) = NHOP4X
             BDBAS (3,BDP4X) = NHOP4Y
             BDBAS (1,BDP4Y) = NHOS4
             BDBAS (2,BDP4Y) = NHOP4X
             BDBAS (3,BDP4Y) = NHOP4Y

         ELSE
C
C
C             ...no sp2 hybridization possible. Print info.
C
C
             WRITE (*,*) ' Cannot perform sp2 hybridization! '
             WRITE (*,*) ' T-symmetry reduction will result! '

         END IF

C
C
C             ...ready!
C
C
         RETURN
         END
