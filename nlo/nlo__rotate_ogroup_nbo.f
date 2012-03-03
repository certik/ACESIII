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
         SUBROUTINE  NLO__ROTATE_OGROUP_NBO
     +
     +                    ( NBAS,NATOM,NSYM,
     +                      NBOSIZE,
     +                      NXB,NXA,NXBA,
     +                      ATXIDX,ATXOFF,
     +                      DSYMACC,LSYMACC,QSYMACC,
     +                      NPAIR,
     +                      IDXSYM,
     +                      C3AXIS,C4AXIS,
     +                      NBOBD,
     +                      P,
     +                      LOCAL,
     +                      W,Q,
     +                      DOMAT,
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
C  OPERATION   : NLO__ROTATE_OGROUP_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine rotates an octahedral set (six or eight)
C                of two-fold interaction order degenerate atomic NBOs,
C                such that they become equivalent under operations of
C                the octahedral groups O and Oh.
C
C                The first step is to divide the set into three (if on
C                C4 axis) or four (if on C3 axis) pairs and rotating
C                each pair of degenerate atomic NBOs such that their
C                lobes match and form rectangular angles between them
C                (for the mathematics on how this is done, please see
C                the corresponding NBO pair rotation routine):
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
C                In the octahedral case such an arrangement of the lobes
C                is well defined, as long as the pair of degenerate
C                NBOs is not sitting opposite to each other. If they
C                would be opposite then there is obviously an infinite
C                amount of such arrangements in space.
C
C                For both cases (C3 and C4 axis) we have to choose
C                the pairs such that they all are on neighboring atomic
C                centers, which can be done by using the 'distance'
C                criterion transmitted via the distance order matrix.
C
C                Consider the C3 case. Here we know that the 8 atomic
C                centers are sitting on the vertices of a cube:
C
C
C
C                                   1---------2
C                                  /|        /|
C                                 / |       / |
C                                /  |      /  |
C                               5---------6   |
C                               |   |     |   |
C                               |   3-----|---4
C                               |  /      |  /
C                               | /       | /
C                               |/        |/
C                               7---------8
C
C
C                Permissible sets of atomic center pairs would be:
C
C                        (1 2) , (3 4) , (5 6) , (7 8)
C
C                                      or
C
C                        (1 5) , (2 6) , (3 7) , (4 8)
C
C                The algorithm used here to establish such a set is
C                as follows:
C
C                    i) Pick any center, say  # 4
C                   ii) Find one nearest neighbor, say # 8
C                  iii) Find all nearest neighbors for the pair (4 8)
C                       -> # 2,3,6,7
C                   iv) Among the set 2,3,6,7 find the nearest neighbor
C                       to # 2 -> # 6. Thus we have the pairs (2 6)
C                       and (3 7)
C                    v) The remaining # 1 and # 5 form the last pair
C                       (1 5)
C
C
C                For the C4 case we know that the 6 atomic centers are
C                sitting on the vertices of an octahedron (looking
C                from above, # 6 hidden):
C
C
C                               1-------------------2
C                               |*        |        *|
C                               |  *      |      *  |
C                               |    *    |    *    |
C                               |      *  |  *      |
C                               |        *|*        |
C                               |---------3---------|
C                               |        *|*        |
C                               |      *  |  *      |
C                               |    *    |    *    |
C                               |  *      |      *  |
C                               |*        |        *|
C                               4-------------------5
C
C
C                A slightly different algorithm must be used to
C                establish the 3 nearest neighbor pairs:
C
C                    i) Pick any center, say # 3
C                   ii) Find one nearest neighbor, say # 5,
C                       thus having the first pair (3 5)
C                  iii) Choose any nearest neighbor of # 5
C                       but excluding the center opposite (farthest
C                       away) to # 3, say # 2
C                   iv) Find any nearest neighbor of # 2, say # 6,
C                       thus having the second pair (2 6)
C                    v) The remaining # 1 and # 5 form the last pair
C                       (1 4)
C      
C
C                For the C4 case the story ends here, as the rotated
C                pairs of degenerate atomic NBOs transform properly
C                under C4. For the C3 case however we need to proceed
C                further and try to construct sp2 hybrids that will
C                transform properly under C3.
C
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    NSYM         =  # of atomic centers related
C                                    by symmetry
C                    NBOSIZE      =  maximum # of NHOs that will
C                                    participate in forming a NBO.
C                    NXB          =  total # of atomic x-type NBOs.
C                    NXA          =  total # of atomic x-type NBO atoms.
C                    NXBA         =  # of atomic x-type NBOs on one of
C                                    the octahedrally related atoms.
C                    ATXIDX (A)   =  atomic index for x-type NBO atom A.
C                    ATXOFF (A)   =  index offset for atomic NBOs for
C                                    x-type NBO atom A. This index is
C                                    equal to the total number of atomic
C                                    x-type NBOs on all x-type NBO atoms
C                                    preceeding x-type NBO atom A.
C                    DSYMACC      =  symmetry accuracy for distance
C                                    order matrix values
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
C                    IDXSYM       =  array that contains the six or
C                                    eight octahedrally related x-type
C                                    NBO atomic indices.
C                    CxAXIS       =  is true, if the symmetry related
C                                    atomic centers are sitting on a
C                                    rotational Cx axis (x=3 or 4)
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    for the I-th atomic x-type NBO.
C                                    This array is the handle for
C                                    accessing and modifying info of
C                                    the NHO bonds sitting in the
C                                    arrays BDNCEN,BDCEN,BDBAS and
C                                    BDOCC.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    LOCAL        =  atomic localization content for
C                                    all NBOs on one of the 6 or 8
C                                    atoms.
C                    W            =  NBO weight vector at present stage
C                                    of symmetrization for all NBOs
C                                    on one of the 6 or 8 atoms.
C                    Q            =  NBO interaction order vector at
C                                    present stage of symmetrization
C                                    for all NBOs on one of the 6 or
C                                    8 atoms.
C                    DOMAT (A,B)  =  distance order matrix containing
C                                    'distance' info between atoms A
C                                    and B.
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

         LOGICAL     ALARM
         LOGICAL     C3AXIS,C4AXIS
         LOGICAL     HYBRID

         INTEGER     ATOM1,ATOM2,ATOM3,ATOM4,ATOM5,ATOM6,ATOM7,ATOM8
         INTEGER     BDS1,BDS2,BDS3,BDS4,BDS5,BDS6,BDS7,BDS8
         INTEGER     BDP1X,BDP2X,BDP3X,BDP4X,BDP5X,BDP6X,BDP7X,BDP8X
         INTEGER     BDP1Y,BDP2Y,BDP3Y,BDP4Y,BDP5Y,BDP6Y,BDP7Y,BDP8Y
         INTEGER     I
         INTEGER     INDEX
         INTEGER     NBAS,NATOM,NSYM
         INTEGER     NBOSIZE
         INTEGER     NHOS1,NHOS2,NHOS3,NHOS4,NHOS5,NHOS6,NHOS7,NHOS8
         INTEGER     NHOP1X,NHOP2X,NHOP3X,NHOP4X,
     +               NHOP5X,NHOP6X,NHOP7X,NHOP8X
         INTEGER     NHOP1Y,NHOP2Y,NHOP3Y,NHOP4Y,
     +               NHOP5Y,NHOP6Y,NHOP7Y,NHOP8Y
         INTEGER     NLEFT
         INTEGER     NLO__ANGULAR_MOMENTUM_NBO
         INTEGER     NPAIR
         INTEGER     NXB,NXA,NXBA
         INTEGER     OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7,OFF8
         INTEGER     P1,P2,P3,P4,P5,P6,P7,P8
         INTEGER     S1,S2,S3,S4,S5,S6,S7,S8
         INTEGER     TOTSYM

         INTEGER     ATXIDX  (1:NXA )
         INTEGER     ATXOFF  (1:NXA )
         INTEGER     BDNBAS  (1:NBAS)
         INTEGER     IDXSYM  (1:NSYM)
         INTEGER     LEFT    (1:4   )
         INTEGER     NBOBD   (1:NXB )
         INTEGER     SLTYPE  (1:NXB )
         INTEGER     SYMNBO  (1:NXB )

         INTEGER     BDBAS   (1:NBOSIZE,1:NBAS)

         DOUBLE PRECISION  COS1,COS2,COS3,COS4,COS5,COS6,COS7,COS8
         DOUBLE PRECISION  D1,D2,D3,D4
         DOUBLE PRECISION  DMAX,DMIN,DIFF
         DOUBLE PRECISION  DSYMACC,LSYMACC,QSYMACC
         DOUBLE PRECISION  SIN1,SIN2,SIN3,SIN4,SIN5,SIN6,SIN7,SIN8
         DOUBLE PRECISION  ZERO,SQR2TH

         DOUBLE PRECISION  LOCAL (1:NXBA)
         DOUBLE PRECISION  Q     (1:NXBA)
         DOUBLE PRECISION  W     (1:NXBA)
         DOUBLE PRECISION  XVEC  (1:NBAS)

         DOUBLE PRECISION  HYB1  (1:3,1:3)
         DOUBLE PRECISION  HYB2  (1:3,1:3)
         DOUBLE PRECISION  HYB3  (1:3,1:3)
         DOUBLE PRECISION  HYB4  (1:3,1:3)
         DOUBLE PRECISION  HYB5  (1:3,1:3)
         DOUBLE PRECISION  HYB6  (1:3,1:3)
         DOUBLE PRECISION  HYB7  (1:3,1:3)
         DOUBLE PRECISION  HYB8  (1:3,1:3)

         DOUBLE PRECISION  B     (1:NBOSIZE,1:NXB  )
         DOUBLE PRECISION  C     (1:NBAS   ,1:NXB  )
         DOUBLE PRECISION  P     (1:NBAS   ,1:NBAS )
         DOUBLE PRECISION  DOMAT (1:NATOM  ,1:NATOM)
         DOUBLE PRECISION  XMAT  (1:NBAS   ,1:2    )

         PARAMETER  (ZERO   = 0.D0)
         PARAMETER  (SQR2TH = 0.70710678118654752440084436210484909D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...decide on the case.
C
C
         IF (C3AXIS) THEN
C
C
C             ...the C3 axis (cubical background) case. Establish the
C                four nearest neighbor atomic center pairs (P1 P2),
C                (P3 P4), (P5 P6), (P7 P8).
C
C
             P1 = IDXSYM (1)
             P2 = 0
             P3 = 0
             P4 = 0
             P5 = 0
             P6 = 0
             P7 = 0
             P8 = 0

             ATOM1 = ATXIDX (P1)

             DMAX = ZERO
             DMIN = DOMAT (ATOM1,ATOM1)
             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                ATOM2 = ATXIDX (INDEX)
                DMAX = MAX (DOMAT (ATOM2,ATOM1),DMAX)
                DMIN = MIN (DOMAT (ATOM2,ATOM1),DMIN)
             END DO

             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                ATOM2 = ATXIDX (INDEX)
                DIFF = ABS (DMAX - DOMAT (ATOM2,ATOM1))
                IF (DIFF.LT.DSYMACC) THEN
                    P2 = INDEX
                END IF
             END DO

             ATOM2 = ATXIDX (P2)

             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                IF (INDEX.NE.P2) THEN
                    ATOM3 = ATXIDX (INDEX)
                    DIFF = ABS (DMIN - DOMAT (ATOM3,ATOM1))
                    IF (DIFF.LT.DSYMACC) THEN
                        P7 = INDEX
                    END IF
                    DIFF = ABS (DMIN - DOMAT (ATOM3,ATOM2))
                    IF (DIFF.LT.DSYMACC) THEN
                        P8 = INDEX
                    END IF
                END IF
             END DO

             ATOM7 = ATXIDX (P7)
             ATOM8 = ATXIDX (P8)

             NLEFT = 0
             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                IF (INDEX.NE.P2.AND.INDEX.NE.P7.AND.INDEX.NE.P8) THEN
                    NLEFT = NLEFT + 1
                    LEFT (NLEFT) = INDEX
                END IF
             END DO

             P3 = LEFT (1)
             ATOM3 = ATXIDX (P3)

             DO I = 2,NLEFT
                INDEX = LEFT (I)
                ATOM4 = ATXIDX (INDEX)
                DIFF = ABS (DMAX - DOMAT (ATOM4,ATOM3))
                IF (DIFF.LT.DSYMACC) THEN
                    P4 = INDEX
                END IF
             END DO

             ATOM4 = ATXIDX (P4)

             DO I = 2,NLEFT
                INDEX = LEFT (I)
                IF (INDEX.NE.P4) THEN
                    ATOM5 = ATXIDX (INDEX)
                    DIFF = ABS (DMIN - DOMAT (ATOM5,ATOM3))
                    IF (DIFF.LT.DSYMACC) THEN
                        P5 = INDEX
                    END IF
                    DIFF = ABS (DMIN - DOMAT (ATOM5,ATOM4))
                    IF (DIFF.LT.DSYMACC) THEN
                        P6 = INDEX
                    END IF
                END IF
             END DO

             ATOM5 = ATXIDX (P5)
             ATOM6 = ATXIDX (P6)

             D1 = DOMAT (ATOM1,ATOM2)
             D2 = DOMAT (ATOM3,ATOM4)
             D3 = DOMAT (ATOM5,ATOM6)
             D4 = DOMAT (ATOM7,ATOM8)

             ALARM = (P1.EQ.0) .OR. (P2.EQ.0) .OR.
     +               (P3.EQ.0) .OR. (P4.EQ.0) .OR.
     +               (P5.EQ.0) .OR. (P6.EQ.0) .OR.
     +               (P7.EQ.0) .OR. (P8.EQ.0) .OR.
     +               (ABS (DMAX - D1).GT.DSYMACC) .OR.
     +               (ABS (DMAX - D2).GT.DSYMACC) .OR.
     +               (ABS (DMAX - D3).GT.DSYMACC) .OR.
     +               (ABS (DMAX - D4).GT.DSYMACC)

             IF (ALARM) THEN
                 WRITE (*,*) ' Problems finding octahedral C3 pairs! '
                 WRITE (*,*) ' P1,P2,P3,P4,P5,P6,P7,P8 = ',
     +                         P1,P2,P3,P4,P5,P6,P7,P8
                 WRITE (*,*) ' Overlap orders #1,2,3,4,max = ',
     +                         D1,D2,D3,D4,DMAX
                 WRITE (*,*) ' nlo__rotate_ogroup_nbo '
                 WRITE (1,*) ' Problems finding octahedral C3 pairs! '
                 WRITE (1,*) ' P1,P2,P3,P4,P5,P6,P7,P8 = ',
     +                         P1,P2,P3,P4,P5,P6,P7,P8
                 WRITE (1,*) ' Overlap orders #1,2,3,4,max = ',
     +                         D1,D2,D3,D4,DMAX
                 WRITE (1,*) ' nlo__rotate_ogroup_nbo '
                 STOP
             END IF
C
C
C             ...rotate the four pairs into appropriate positions
C                by calling four times the pair rotation routine.
C
C
             OFF1 = ATXOFF (P1)
             OFF2 = ATXOFF (P2)
             OFF3 = ATXOFF (P3)
             OFF4 = ATXOFF (P4)
             OFF5 = ATXOFF (P5)
             OFF6 = ATXOFF (P6)
             OFF7 = ATXOFF (P7)
             OFF8 = ATXOFF (P8)

             P1 = OFF1 + NPAIR
             P2 = OFF2 + NPAIR
             P3 = OFF3 + NPAIR
             P4 = OFF4 + NPAIR
             P5 = OFF5 + NPAIR
             P6 = OFF6 + NPAIR
             P7 = OFF7 + NPAIR
             P8 = OFF8 + NPAIR

             BDP1X = NBOBD (P1  )
             BDP1Y = NBOBD (P1+1)
             BDP2X = NBOBD (P2  )
             BDP2Y = NBOBD (P2+1)
             BDP3X = NBOBD (P3  )
             BDP3Y = NBOBD (P3+1)
             BDP4X = NBOBD (P4  )
             BDP4Y = NBOBD (P4+1)
             BDP5X = NBOBD (P5  )
             BDP5Y = NBOBD (P5+1)
             BDP6X = NBOBD (P6  )
             BDP6Y = NBOBD (P6+1)
             BDP7X = NBOBD (P7  )
             BDP7Y = NBOBD (P7+1)
             BDP8X = NBOBD (P8  )
             BDP8Y = NBOBD (P8+1)

             NHOP1X = BDBAS (1,BDP1X)
             NHOP1Y = BDBAS (1,BDP1Y)
             NHOP2X = BDBAS (1,BDP2X)
             NHOP2Y = BDBAS (1,BDP2Y)
             NHOP3X = BDBAS (1,BDP3X)
             NHOP3Y = BDBAS (1,BDP3Y)
             NHOP4X = BDBAS (1,BDP4X)
             NHOP4Y = BDBAS (1,BDP4Y)
             NHOP5X = BDBAS (1,BDP5X)
             NHOP5Y = BDBAS (1,BDP5Y)
             NHOP6X = BDBAS (1,BDP6X)
             NHOP6Y = BDBAS (1,BDP6Y)
             NHOP7X = BDBAS (1,BDP7X)
             NHOP7Y = BDBAS (1,BDP7Y)
             NHOP8X = BDBAS (1,BDP8X)
             NHOP8Y = BDBAS (1,BDP8Y)

             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS1,SIN1,
     +                           COS2,SIN2,
     +                           C (1,P1),
     +                           C (1,P2) )
     +
     +
             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS3,SIN3,
     +                           COS4,SIN4,
     +                           C (1,P3),
     +                           C (1,P4) )
     +
     +
             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS5,SIN5,
     +                           COS6,SIN6,
     +                           C (1,P5),
     +                           C (1,P6) )
     +
     +
             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS7,SIN7,
     +                           COS8,SIN8,
     +                           C (1,P7),
     +                           C (1,P8) )
     +
     +
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
     +                     ( NXBA,
     +                       SLTYPE (OFF1+1),
     +                       NPAIR,NPAIR+1,
     +                       0,
     +                       1,
     +                       QSYMACC,LSYMACC,
     +                       W,Q,
     +                       LOCAL )
     +
     +
             HYBRID = TOTSYM .NE. 0

             IF (HYBRID) THEN

                 IF (NBOSIZE.LT.3) THEN
                     WRITE (*,*) ' Cannot form sp2 hybrids! '
                     WRITE (*,*) ' NBOSIZE = ',NBOSIZE,' (> 2 required)'
                     WRITE (*,*) ' nlo__rotate_ogroup_nbo '
                     WRITE (1,*) ' Cannot form sp2 hybrids! '
                     WRITE (1,*) ' NBOSIZE = ',NBOSIZE,' (> 2 required)'
                     WRITE (1,*) ' nlo__rotate_ogroup_nbo '
                     STOP
                 END IF

                 S1 = OFF1 + TOTSYM
                 S2 = OFF2 + TOTSYM
                 S3 = OFF3 + TOTSYM
                 S4 = OFF4 + TOTSYM
                 S5 = OFF5 + TOTSYM
                 S6 = OFF6 + TOTSYM
                 S7 = OFF7 + TOTSYM
                 S8 = OFF8 + TOTSYM

                 BDS1 = NBOBD (S1)
                 BDS2 = NBOBD (S2)
                 BDS3 = NBOBD (S3)
                 BDS4 = NBOBD (S4)
                 BDS5 = NBOBD (S5)
                 BDS6 = NBOBD (S6)
                 BDS7 = NBOBD (S7)
                 BDS8 = NBOBD (S8)

                 NHOS1 = BDBAS (1,BDS1)
                 NHOS2 = BDBAS (1,BDS2)
                 NHOS3 = BDBAS (1,BDS3)
                 NHOS4 = BDBAS (1,BDS4)
                 NHOS5 = BDBAS (1,BDS5)
                 NHOS6 = BDBAS (1,BDS6)
                 NHOS7 = BDBAS (1,BDS7)
                 NHOS8 = BDBAS (1,BDS8)

                 CALL  NLO__FORM_SP2_PAIR_NBO
     +
     +                      ( NBAS,
     +                        P,
     +                        XVEC,
     +
     +                               HYB1,HYB2,
     +                               C (1,S1),
     +                               C (1,P1),
     +                               C (1,S2),
     +                               C (1,P2) )
     +
     +
                 CALL  NLO__FORM_SP2_PAIR_NBO
     +
     +                      ( NBAS,
     +                        P,
     +                        XVEC,
     +
     +                               HYB3,HYB4,
     +                               C (1,S3),
     +                               C (1,P3),
     +                               C (1,S4),
     +                               C (1,P4) )
     +
     +
                 CALL  NLO__FORM_SP2_PAIR_NBO
     +
     +                      ( NBAS,
     +                        P,
     +                        XVEC,
     +
     +                               HYB5,HYB6,
     +                               C (1,S5),
     +                               C (1,P5),
     +                               C (1,S6),
     +                               C (1,P6) )
     +
     +
                 CALL  NLO__FORM_SP2_PAIR_NBO
     +
     +                      ( NBAS,
     +                        P,
     +                        XVEC,
     +
     +                               HYB7,HYB8,
     +                               C (1,S7),
     +                               C (1,P7),
     +                               C (1,S8),
     +                               C (1,P8) )
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

                 B (1,S5  ) =          HYB5 (1,1)
                 B (2,S5  ) =   COS5 * HYB5 (2,1) + SIN5 * HYB5 (3,1)
                 B (3,S5  ) = - SIN5 * HYB5 (2,1) + COS5 * HYB5 (3,1)
                 B (1,P5  ) =          HYB5 (1,2)
                 B (2,P5  ) =   COS5 * HYB5 (2,2) + SIN5 * HYB5 (3,2)
                 B (3,P5  ) = - SIN5 * HYB5 (2,2) + COS5 * HYB5 (3,2)
                 B (1,P5+1) =          HYB5 (1,3)
                 B (2,P5+1) =   COS5 * HYB5 (2,3) + SIN5 * HYB5 (3,3)
                 B (3,P5+1) = - SIN5 * HYB5 (2,3) + COS5 * HYB5 (3,3)

                 B (1,S6  ) =          HYB6 (1,1)
                 B (2,S6  ) =   COS6 * HYB6 (2,1) + SIN6 * HYB6 (3,1)
                 B (3,S6  ) = - SIN6 * HYB6 (2,1) + COS6 * HYB6 (3,1)
                 B (1,P6  ) =          HYB6 (1,2)
                 B (2,P6  ) =   COS6 * HYB6 (2,2) + SIN6 * HYB6 (3,2)
                 B (3,P6  ) = - SIN6 * HYB6 (2,2) + COS6 * HYB6 (3,2)
                 B (1,P6+1) =          HYB6 (1,3)
                 B (2,P6+1) =   COS6 * HYB6 (2,3) + SIN6 * HYB6 (3,3)
                 B (3,P6+1) = - SIN6 * HYB6 (2,3) + COS6 * HYB6 (3,3)

                 B (1,S7  ) =          HYB7 (1,1)
                 B (2,S7  ) =   COS7 * HYB7 (2,1) + SIN7 * HYB7 (3,1)
                 B (3,S7  ) = - SIN7 * HYB7 (2,1) + COS7 * HYB7 (3,1)
                 B (1,P7  ) =          HYB7 (1,2)
                 B (2,P7  ) =   COS7 * HYB7 (2,2) + SIN7 * HYB7 (3,2)
                 B (3,P7  ) = - SIN7 * HYB7 (2,2) + COS7 * HYB7 (3,2)
                 B (1,P7+1) =          HYB7 (1,3)
                 B (2,P7+1) =   COS7 * HYB7 (2,3) + SIN7 * HYB7 (3,3)
                 B (3,P7+1) = - SIN7 * HYB7 (2,3) + COS7 * HYB7 (3,3)

                 B (1,S8  ) =          HYB8 (1,1)
                 B (2,S8  ) =   COS8 * HYB8 (2,1) + SIN8 * HYB8 (3,1)
                 B (3,S8  ) = - SIN8 * HYB8 (2,1) + COS8 * HYB8 (3,1)
                 B (1,P8  ) =          HYB8 (1,2)
                 B (2,P8  ) =   COS8 * HYB8 (2,2) + SIN8 * HYB8 (3,2)
                 B (3,P8  ) = - SIN8 * HYB8 (2,2) + COS8 * HYB8 (3,2)
                 B (1,P8+1) =          HYB8 (1,3)
                 B (2,P8+1) =   COS8 * HYB8 (2,3) + SIN8 * HYB8 (3,3)
                 B (3,P8+1) = - SIN8 * HYB8 (2,3) + COS8 * HYB8 (3,3)

                 SYMNBO (S1)   = 2
                 SYMNBO (P1  ) = 2
                 SYMNBO (P1+1) = 2
                 SYMNBO (S2)   = 2
                 SYMNBO (P2  ) = 2
                 SYMNBO (P2+1) = 2
                 SYMNBO (S3)   = 2
                 SYMNBO (P3  ) = 2
                 SYMNBO (P3+1) = 2
                 SYMNBO (S4)   = 2
                 SYMNBO (P4  ) = 2
                 SYMNBO (P4+1) = 2
                 SYMNBO (S5)   = 2
                 SYMNBO (P5  ) = 2
                 SYMNBO (P5+1) = 2
                 SYMNBO (S6)   = 2
                 SYMNBO (P6  ) = 2
                 SYMNBO (P6+1) = 2
                 SYMNBO (S7)   = 2
                 SYMNBO (P7  ) = 2
                 SYMNBO (P7+1) = 2
                 SYMNBO (S8)   = 2
                 SYMNBO (P8  ) = 2
                 SYMNBO (P8+1) = 2

                 SLTYPE (S1)   = -1
                 SLTYPE (P1  ) = -1
                 SLTYPE (P1+1) = -1
                 SLTYPE (S2)   = -1
                 SLTYPE (P2  ) = -1
                 SLTYPE (P2+1) = -1
                 SLTYPE (S3)   = -1
                 SLTYPE (P3  ) = -1
                 SLTYPE (P3+1) = -1
                 SLTYPE (S4)   = -1
                 SLTYPE (P4  ) = -1
                 SLTYPE (P4+1) = -1
                 SLTYPE (S5)   = -1
                 SLTYPE (P5  ) = -1
                 SLTYPE (P5+1) = -1
                 SLTYPE (S6)   = -1
                 SLTYPE (P6  ) = -1
                 SLTYPE (P6+1) = -1
                 SLTYPE (S7)   = -1
                 SLTYPE (P7  ) = -1
                 SLTYPE (P7+1) = -1
                 SLTYPE (S8)   = -1
                 SLTYPE (P8  ) = -1
                 SLTYPE (P8+1) = -1

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
                 BDNBAS (BDS5 ) = 3
                 BDNBAS (BDP5X) = 3
                 BDNBAS (BDP5Y) = 3
                 BDNBAS (BDS6 ) = 3
                 BDNBAS (BDP6X) = 3
                 BDNBAS (BDP6Y) = 3
                 BDNBAS (BDS7 ) = 3
                 BDNBAS (BDP7X) = 3
                 BDNBAS (BDP7Y) = 3
                 BDNBAS (BDS8 ) = 3
                 BDNBAS (BDP8X) = 3
                 BDNBAS (BDP8Y) = 3

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
                 BDBAS (1,BDS5 ) = NHOS5
                 BDBAS (2,BDS5 ) = NHOP5X
                 BDBAS (3,BDS5 ) = NHOP5Y
                 BDBAS (1,BDP5X) = NHOS5
                 BDBAS (2,BDP5X) = NHOP5X
                 BDBAS (3,BDP5X) = NHOP5Y
                 BDBAS (1,BDP5Y) = NHOS5
                 BDBAS (2,BDP5Y) = NHOP5X
                 BDBAS (3,BDP5Y) = NHOP5Y
                 BDBAS (1,BDS6 ) = NHOS6
                 BDBAS (2,BDS6 ) = NHOP6X
                 BDBAS (3,BDS6 ) = NHOP6Y
                 BDBAS (1,BDP6X) = NHOS6
                 BDBAS (2,BDP6X) = NHOP6X
                 BDBAS (3,BDP6X) = NHOP6Y
                 BDBAS (1,BDP6Y) = NHOS6
                 BDBAS (2,BDP6Y) = NHOP6X
                 BDBAS (3,BDP6Y) = NHOP6Y
                 BDBAS (1,BDS7 ) = NHOS7
                 BDBAS (2,BDS7 ) = NHOP7X
                 BDBAS (3,BDS7 ) = NHOP7Y
                 BDBAS (1,BDP7X) = NHOS7
                 BDBAS (2,BDP7X) = NHOP7X
                 BDBAS (3,BDP7X) = NHOP7Y
                 BDBAS (1,BDP7Y) = NHOS7
                 BDBAS (2,BDP7Y) = NHOP7X
                 BDBAS (3,BDP7Y) = NHOP7Y
                 BDBAS (1,BDS8 ) = NHOS8
                 BDBAS (2,BDS8 ) = NHOP8X
                 BDBAS (3,BDS8 ) = NHOP8Y
                 BDBAS (1,BDP8X) = NHOS8
                 BDBAS (2,BDP8X) = NHOP8X
                 BDBAS (3,BDP8X) = NHOP8Y
                 BDBAS (1,BDP8Y) = NHOS8
                 BDBAS (2,BDP8Y) = NHOP8X
                 BDBAS (3,BDP8Y) = NHOP8Y

             ELSE
C
C
C             ...no sp2 hybridization possible. Print info.
C
C
                 WRITE (*,*) ' Cannot perform sp2 hybridization! '
                 WRITE (*,*) ' O-symmetry reduction will result! '

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
                 B (1,P5  ) =   COS5
                 B (2,P5  ) = - SIN5
                 B (1,P5+1) =   SIN5
                 B (2,P5+1) =   COS5
                 B (1,P6  ) =   COS6
                 B (2,P6  ) = - SIN6
                 B (1,P6+1) =   SIN6
                 B (2,P6+1) =   COS6
                 B (1,P7  ) =   COS7
                 B (2,P7  ) = - SIN7
                 B (1,P7+1) =   SIN7
                 B (2,P7+1) =   COS7
                 B (1,P8  ) =   COS8
                 B (2,P8  ) = - SIN8
                 B (1,P8+1) =   SIN8
                 B (2,P8+1) =   COS8

                 SYMNBO (P1  ) = 2
                 SYMNBO (P1+1) = 2
                 SYMNBO (P2  ) = 2
                 SYMNBO (P2+1) = 2
                 SYMNBO (P3  ) = 2
                 SYMNBO (P3+1) = 2
                 SYMNBO (P4  ) = 2
                 SYMNBO (P4+1) = 2
                 SYMNBO (P5  ) = 2
                 SYMNBO (P5+1) = 2
                 SYMNBO (P6  ) = 2
                 SYMNBO (P6+1) = 2
                 SYMNBO (P7  ) = 2
                 SYMNBO (P7+1) = 2
                 SYMNBO (P8  ) = 2
                 SYMNBO (P8+1) = 2

                 BDNBAS (BDP1X) = 2
                 BDNBAS (BDP1Y) = 2
                 BDNBAS (BDP2X) = 2
                 BDNBAS (BDP2Y) = 2
                 BDNBAS (BDP3X) = 2
                 BDNBAS (BDP3Y) = 2
                 BDNBAS (BDP4X) = 2
                 BDNBAS (BDP4Y) = 2
                 BDNBAS (BDP5X) = 2
                 BDNBAS (BDP5Y) = 2
                 BDNBAS (BDP6X) = 2
                 BDNBAS (BDP6Y) = 2
                 BDNBAS (BDP7X) = 2
                 BDNBAS (BDP7Y) = 2
                 BDNBAS (BDP8X) = 2
                 BDNBAS (BDP8Y) = 2

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
                 BDBAS (1,BDP5X) = NHOP5X
                 BDBAS (2,BDP5X) = NHOP5Y
                 BDBAS (1,BDP5Y) = NHOP5X
                 BDBAS (2,BDP5Y) = NHOP5Y
                 BDBAS (1,BDP6X) = NHOP6X
                 BDBAS (2,BDP6X) = NHOP6Y
                 BDBAS (1,BDP6Y) = NHOP6X
                 BDBAS (2,BDP6Y) = NHOP6Y
                 BDBAS (1,BDP7X) = NHOP7X
                 BDBAS (2,BDP7X) = NHOP7Y
                 BDBAS (1,BDP7Y) = NHOP7X
                 BDBAS (2,BDP7Y) = NHOP7Y
                 BDBAS (1,BDP8X) = NHOP8X
                 BDBAS (2,BDP8X) = NHOP8Y
                 BDBAS (1,BDP8Y) = NHOP8X
                 BDBAS (2,BDP8Y) = NHOP8Y

             END IF

         ELSE IF (C4AXIS) THEN
C
C
C             ...the C4 axis case. Establish the three nearest neighbor
C                atomic center pairs (P1 P2), (P3 P4), (P5 P6).
C
C
             P1 = IDXSYM (1)
             P2 = 0
             P3 = 0
             P4 = 0
             P5 = 0
             P6 = 0

             ATOM1 = ATXIDX (P1)

             DMAX = ZERO
             DMIN = DOMAT (ATOM1,ATOM1)
             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                ATOM2 = ATXIDX (INDEX)
                DMAX = MAX (DOMAT (ATOM2,ATOM1),DMAX)
                DMIN = MIN (DOMAT (ATOM2,ATOM1),DMIN)
             END DO

             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                ATOM2 = ATXIDX (INDEX)
                DIFF = ABS (DMAX - DOMAT (ATOM2,ATOM1))
                IF (DIFF.LT.DSYMACC) THEN
                    P2 = INDEX
                END IF
                DIFF = ABS (DMIN - DOMAT (ATOM2,ATOM1))
                IF (DIFF.LT.DSYMACC) THEN
                    ATOM4 = ATOM2
                END IF
             END DO

             ATOM2 = ATXIDX (P2)

             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                ATOM3 = ATXIDX (INDEX)
                IF (INDEX.NE.P2 .AND. ATOM3.NE.ATOM4) THEN
                    DIFF = ABS (DMAX - DOMAT (ATOM3,ATOM2))
                    IF (DIFF.LT.DSYMACC) THEN
                        P3 = INDEX
                    END IF
                END IF
             END DO

             ATOM3 = ATXIDX (P3)

             NLEFT = 0
             DO I = 2,NSYM
                INDEX = IDXSYM (I)
                IF (INDEX.NE.P2.AND.INDEX.NE.P3) THEN
                    NLEFT = NLEFT + 1
                    LEFT (NLEFT) = INDEX
                END IF
             END DO

             DO I = 1,NLEFT
                INDEX = LEFT (I)
                ATOM4 = ATXIDX (INDEX)
                DIFF = ABS (DMAX - DOMAT (ATOM4,ATOM3))
                IF (DIFF.LT.DSYMACC) THEN
                    P4 = INDEX
                END IF
             END DO

             DO I = 1,NLEFT
                INDEX = LEFT (I)
                IF (INDEX.NE.P4) THEN
                    P5 = INDEX
                END IF
             END DO

             DO I = 1,NLEFT
                INDEX = LEFT (I)
                IF (INDEX.NE.P4.AND.INDEX.NE.P5) THEN
                    P6 = INDEX
                END IF
             END DO

             ATOM4 = ATXIDX (P4)
             ATOM5 = ATXIDX (P5)
             ATOM6 = ATXIDX (P6)

             D1 = DOMAT (ATOM1,ATOM2)
             D2 = DOMAT (ATOM3,ATOM4)
             D3 = DOMAT (ATOM5,ATOM6)

             ALARM = (P1.EQ.0) .OR. (P2.EQ.0) .OR.
     +               (P3.EQ.0) .OR. (P4.EQ.0) .OR.
     +               (P5.EQ.0) .OR. (P6.EQ.0) .OR.
     +               (ABS (DMAX - D1).GT.DSYMACC) .OR.
     +               (ABS (DMAX - D2).GT.DSYMACC) .OR.
     +               (ABS (DMAX - D3).GT.DSYMACC)

             IF (ALARM) THEN
                 WRITE (*,*) ' Problems finding octahedral C4 pairs! '
                 WRITE (*,*) ' P1,P2,P3,P4,P5,P6 = ',
     +                         P1,P2,P3,P4,P5,P6
                 WRITE (*,*) ' Overlap orders #1,2,3,max = ',
     +                         D1,D2,D3,DMAX
                 WRITE (*,*) ' nlo__rotate_ogroup_nbo '
                 WRITE (1,*) ' Problems finding octahedral C4 pairs! '
                 WRITE (1,*) ' P1,P2,P3,P4,P5,P6 = ',
     +                         P1,P2,P3,P4,P5,P6
                 WRITE (1,*) ' Overlap orders #1,2,3,max = ',
     +                         D1,D2,D3,DMAX
                 WRITE (1,*) ' nlo__rotate_ogroup_nbo '
                 STOP
             END IF
C
C
C             ...rotate the three pairs into appropriate positions
C                by calling three times the pair rotation routine.
C
C
             OFF1 = ATXOFF (P1)
             OFF2 = ATXOFF (P2)
             OFF3 = ATXOFF (P3)
             OFF4 = ATXOFF (P4)
             OFF5 = ATXOFF (P5)
             OFF6 = ATXOFF (P6)

             P1 = OFF1 + NPAIR
             P2 = OFF2 + NPAIR
             P3 = OFF3 + NPAIR
             P4 = OFF4 + NPAIR
             P5 = OFF5 + NPAIR
             P6 = OFF6 + NPAIR

             BDP1X = NBOBD (P1  )
             BDP1Y = NBOBD (P1+1)
             BDP2X = NBOBD (P2  )
             BDP2Y = NBOBD (P2+1)
             BDP3X = NBOBD (P3  )
             BDP3Y = NBOBD (P3+1)
             BDP4X = NBOBD (P4  )
             BDP4Y = NBOBD (P4+1)
             BDP5X = NBOBD (P5  )
             BDP5Y = NBOBD (P5+1)
             BDP6X = NBOBD (P6  )
             BDP6Y = NBOBD (P6+1)

             NHOP1X = BDBAS (1,BDP1X)
             NHOP1Y = BDBAS (1,BDP1Y)
             NHOP2X = BDBAS (1,BDP2X)
             NHOP2Y = BDBAS (1,BDP2Y)
             NHOP3X = BDBAS (1,BDP3X)
             NHOP3Y = BDBAS (1,BDP3Y)
             NHOP4X = BDBAS (1,BDP4X)
             NHOP4Y = BDBAS (1,BDP4Y)
             NHOP5X = BDBAS (1,BDP5X)
             NHOP5Y = BDBAS (1,BDP5Y)
             NHOP6X = BDBAS (1,BDP6X)
             NHOP6Y = BDBAS (1,BDP6Y)

             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS1,SIN1,
     +                           COS2,SIN2,
     +                           C (1,P1),
     +                           C (1,P2) )
     +
     +
             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS3,SIN3,
     +                           COS4,SIN4,
     +                           C (1,P3),
     +                           C (1,P4) )
     +
     +
             CALL  NLO__ROTATE_PAIR_DEGENERATE_NBO
     +
     +                  ( NBAS,
     +                    P,
     +                    XVEC,
     +                    XMAT,
     +
     +                           COS5,SIN5,
     +                           COS6,SIN6,
     +                           C (1,P5),
     +                           C (1,P6) )
     +
     +
C
C
C             ...rotate the three pairs 45 degrees (not needed!).
C
C
C             DO N = 1,NBAS
C                X = SQR2TH * C (N,P1) - SQR2TH * C (N,P1+1)
C                Y = SQR2TH * C (N,P1) + SQR2TH * C (N,P1+1)
C                C (N,P1  ) = X
C                C (N,P1+1) = Y
C                X = SQR2TH * C (N,P2) - SQR2TH * C (N,P2+1)
C                Y = SQR2TH * C (N,P2) + SQR2TH * C (N,P2+1)
C                C (N,P2  ) = X
C                C (N,P2+1) = Y
C                X = SQR2TH * C (N,P3) - SQR2TH * C (N,P3+1)
C                Y = SQR2TH * C (N,P3) + SQR2TH * C (N,P3+1)
C                C (N,P3  ) = X
C                C (N,P3+1) = Y
C                X = SQR2TH * C (N,P4) - SQR2TH * C (N,P4+1)
C                Y = SQR2TH * C (N,P4) + SQR2TH * C (N,P4+1)
C                C (N,P4  ) = X
C                C (N,P4+1) = Y
C                X = SQR2TH * C (N,P5) - SQR2TH * C (N,P5+1)
C                Y = SQR2TH * C (N,P5) + SQR2TH * C (N,P5+1)
C                C (N,P5  ) = X
C                C (N,P5+1) = Y
C                X = SQR2TH * C (N,P6) - SQR2TH * C (N,P6+1)
C                Y = SQR2TH * C (N,P6) + SQR2TH * C (N,P6+1)
C                C (N,P6  ) = X
C                C (N,P6+1) = Y
C             END DO
C
C             B (1,P1  ) =   SQR2TH * (COS1 - SIN1)
C             B (2,P1  ) = - SQR2TH * (SIN1 + COS1)
C             B (1,P1+1) =   SQR2TH * (COS1 + SIN1)
C             B (2,P1+1) =   SQR2TH * (COS1 - SIN1)
C
C             B (1,P2  ) =   SQR2TH * (COS2 - SIN2)
C             B (2,P2  ) = - SQR2TH * (SIN2 + COS2)
C             B (1,P2+1) =   SQR2TH * (COS2 + SIN2)
C             B (2,P2+1) =   SQR2TH * (COS2 - SIN2)
C
C             B (1,P3  ) =   SQR2TH * (COS3 - SIN3)
C             B (2,P3  ) = - SQR2TH * (SIN3 + COS3)
C             B (1,P3+1) =   SQR2TH * (COS3 + SIN3)
C             B (2,P3+1) =   SQR2TH * (COS3 - SIN3)
C
C             B (1,P4  ) =   SQR2TH * (COS4 - SIN4)
C             B (2,P4  ) = - SQR2TH * (SIN4 + COS4)
C             B (1,P4+1) =   SQR2TH * (COS4 + SIN4)
C             B (2,P4+1) =   SQR2TH * (COS4 - SIN4)
C
C             B (1,P5  ) =   SQR2TH * (COS5 - SIN5)
C             B (2,P5  ) = - SQR2TH * (SIN5 + COS5)
C             B (1,P5+1) =   SQR2TH * (COS5 + SIN5)
C             B (2,P5+1) =   SQR2TH * (COS5 - SIN5)
C
C             B (1,P6  ) =   SQR2TH * (COS6 - SIN6)
C             B (2,P6  ) = - SQR2TH * (SIN6 + COS6)
C             B (1,P6+1) =   SQR2TH * (COS6 + SIN6)
C             B (2,P6+1) =   SQR2TH * (COS6 - SIN6)

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
             B (1,P5  ) =   COS5
             B (2,P5  ) = - SIN5
             B (1,P5+1) =   SIN5
             B (2,P5+1) =   COS5
             B (1,P6  ) =   COS6
             B (2,P6  ) = - SIN6
             B (1,P6+1) =   SIN6
             B (2,P6+1) =   COS6

             SYMNBO (P1  ) = 2
             SYMNBO (P1+1) = 2
             SYMNBO (P2  ) = 2
             SYMNBO (P2+1) = 2
             SYMNBO (P3  ) = 2
             SYMNBO (P3+1) = 2
             SYMNBO (P4  ) = 2
             SYMNBO (P4+1) = 2
             SYMNBO (P5  ) = 2
             SYMNBO (P5+1) = 2
             SYMNBO (P6  ) = 2
             SYMNBO (P6+1) = 2

             BDNBAS (BDP1X) = 2
             BDNBAS (BDP1Y) = 2
             BDNBAS (BDP2X) = 2
             BDNBAS (BDP2Y) = 2
             BDNBAS (BDP3X) = 2
             BDNBAS (BDP3Y) = 2
             BDNBAS (BDP4X) = 2
             BDNBAS (BDP4Y) = 2
             BDNBAS (BDP5X) = 2
             BDNBAS (BDP5Y) = 2
             BDNBAS (BDP6X) = 2
             BDNBAS (BDP6Y) = 2

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
             BDBAS (1,BDP5X) = NHOP5X
             BDBAS (2,BDP5X) = NHOP5Y
             BDBAS (1,BDP5Y) = NHOP5X
             BDBAS (2,BDP5Y) = NHOP5Y
             BDBAS (1,BDP6X) = NHOP6X
             BDBAS (2,BDP6X) = NHOP6Y
             BDBAS (1,BDP6Y) = NHOP6X
             BDBAS (2,BDP6Y) = NHOP6Y

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
