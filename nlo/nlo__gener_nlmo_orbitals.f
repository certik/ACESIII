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
         SUBROUTINE  NLO__GENER_NLMO_ORBITALS
     +
     +                    ( NBAS,NOCC,
     +                      MAXOCC,
     +                      USED,
     +                      PSYMACC,QSYMACC,
     +                      P,
     +                      XMAT,
     +
     +                             W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NLMO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the natural localized molecular
C                orbitals (NLMOs) from the previously found set of
C                NBOs.
C
C                Procedure:
C
C                The NLMO expansion coefficients C in terms of NBOs
C                are found by carefully zeroing the offdiagonal NBO
C                occupation matrix block P(AB,NBO), where A represents
C                the high occupancy block (core, lone-pair, bond NBOs)
C                and B the low occupancy block (anti-bond, Rydberg
C                NBOs). The NLMO expansion coefficients C are thus
C                defined as:
C
C
C
C    |            |            |          |           |           |
C    | P(AA,NLMO) |      0     |          | P(AA,NBO) | P(AB,NBO) |
C    |            |            |          |           |           |
C    |------------|------------|  =  C(T) |-----------|-----------| C
C    |            |            |          |           |           |
C    |     0      | P(BB,NLMO) |          | P(BA,NBO) | P(BB,NBO) |
C    |            |            |          |           |           |
C
C
C
C                Note that if the one-electron occupation matrix used
C                for finding the NAO and NBO sets corresponds to the
C                closed shell SCF density, then P(AA,NLMO) will be
C                equal to the diagonal matrix 2I (with I = unit matrix)
C                and P(BB,NLMO) will become a zero matrix. Only if
C                we use originally a correlated one-electron occupation
C                matrix we will have non-diagonal P(AA,NLMO) and
C                P(BB,NLMO) blocks.
C
C                To find the NLMO expansion coefficient matrix, a
C                sequence of well controlled symmetrized Jacobi
C                rotations is performed on P(NBO):
C
C                   1) Find largest element Pij of P(AB,NBO)
C                      in magnitude. If |Pij| below certain
C                      threshold => stop.
C
C                   2) Find all elements Pkl of P(AB,NBO) for
C                      which |Pkl| = |Pij|, Pkk = Pii and Pll = Pjj
C                      (within certain accuracy limits). All these
C                      2x2 subblocks of P(AB,NBO) are related by
C                      symmetry and should be treated all at once
C                      in order not to break the symmetry.
C
C                   3) Perform a symmetrized Jacobi rotation on
C                      all the 2x2 subblocks found in 2).
C
C                   4) Return to 1).
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis.
C                    NOCC         =  # of occupied (large weight) NBOs.
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    USED         =  int vector that will be used
C                                    for the symmetrized Jacobi routine.
C                    PSYMACC      =  symmetry accuracy for occupation
C                                    matrix values (including weights)
C                    QSYMACC      =  symmetry accuracy for orbital
C                                    interaction order values (sum of
C                                    squares of offdiagonal occupation
C                                    matrix elements for one orbital)
C                    P            =  full NBAS x NBAS occupation matrix
C                                    in AO basis.
C                    XMAT         =  flp scratch array of matrix type.
C                    C            =  NBO coefficient matrix in AO basis.
C
C
C                  Output:
C
C                    W            =  NLMO weight vector in original NBO
C                                    order.
C                    C            =  NLMO coefficient matrix in AO basis
C                                    with columns in original NBO order.
C
C
C                  Note that the occupation matrix P in AO basis will
C                  be destroyed here.
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

         LOGICAL     LTRG,UTRG
         LOGICAL     REVERS
         LOGICAL     SAVEP

         INTEGER     NBAS
         INTEGER     NITER
         INTEGER     NOCC

         INTEGER     USED  (1:NBAS)

         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  PSYMACC,QSYMACC
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  W (1:NBAS)

         DOUBLE PRECISION  C     (1:NBAS,1:NBAS)
         DOUBLE PRECISION  P     (1:NBAS,1:NBAS)
         DOUBLE PRECISION  XMAT  (1:NBAS,1:NBAS)

         DATA  ZERO   /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...set up the lower triangle of the ocupation matrix
C                in NBO basis.
C
C
         LTRG  = .TRUE.
         UTRG  = .TRUE.
         SAVEP = .FALSE.

         CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,
     +                NBAS,NBAS,
     +                0,
     +                SAVEP,LTRG,UTRG,
     +                C,
     +                P,
     +                W,
     +
     +                        XMAT )
     +
     +
         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                XMAT,
     +
     +                        P )
     +
     +
C
C
C             ...refine the ocupation matrix. By this we mean to
C                set the offdiagonal elements of the occupation
C                matrix corresponding to zero and MAXOCC diagonal
C                elements equal to exactly zero, thereby eliminating
C                any chance of (possibly incorrect) symmetry mixing.
C
C
         CALL  NLO__REFINE_OCCUPATION_MATRIX
     +
     +              ( NBAS,NBAS,
     +                NBAS,
     +                MAXOCC,
     +                QSYMACC,
     +
     +                         P )
     +
     +
C
C
C             ...balance the ocupation matrix. By this we mean to
C                set the offdiagonal elements of the occupation
C                matrix equal if they differ only by the symmetry
C                accuracy.
C
C
         CALL  MAT__BALANCE_SYMMETRIC_MATRIX
     +
     +              ( NBAS,NBAS,
     +                NBAS,
     +                NBAS,
     +                PSYMACC,
     +                W,
     +
     +                         P )
     +
     +
         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                ( 6,
     +                  ' P matrix ',
     +                  NBAS,NBAS,
     +                  NBAS,0,
     +                  P )
     +
     +
C         CALL  MAT__C_EQ_ZERO_FLOAT
C     +
C     +              ( NBAS,NBAS,
C     +                NBAS,NBAS,
C     +                XMAT )
C     +
C     +
C         CALL  MAT__C_EQ_A_TIMES_B_FLOAT
C     +
C     +              ( NBAS,NBAS,
C     +                NBAS,NBAS,
C     +                NBAS,NBAS,
C     +                NBAS,NBAS,NBAS,
C     +                P,P,
C     +
C     +                         XMAT )
C     +
C     +
C         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
C     +
C     +                ( 6,
C     +                  ' P * P ',
C     +                  NBAS,NBAS,
C     +                  NBAS,0,
C     +                  XMAT )
C     +
C     +
C
C
C             ...invoke the symmetrized Jacobi rotations.
C
C
         REVERS = .FALSE.

         CALL  NLO__SYMMETRIC_JACOBI
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,
     +                NBAS,
     +                NOCC+1,NOCC,
     +                REVERS,
     +                USED,
     +                P,
     +
     +                         NITER,
     +                         W,
     +                         XMAT )
     +
     +
         WRITE (*,*) ' NLMO converged after ',NITER,' iterations! '
         WRITE (1,*) ' NLMO converged after ',NITER,' iterations! '
C
C
C             ...generate the NLMO coefficient matrix in AO basis.
C
C
         CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,NBAS,
     +                C,XMAT,
     +
     +                         P )
     +
     +
         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                P,
     +
     +                        C )
     +
     +
         CALL    MAT__PRINT_V_FLOAT_12_NOZEROS
     +
     +                ( 1,
     +                  ' NLMO weight vector ',
     +                  NBAS,
     +                  NBAS,
     +                  W )
     +
     +
         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                ( 1,
     +                  ' C (NLMO) in AO basis ',
     +                  NBAS,NBAS,
     +                  NBAS,NBAS,
     +                  C )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
