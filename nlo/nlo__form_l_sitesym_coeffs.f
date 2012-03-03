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
         SUBROUTINE  NLO__FORM_L_SITESYM_COEFFS
     +
     +                    ( NBAS,MXLSIZE,
     +                      LTOT,LDIM,NAL,
     +                      TSYM,
     +                      P,
     +                      XVEC,
     +                      XMAT,
     +
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_L_SITESYM_COEFFS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms occupation matrix L-blocks from
C                an input set of coefficients and diagonalizes them,
C                thus reflecting the site symmetry of the atom. The
C                total number of coefficient vectors in matrix C is
C                LTOT, and LDIM of these correspond to a complete m-set
C                for a particular L block. Hence the size of each
C                occupation matrix L-block is LDIM and there are a
C                total of NAL = LTOT/LDIM of these blocks.
C
C                As an example, suppose L = 1 (i.e. p-functions) and
C                thus LDIM = 3 corresponding to the three px,py,pz.
C                If LTOT = 9, then there are 3 p-blocks to consider
C                for diagonalization, each representing one px,py,pz
C                set.
C
C                The order of the coefficient matrix at this stage
C                is L-block wise.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    MXLSIZE      =  largest l-size value, i.e. maximum
C                                    m-degeneracy.
C                    LTOT         =  size of l-space
C                    LDIM         =  size of m-space
C                    NAL          =  # of l-shells
C                    TSYM         =  MXLSIZE x LDIM matrix that will
C                                    be used to hold the l-sitesymmetry
C                                    transformation matrix.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    C            =  current NAO coefficient matrix
C                                    of size NBAS x LTOT before site
C                                    symmetrization.
C
C                  Output:
C
C                    C            =  l-sitesymmetrized NAO coefficient
C                                    matrix
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

         INTEGER     BASNR
         INTEGER     L
         INTEGER     LTOT,LDIM
         INTEGER     MXLSIZE
         INTEGER     NAL
         INTEGER     NBAS

         DOUBLE PRECISION  BALANCE

         DOUBLE PRECISION  XVEC (1:LDIM)

         DOUBLE PRECISION  C    (1:NBAS   ,1:LTOT)
         DOUBLE PRECISION  P    (1:NBAS   ,1:NBAS)
         DOUBLE PRECISION  TSYM (1:MXLSIZE,1:LDIM)
         DOUBLE PRECISION  XMAT (1:NBAS   ,1:LDIM)

         DATA  REVERS   /.TRUE./

         PARAMETER  (BALANCE = 1.D-12)
C
C
C------------------------------------------------------------------------
C
C
C             ...loop over all L-blocks.
C
C
         LTRG  = .TRUE.
         UTRG  = .FALSE.
         SAVEP = .TRUE.

         BASNR = 0

         DO L = 1,NAL

            CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                 ( NBAS,LDIM,
     +                   NBAS,NBAS,
     +                   NBAS,LDIM,
     +                   LDIM,
     +                   LDIM,NBAS,
     +                   0,
     +                   SAVEP,LTRG,UTRG,
     +                   C (1,BASNR+1),
     +                   P,
     +                   XVEC,
     +
     +                           XMAT )
     +
     +
            CALL  MAT__C_EQ_A_FLOAT
     +
     +                 ( NBAS,LDIM,
     +                   MXLSIZE,LDIM,
     +                   LDIM,LDIM,
     +                   XMAT,
     +
     +                           TSYM )
     +
     +
C            CALL  MAT__PRINT_A_FLOAT_18_NOZEROS
C     +
C     +                 ( 6,
C     +                   ' TSYM matrix before balance (ltrg) ',
C     +                   MXLSIZE,LDIM,
C     +                   LDIM,0,
C     +                   TSYM )
C     +
C     +
            CALL  MAT__BALANCE_SYMMETRIC_MATRIX
     +
     +                 ( MXLSIZE,LDIM,
     +                   LDIM,
     +                   LDIM,
     +                   BALANCE,
     +                   XVEC,
     +
     +                           TSYM )
     +
     +
C            CALL  MAT__PRINT_A_FLOAT_18_NOZEROS
C     +
C     +                 ( 6,
C     +                   ' TSYM matrix after balance (ltrg) ',
C     +                   MXLSIZE,LDIM,
C     +                   LDIM,0,
C     +                   TSYM )
C     +
C     +
            CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                 ( MXLSIZE,LDIM,LDIM,
     +                   LDIM,
     +                   REVERS,
     +
     +                           XVEC,
     +                           TSYM )
     +
     +
C            CALL  MAT__PRINT_A_FLOAT_18_NOZEROS
C     +
C     +                ( 6,
C     +                  ' TSYM eigenvalues ',
C     +                  LDIM,1,
C     +                  LDIM,1,
C     +                  XVEC )
C     +
C     +
C            CALL  MAT__PRINT_A_FLOAT_18_NOZEROS
C     +
C     +                 ( 6,
C     +                   ' TSYM eigenvectors ',
C     +                   MXLSIZE,LDIM,
C     +                   LDIM,LDIM,
C     +                   TSYM )
C     +
C     +
            CALL  MAT__C_EQ_A_FLOAT
     +
     +                 ( NBAS,LDIM,
     +                   NBAS,LDIM,
     +                   NBAS,LDIM,
     +                   C (1,BASNR+1),
     +
     +                           XMAT )
     +
     +
            CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                 ( NBAS,LDIM,
     +                   MXLSIZE,LDIM,
     +                   NBAS,LDIM,
     +                   NBAS,LDIM,LDIM,
     +                   XMAT,TSYM,
     +
     +                           C (1,BASNR+1) )
     +
     +
            BASNR = BASNR + LDIM

         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
