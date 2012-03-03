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
         SUBROUTINE  NLO__FORM_ATOMIC_NONVALENCE_NBO
     +
     +                    ( NBAS,MXNBA,
     +                      NXB,NXA,
     +                      ATNXB,ATXOFF,
     +                      P,
     +                      H,
     +                      XMAT,
     +
     +                              W,
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__FORM_ATOMIC_NONVALENCE_NBO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine establishes the atomic nonvalence (Core
C                and Rydberg) NBO coefficient matrix in terms of the
C                AO basis and the corresponding weights. The NHO info
C                sitting in array H is simply transferred to the
C                coefficient matrix C.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    MXNBA        =  maximum # of NHOs per atom.
C                    NXB          =  total # of nonvalence NBOs.
C                    NXA          =  total # of nonvalence NBO atoms.
C                    ATNXB (A)    =  # of atomic nonvalence NBOs on
C                                    atom A having nonvalence NBOs.
C                    ATXOFF (A)   =  index offset for atomic nonvalence
C                                    NBOs for atom A having nonvalence
C                                    NBOs. This index is equal to the
C                                    total number of atomic nonvalence
C                                    NBOs on all nonvalence NBO atoms
C                                    preceeding nonvalence NBO atom A.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    H (I,J)      =  MXNBA x NXB matrix containing the
C                                    atomic nonvalence NHOs. I is the
C                                    local atomic index labeling the
C                                    atomic NAOs from which the
C                                    nonvalence NHOs were constructed.
C                                    J is the global NHO index running
C                                    over all nonvalence NXB space
C                                    NHOs, with all NHOs belonging to
C                                    a specific atomic center being
C                                    grouped together.
C                    XMAT         =  flp scratch array of matrix type.
C                    W            =  section of the NAO weights vector
C                                    corresponding to the nonvalence
C                                    NAOs treated here.
C                    C            =  section of the NAO coefficient
C                                    matrix in AO basis corresponding
C                                    to the nonvalence NAOs treated
C                                    here.
C
C
C                  Output:
C
C                    W            =  section of the NBO weights vector
C                                    corresponding to the nonvalence
C                                    NBOs treated here.
C                    C            =  section of the NBO coefficient
C                                    matrix in AO basis corresponding
C                                    to the nonvalence NBOs treated
C                                    here.
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
         LOGICAL     SAVEP

         INTEGER     ATOM
         INTEGER     I,J,N
         INTEGER     MXNBA
         INTEGER     NAO,NBO
         INTEGER     NBAS
         INTEGER     NXB,NXA
         INTEGER     NXBA
         INTEGER     OFFA

         INTEGER     ATNXB   (1:NXA)
         INTEGER     ATXOFF  (1:NXA)

         DOUBLE PRECISION  X

         DOUBLE PRECISION  W  (1:NXB)

         DOUBLE PRECISION  C     (1:NBAS ,1:NXB  )
         DOUBLE PRECISION  H     (1:MXNBA,1:NXB  )
         DOUBLE PRECISION  P     (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  XMAT  (1:NBAS ,1:MXNBA)
C
C
C------------------------------------------------------------------------
C
C
C             ...establish the coefficient matrix C and weights
C                using the NHO info transmitted in array H.
C
C
         LTRG = .FALSE.
         UTRG = .FALSE.
         SAVEP = .TRUE.

         DO ATOM = 1,NXA

            OFFA = ATXOFF (ATOM)
            NXBA = ATNXB  (ATOM)

            CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +                 ( NBAS,NXBA,
     +                   NBAS,NXBA,
     +
     +                          XMAT )
     +
     +
            DO J = 1,NXBA
               NBO = OFFA + J
               DO I = 1,NXBA
                  NAO = OFFA + I
                  X = H (I,NBO)
                  DO N = 1,NBAS
                     XMAT (N,J) = XMAT (N,J) + X * C (N,NAO)
                  END DO
               END DO
            END DO

            CALL  MAT__C_EQ_A_FLOAT
     +
     +                 ( NBAS,NXBA,
     +                   NBAS,NXBA,
     +                   NBAS,NXBA,
     +                   XMAT,
     +
     +                           C (1,OFFA+1) )
     +
     +
            CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                 ( NBAS,NXBA,
     +                   NBAS,NBAS,
     +                   NBAS,NXBA,
     +                   NXBA,
     +                   NXBA,NBAS,
     +                   0,
     +                   SAVEP,LTRG,UTRG,
     +                   C (1,OFFA+1),
     +                   P,
     +
     +                           W (OFFA+1),
     +
     +                   XMAT )
     +
     +
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
