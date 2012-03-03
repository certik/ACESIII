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
         SUBROUTINE  NLO__PRINT_NAO_RESULTS
     +
     +                    ( UNITID,
     +                      NBAS,NATOM,
     +                      MXSHELL,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      ZATOM,
     +                      LSIZE,
     +                      POPCOR,POPVAL,POPRYD,POPSUM,
     +                      NAOTYP,
     +                      ANGSYM,
     +                      LOCAL,
     +                      W )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NAO_RESULTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints the final results of the natural
C                atomic population analysis to the output file specified
C                by its unit identification number.
C
C                The following is printed:
C
C                   1) The occupation numbers (weights) of the
C                      final NAO's and their characterization as
C                      Core, Valence and Rydberg type orbitals.
C
C                   2) A population chart indicating the natural
C                      atomic populations broken down in Core,
C                      Valence and Rydberg populations, together
C                      with the overall natural atomic charges
C                      for each atomic center.
C
C                   3) The NAO atomic localization map.
C
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    MXSHELL      =  largest l-shell value
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A (in increasing order!).
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    LSIZE (I)    =  I-th l-shell size
C                    POPxxx (A)   =  # of electrons populating the
C                                    Core, Valence, Rydberg and all
C                                    NAOs on atom A (xxx = COR,VAL,
C                                    RYD and SUM, respectively).
C                    NAOTYP (I)   =  integer 2,1 or 0, identifying
C                                    the type of NAO encountered for
C                                    the I-th NAO (= 2 for Core, = 1
C                                    for Valence and = 0 for Rydberg).
C                    ANGSYM (I)   =  angular symmetry of I-th NAO.
C                                    This is equal to the lowest l-shell
C                                    value that can mix with the I-th
C                                    NAO and is dependent on the site
C                                    symmetry of the atom to which the
C                                    I-th NAO belongs.
C                    LOCAL (A,I)  =  NAO atom locality map for each
C                                    atom A and I-th NAO.
C                    W            =  NAO weight vector.
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

         LOGICAL     FIRST

         CHARACTER*2   ATCHAR
         CHARACTER*1   BAR,DASH
         CHARACTER*3   LSTATE,LSYM
         CHARACTER*7   NAOCHAR

         CHARACTER*2   ATSYMB  (1:104)
         CHARACTER*3   LSYMB   (0:7)
         CHARACTER*7   NAOSYMB (0:2)

         INTEGER     ATOM
         INTEGER     BASNR
         INTEGER     I,L,N
         INTEGER     LDIM,LTOT,LTYPE
         INTEGER     MXSHELL
         INTEGER     NBAS,NATOM
         INTEGER     NAL
         INTEGER     NLTYPE
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     ANGSYM  (1:NBAS)
         INTEGER     LSIZE   (0:MXSHELL)
         INTEGER     NAOTYP  (1:NBAS)
         INTEGER     NSHELLS (1:NATOM)
         INTEGER     ZATOM   (1:NATOM)

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)

         DOUBLE PRECISION  CHARGE,WEIGHT
         DOUBLE PRECISION  TOTCHG,TOTCOR,TOTVAL,TOTRYD,TOTSUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  POPCOR  (1:NATOM)
         DOUBLE PRECISION  POPVAL  (1:NATOM)
         DOUBLE PRECISION  POPRYD  (1:NATOM)
         DOUBLE PRECISION  POPSUM  (1:NATOM)
         DOUBLE PRECISION  W       (1:NBAS)

         DOUBLE PRECISION  LOCAL (1:NATOM,1:NBAS)

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
C             ...print out header.
C
C
         WRITE (UNITID,8000) 'Natural atomic orbital analysis'
 8000    FORMAT (//,20X,A31,//)
C
C
C             ...print out header of occupation analysis.
C
C
         WRITE (UNITID,9000) 'NAO #',BAR,'Atom',BAR,'Atom #',BAR,
     +                       'L-state',BAR,'L-sym',BAR,'NAO-type',BAR,
     +                       'Occupancy'
         WRITE (UNITID,9010) (DASH,I=1,72)
C
C
C             ...print out occupation analysis.
C
C
         BASNR = 0
         DO 1000 ATOM = 1,NATOM
            ZVAL = ZATOM (ATOM)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            FIRST = .TRUE.
            NLTYPE = NSHELLS (ATOM)
            DO 2000 N = 1,NLTYPE

               NAL = NBASAL (N,ATOM)
               LTYPE = SHELLS (N,ATOM)
               LDIM = LSIZE (LTYPE)
               LTOT = NAL * LDIM

               IF (LTYPE.GT.6) THEN
                   LSTATE = LSYMB (7)
               ELSE
                   LSTATE = LSYMB (LTYPE)
               END IF

               LTYPE = ANGSYM (BASNR+1)

               IF (LTYPE.GT.6) THEN
                   LSYM = LSYMB (7)
               ELSE
                   LSYM = LSYMB (LTYPE)
               END IF

               NAOCHAR = NAOSYMB (NAOTYP (BASNR+1))
               WEIGHT = W (BASNR+1)

               IF (FIRST) THEN
                   WRITE (UNITID,9020) BASNR+1,BAR,ATCHAR,BAR,ATOM,BAR,
     +                                 LSTATE,BAR,LSYM,BAR,NAOCHAR,BAR,
     +                                 WEIGHT
                   FIRST = .FALSE.
               ELSE
                   WRITE (UNITID,9030) BASNR+1,BAR,BAR,BAR,
     +                                 LSTATE,BAR,LSYM,BAR,NAOCHAR,BAR,
     +                                 WEIGHT
               END IF

               DO 100 L = 2,LTOT
                  LTYPE = ANGSYM (BASNR+L)
                  IF (LTYPE.GT.6) THEN
                      LSYM = LSYMB (7)
                  ELSE
                      LSYM = LSYMB (LTYPE)
                  END IF

                  NAOCHAR = NAOSYMB (NAOTYP (BASNR+L))
                  WEIGHT = W (BASNR+L)

                  WRITE (UNITID,9030) BASNR+L,BAR,BAR,BAR,
     +                                LSTATE,BAR,LSYM,BAR,NAOCHAR,BAR,
     +                                WEIGHT
  100          CONTINUE

               BASNR = BASNR + LTOT

 2000       CONTINUE
            WRITE (UNITID,9010) (DASH,I=1,72)
 1000    CONTINUE
C
C
C             ...formats for printing occupation analysis.
C
C
 9000    FORMAT (1X,A5,1X, A1, 1X,A4,1X, A1, 1X,A6,1X, A1,
     +           1X,A7,1X, A1, 1X,A5,1X, A1, 1X,A8,1X, A1,
     +           5X,A9,5X)
 9010    FORMAT (1X,90A1)
 9020    FORMAT (2X,I3,2X, A1, 2X,A2,2X, A1, 2X,I3,3X, A1,
     +           3X,A3,3X, A1, 2X,A3,2X, A1, 2X,A7,1X, A1,
     +           1X,F16.14)
 9030    FORMAT (2X,I3,2X, A1, 2X,2X,2X, A1, 2X,3X,3X, A1,
     +           3X,A3,3X, A1, 2X,A3,2X, A1, 2X,A7,1X, A1,
     +           1X,F16.14)
C
C
C             ...print out header for the population analysis.
C
C
         WRITE (UNITID,9040) 'Population Analysis'
 9040    FORMAT (//,20X,A19,//)

         WRITE (UNITID,9050) 'Atom',BAR,'Atom #',BAR,'Charge',BAR,
     +                       'Core',BAR,'Valence',BAR,
     +                       'Rydberg',BAR,'Sum (C+V+R)'
         WRITE (UNITID,9060) (DASH,I=1,90)
C
C
C             ...print out the population analysis.
C
C
         TOTCHG = ZERO
         TOTCOR = ZERO
         TOTVAL = ZERO
         TOTRYD = ZERO
         TOTSUM = ZERO

         DO 3000 ATOM = 1,NATOM
            ZVAL = ZATOM (ATOM)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            CHARGE = DFLOAT (ZVAL) - POPSUM (ATOM)

            TOTCHG = TOTCHG + CHARGE
            TOTCOR = TOTCOR + POPCOR (ATOM)
            TOTVAL = TOTVAL + POPVAL (ATOM)
            TOTRYD = TOTRYD + POPRYD (ATOM)
            TOTSUM = TOTSUM + POPSUM (ATOM)

            WRITE (UNITID,9070) ATCHAR,BAR,ATOM,BAR,CHARGE,BAR,
     +                          POPCOR (ATOM),BAR,POPVAL (ATOM),BAR,
     +                          POPRYD (ATOM),BAR,POPSUM (ATOM)
 3000    CONTINUE
         WRITE (UNITID,9060) (DASH,I=1,90)
         WRITE (UNITID,9080) 'Total',BAR,TOTCHG,BAR,TOTCOR,BAR,
     +                       TOTVAL,BAR,TOTRYD,BAR,TOTSUM

 9050    FORMAT (1X,A4,1X,     A1, 1X,A6,1X,     A1, 4X,A6,4X,     A1,
     +           5X,A4,5X,     A1, 4X,A7,3X,     A1,
     +           4X,A7,3X,     A1, 2X,A11)
 9060    FORMAT (1X,90A1)
 9070    FORMAT (2X,A2,2X,     A1, 2X,I3,3X,     A1, 1X,F12.8,1X,  A1,
     +           1X,F12.8,1X,  A1, 1X,F12.8,1X,  A1,
     +           1X,F12.8,1X,  A1, 1X,F12.8)
 9080    FORMAT (           5X,A5,5X,            A1, 1X,F12.8,1X,  A1,
     +           1X,F12.8,1X,  A1, 1X,F12.8,1X,  A1,
     +           1X,F12.8,1X,  A1, 1X,F12.8)
C
C
C             ...print out the NAO atomic localization map.
C
C
         WRITE (UNITID,9090) 'NAO Atomic Localization Map'
 9090    FORMAT (//,20X,A27)

         CALL    MAT__PRINT_A_FLOAT_3_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows = Atoms ; Columns = NAOs ',
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
