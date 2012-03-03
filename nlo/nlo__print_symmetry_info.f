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
         SUBROUTINE  NLO__PRINT_SYMMETRY_INFO
     +
     +                    ( UNITID,
     +                      NATOM,
     +                      ZATOM,
     +                      SYMMAP,
     +                      RING,NRING,RINGSZ,
     +                      PLATO )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_SYMMETRY_INFO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints the symmetry info, which consists
C                of the symmetry map, the symmetry rings and its sizes
C                and the atoms related by platonic symmetry, to the
C                output file specified by its unit identification
C                number.
C
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NATOM        =  total # of atoms
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    SYMMAP (A,B) =  has been set equal to 1 if atoms
C                                    A and B were found to be symmetry
C                                    related. A value of 0 indicates
C                                    no symmetry relation.
C                    RING (A,B)   =  integer array containing all atomic
C                                    indices A of those symmetry related
C                                    sets that have equal distance from
C                                    center B and are of maximum size
C                    NRING (B)    =  integer vector containing the # of
C                                    sets of maximum size related to
C                                    center B.
C                    RINGSZ (B)   =  integer vector containing the
C                                    maximum size of the above found
C                                    sets related to center B.
C                    PLATO (A,P)  =  integer array containing all atomic
C                                    indices A for specific platonic
C                                    symmetry types as characterized by
C                                    the index P: P = 1 (four atoms
C                                    tetrahedrally arranged), P = 2
C                                    (eight atoms arranged in a cube),
C                                    P = 3 (six octahedral atoms), P = 4
C                                    (twenty atoms in a dodecahedron),
C                                    P = 5 (twelve icosahedral atoms)
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

         CHARACTER*2   ATCHAR,ATTYPE
         CHARACTER*1   BAR,DASH
         CHARACTER*12  SYMTYP

         CHARACTER*2   ATSYMB  (1:104)

         INTEGER     ATOM
         INTEGER     ATRING
         INTEGER     I,M,N
         INTEGER     NATOM
         INTEGER     NSET
         INTEGER     NSTEP
         INTEGER     OFFR
         INTEGER     REST
         INTEGER     RINGS
         INTEGER     SIZE
         INTEGER     STEP
         INTEGER     SYM
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     NRING   (1:NATOM)
         INTEGER     RINGSZ  (1:NATOM)
         INTEGER     ZATOM   (1:NATOM)

         INTEGER     PLATO   (1:NATOM,1:5    )
         INTEGER     RING    (1:NATOM,1:NATOM)
         INTEGER     SYMMAP  (1:NATOM,1:NATOM)

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
C
C
C------------------------------------------------------------------------
C
C
C             ...print out the atomic symmetry map.
C
C
         WRITE (UNITID,9000) 'Atomic Symmetry Map'
 9000    FORMAT (///,20X,A19)

         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows,Columns = Atoms ; 1 = Symmetry related ',
     +                  NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  SYMMAP )
     +
     +
C
C
C             ...print out the atomic symmetry rings info.
C
C
         WRITE (UNITID,9010) 'Ring Symmetry Info'
 9010    FORMAT (///,20X,A18,///)
C
C
C             ...print out header of atomic symmetry rings.
C
C
         WRITE (UNITID,9020) 'Atom',BAR,'Atom #',BAR,
     +                       'Ring size',BAR,'Ring #',BAR,
     +                       'Ring atom type',BAR,'Ring atomic indices'
         WRITE (UNITID,9030) (DASH,I=1,94)
C
C
C             ...print out atomic symmetry rings.
C
C
         NSTEP = 10

         DO ATOM = 1,NATOM
            ZVAL = ZATOM (ATOM)
            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            RINGS = NRING (ATOM)
            IF (RINGS.GT.0) THEN
                SIZE = RINGSZ (ATOM)
                STEP = (SIZE / NSTEP) - 1
                REST = MOD (SIZE,NSTEP)
                OFFR = 0

                DO N = 1,RINGS
                   ATRING = RING (OFFR+1,ATOM)
                   ZVAL = ZATOM (ATRING)
                   IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                       ATTYPE = ATSYMB (104)
                   ELSE
                       ATTYPE = ATSYMB (ZVAL)
                   END IF

                   M = MIN (SIZE,NSTEP)
                   IF (N.EQ.1) THEN
                       WRITE (UNITID,9040) ATCHAR,BAR,ATOM,BAR,
     +                                     SIZE,BAR,N,BAR,ATTYPE,BAR,
     +                                     (RING (OFFR+I,ATOM),I=1,M)
                   ELSE
                       WRITE (UNITID,9050) BAR,BAR,BAR,N,BAR,ATTYPE,BAR,
     +                                     (RING (OFFR+I,ATOM),I=1,M)
                   END IF
                   OFFR = OFFR + M
                   IF (SIZE.GT.NSTEP) THEN
                       DO M = 1,STEP
                          WRITE (UNITID,9060) BAR,BAR,BAR,BAR,BAR,
     +                                    (RING (OFFR+I,ATOM),I=1,NSTEP)
                          OFFR = OFFR + NSTEP
                       END DO
                       WRITE (UNITID,9060) BAR,BAR,BAR,BAR,BAR,
     +                                    (RING (OFFR+I,ATOM),I=1,REST)
                       OFFR = OFFR + REST
                   END IF
                END DO
            END IF
         END DO
C
C
C             ...formats for printing atomic symmetry rings.
C
C
 9020    FORMAT (1X,A4,1X,     A1, 1X,A6,1X,     A1,
     +           1X,A9,1X,     A1, 1X,A6,1X,     A1,
     +           1X,A14,1X,    A1, 10X,A19)
 9030    FORMAT (1X,94A1)
 9040    FORMAT (2X,A2,2X,     A1, 2X,I3,3X,     A1,
     +           3X,I3,5X,     A1, 2X,I3,3X,     A1,
     +           7X,A2,7X,     A1, 10(I4))
 9050    FORMAT (2X,2X,2X,     A1, 2X,3X,3X,     A1,
     +           3X,3X,5X,     A1, 2X,I3,3X,     A1,
     +           7X,A2,7X,     A1, 10(I4))
 9060    FORMAT (2X,2X,2X,     A1, 2X,3X,3X,     A1,
     +           3X,3X,5X,     A1, 2X,3X,3X,     A1,
     +           7X,2X,7X,     A1, 10(I4))
C
C
C             ...print out the platonic symmetry info.
C
C
         WRITE (UNITID,9070) 'Platonic Symmetry Info'
 9070    FORMAT (///,20X,A22,///)
C
C
C             ...print out header of platonic symmetry.
C
C
         WRITE (UNITID,9080) 'Symmetry',BAR,'Set #',BAR,
     +                       'Atom type',BAR,'Atomic indices'
         WRITE (UNITID,9090) (DASH,I=1,50)
C
C
C             ...print out platonic symmetry.
C
C
         DO 1000 SYM = 1,5
            IF (SYM.EQ.1) THEN
                SIZE = 4
                SYMTYP = ' Tetrahedral'
            ELSE IF (SYM.EQ.2) THEN
                SIZE = 8
                SYMTYP = '   Cubical  '
            ELSE IF (SYM.EQ.3) THEN
                SIZE = 6
                SYMTYP = ' Octahedral '
            ELSE IF (SYM.EQ.4) THEN
                SIZE = 12
                SYMTYP = 'Dodecahedral'
            ELSE
                SIZE = 20
                SYMTYP = ' Icosahedral'
            END IF

            NSET = 0

            DO N = 1,NATOM,SIZE
               ATOM = PLATO (N,SYM)
               IF (ATOM.EQ.0) THEN
                   IF (NSET.EQ.0) THEN
                       WRITE (UNITID,9100) SYMTYP,BAR,DASH,BAR,
     +                                     DASH,BAR,DASH
                   END IF
                   GOTO 1000
               ELSE
                   ZVAL = ZATOM (ATOM)
                   IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                       ATCHAR = ATSYMB (104)
                   ELSE
                       ATCHAR = ATSYMB (ZVAL)
                   END IF

                   NSET = NSET + 1
                   IF (NSET.EQ.1) THEN
                       WRITE (UNITID,9110) SYMTYP,BAR,NSET,BAR,
     +                                     ATCHAR,BAR,
     +                                     (PLATO (N+I-1,SYM),I=1,SIZE)
                   ELSE
                       WRITE (UNITID,9120) BAR,NSET,BAR,ATCHAR,BAR,
     +                                     (PLATO (N+I-1,SYM),I=1,SIZE)
                   END IF
               END IF
            END DO

 1000    CONTINUE
C
C
C             ...formats for printing platonic symmetry.
C
C
 9080    FORMAT (3X,A8,3X,     A1, 1X,A5,1X,     A1,
     +           1X,A9,1X,     A1, 1X,A14)
 9090    FORMAT (1X,50A1)
 9100    FORMAT (1X,A12,1X,    A1, 3X,A1,3X,     A1,
     +           5X,A1,5X,     A1, 7X,A1)
 9110    FORMAT (1X,A12,1X,    A1, 2X,I2,3X,     A1,
     +           4X,A2,5X,     A1, 20(I4))
 9120    FORMAT (1X,12X,1X,    A1, 2X,I2,3X,     A1,
     +           4X,A2,5X,     A1, 20(I4))
C
C
C             ...ready!
C
C
         RETURN
         END
