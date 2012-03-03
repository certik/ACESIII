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
         SUBROUTINE  NLO__PRINT_MOLECULE_INFO
     +
     +                    ( UNITID,
     +                      NATOM,
     +                      XYZ,ZATOM,
     +                      GPOINT,
     +                      DOMAT,
     +                      SOMAT,
     +                      BOMAT )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_MOLECULE_INFO
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints the molecular info such as atom
C                coordinates, several atomic order matrices, etc to
C                the output file specified by its unit identification
C                number.
C
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NATOM        =  total # of atoms
C                    XYZ          =  x,y,z-coordinates for all the
C                                    atoms in 3 x NATOM matrix format.
C                    ZATOM (I)    =  atomic number for I-th atom.
C                    GPOINT       =  if the center of mass of the
C                                    atomic arrangement coincides with
C                                    one of the atoms, this variable
C                                    will contain the corresponding
C                                    atomic index. If no atom coincided
C                                    with the center of mass location
C                                    it contains a value of 0.
C                    DOMAT        =  distance order matrix
C                    SOMAT        =  atomic overlap order matrix
C                    BOMAT        =  bond order matrix
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

         CHARACTER*3   ANSWER
         CHARACTER*2   ATCHAR
         CHARACTER*1   BAR,DASH

         CHARACTER*2   ATSYMB  (1:104)

         INTEGER     ATOM
         INTEGER     GPOINT
         INTEGER     I
         INTEGER     NATOM
         INTEGER     UNITID
         INTEGER     ZVAL

         INTEGER     ZATOM   (1:NATOM)

         DOUBLE PRECISION  X,Y,Z

         DOUBLE PRECISION  BOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  DOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  SOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  XYZ   (1:3    ,1:NATOM)

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
C             ...print out header.
C
C
         WRITE (UNITID,8000) 'Molecular Info'
 8000    FORMAT (//,20X,A14,//)
C
C
C             ...print out header of atomic coordinates.
C
C
         WRITE (UNITID,9000) 'Atom',BAR,'Atom #',BAR,
     +                       'X coordinate',BAR,
     +                       'Y coordinate',BAR,
     +                       'Z coordinate',BAR,
     +                       'is at center of mass'
         WRITE (UNITID,9010) (DASH,I=1,99)
C
C
C             ...print out atomic coordinates.
C
C
         DO ATOM = 1,NATOM
            X = XYZ (1,ATOM)
            Y = XYZ (2,ATOM)
            Z = XYZ (3,ATOM)
            ZVAL = ZATOM (ATOM)

            IF (ZVAL.LT.1 .OR. ZVAL.GT.103) THEN
                ATCHAR = ATSYMB (104)
            ELSE
                ATCHAR = ATSYMB (ZVAL)
            END IF

            IF (ATOM.EQ.GPOINT) THEN
                ANSWER = 'YES'
            ELSE
                ANSWER = ' NO'
            END IF

            WRITE (UNITID,9020) ATCHAR,BAR,ATOM,BAR,
     +                          X,BAR,Y,BAR,Z,BAR,
     +                          ANSWER
         END DO
C
C
C             ...formats for printing atomic coordinates.
C
C
 9000    FORMAT (1X,A4,1X,     A1, 1X,A6,1X,     A1,
     +           3X,A12,5X,    A1, 3X,A12,5X,    A1, 3X,A12,5X,    A1,
     +           1X,A20)
 9010    FORMAT (1X,99A1)
 9020    FORMAT (2X,A2,2X,     A1, 2X,I3,3X,     A1,
     +           1X,F18.14,1X, A1, 1X,F18.14,1X, A1, 1X,F18.14,1X, A1,
     +           9X,A3,9X)
C
C
C             ...print out the distance order matrix.
C
C
         WRITE (UNITID,8010) 'Distance Order Matrix'
 8010    FORMAT (//,20X,A21)

         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows , Columns = Atom indices ',
     +                  NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  DOMAT )
     +
     +
C
C
C             ...print out the atomic overlap order matrix.
C
C
         WRITE (UNITID,8020) 'Atomic Overlap Order Matrix'
 8020    FORMAT (//,20X,A27)

         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows , Columns = Atom indices ',
     +                  NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  SOMAT )
     +
     +
C
C
C             ...print out the bond order matrix.
C
C
         WRITE (UNITID,8030) 'Bond Order Matrix'
 8030    FORMAT (//,20X,A17)

         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows , Columns = Atom indices ',
     +                  NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  BOMAT )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
