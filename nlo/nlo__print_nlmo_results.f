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
         SUBROUTINE  NLO__PRINT_NLMO_RESULTS
     +
     +                    ( UNITID,
     +                      NBAS,NATOM,
     +                      LOCAL )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PRINT_NLMO_RESULTS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine prints info about the set of obtained
C                natural localized molecular orbitals (NLMOs) to a file
C                specified by its unit identification number. The
C                order of the NLMOs corresponds to the parent NBO
C                ordering.
C
C                At the moment only the atomic localization map of
C                the NLMOs is printed, which, when compared to the
C                one obtained for the NBOs, gives an indication on
C                how much locality was lost due to the NBO -> NLMO
C                transformation.
C
C
C                  Input:
C
C                    UNITID       =  printout unit identification #
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    LOCAL (A,I)  =  NLMO atom locality map for each
C                                    atom A and with NLMOs in NCB/NLB/
C                                    NBB/NEB/NAB/NYB/NRB order (i.e.
C                                    corresponding to their parent NBO
C                                    ordering).
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

         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     UNITID

         DOUBLE PRECISION  LOCAL (1:NATOM,1:NBAS)
C
C
C------------------------------------------------------------------------
C
C
C             ...print out the NLMO atomic localization map.
C
C
         WRITE (UNITID,9000) 'NLMO Atomic Localization Map'
 9000    FORMAT (//,20X,A28)

         CALL    MAT__PRINT_A_FLOAT_3_NOZEROS
     +
     +                ( UNITID,
     +                  ' Rows = Atoms ; Columns = NLMOs ',
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
