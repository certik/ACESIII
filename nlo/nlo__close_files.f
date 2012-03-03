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
         SUBROUTINE  NLO__CLOSE_FILES
     +
     +                    ( UNITID,UNITDB )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__OPEN_FILES
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This operation closes the two sequential files:
C
C                         UNITID : printout file
C                         UNITDB : diagnostic/debug file
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    UNITID,UNITDB
C
C
C------------------------------------------------------------------------
C
C
C             ...close the files.
C
C
         CALL  OFH__CLOSE_UNIT (UNITID)
         CALL  OFH__CLOSE_UNIT (UNITDB)
C
C
C             ...ready!
C
C
         RETURN
         END
