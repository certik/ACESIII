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
         SUBROUTINE  NLO__OPEN_FILES
     +
     +                    ( NAMEID,NAMEDB,
     +
     +                              UNITID,
     +                              UNITDB )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__OPEN_FILES
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This operation opens the two sequential files:
C
C                         NAMEID : printout file
C                         NAMEDB : diagnostic/debug file
C
C                and returns their unit numbers.
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         CHARACTER*(*)      NAMEID,NAMEDB

         INTEGER            UNITID,UNITDB
C
C
C------------------------------------------------------------------------
C
C
C             ...open the two files.
C
C
         CALL  OFH__OPEN_FILE
     +
     +              ( NAMEID,
     +                'UNKNOWN',
     +                'SEQUENTIAL',
     +                'FORMATTED',
     +
     +                        UNITID )
     +
     +
         CALL  OFH__OPEN_FILE
     +
     +              ( NAMEDB,
     +                'UNKNOWN',
     +                'SEQUENTIAL',
     +                'FORMATTED',
     +
     +                        UNITDB )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
