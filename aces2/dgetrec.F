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
      subroutine dgetrec(iFlag,szArchive,szRecName,iRecLen,dDest)
c-----------------------------------------------------------------------------
c   Reads double precision data from JOBARC file.
c
c   Note: The data on JOBARC is stored as integer*8 for compatibility with
c         the serial ACES II program's format.  So this subroutine converts
c         the data read from JOBARC into the internal integer format and
c         returns it in dDest.
c
c   iflag > 0 --> an existing record is retrieved, length stored in iRecLen.
c         = 0 --> record length is returned (or -1 if record doesn't exist).
c         < 0 --> existing record is retrieved or iDest(1:iRecLen) is cleared.
c-----------------------------------------------------------------------------

      implicit none

      integer iFlag, iRecLen
      character*(*) szArchive, szRecName
      double precision dDest(*)

#include "machsp.com"
      integer i, iRecLen_int
      integer iconvert

c-----------------------------------------------------------------------------
c   Convert the record length to the actual number of dp values to be returned.
c-----------------------------------------------------------------------------

      iconvert    = ifltln / 8
      iRecLen_int = iRecLen * iconvert

c-----------------------------------------------------------------------------
c   Call getrec to get the raw data.
c-----------------------------------------------------------------------------

      call getrec(iFlag,szArchive,szRecName,iRecLen_int,dDest)

c-----------------------------------------------------------------------------
c   Change the record length back to a number of ints.
c-----------------------------------------------------------------------------

      if (iFlag .eq. 0) iRecLen = iRecLen_int / iconvert
      return
      end
