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
      subroutine iputrec(iFlag,szArchive,szRecName,iRecLen,iSrc)
c-----------------------------------------------------------------------------
c   Writes integer data to JOBARC file.
c
c   Note: The data on JOBARC is stored as integer*8 for compatibility with
c         the serial ACES II program's format.  So this subroutine converts
c         the data in iDest to 8-byte integers before calling the actual
c         PUTREC routine to write it to JOBARC.
c
c   iflag == 0 --> A record of iRecLen ints is set to 0 and written to JOBARC.
c         != 0 --> Data in iSrc is converted to 8-byte ints, then written 
c                    to JOBARC.
c-----------------------------------------------------------------------------

      implicit none

      integer iFlag, iRecLen, iSrc(*)
      character*(*) szArchive, szRecName

      integer*8 iSrc8(iRecLen)
      integer i

      if (iFlag .ne. 0) then
         
c-----------------------------------------------------------------------------
c   Convert the data.
c-----------------------------------------------------------------------------

         do i = 1, iRecLen
            iSrc8(i) = iSrc(i)
         enddo 
      endif

c----------------------------------------------------------------------------
c   Write the data.
c----------------------------------------------------------------------------

      call putrec(iFlag,szArchive,szRecName,iRecLen,isrc8)
      return
      end
