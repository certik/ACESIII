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

c This routine removes all the records in the JOBARC file at the specified
c record and offset. It does not decrease the size of the JOBARC file, it merely
c resets the lookup table data.

c Example:
c  - erase all records after NEXTGEOM
c    call aces_ja_truncate('NEXTGEOM',1)

c INPUT
c char*(*) SZRECORD : the name of the record to measure the offset from
c int      IOFFSET  : the offset from the record index to truncate at

      subroutine aces_ja_truncate(szRecord,iOffset)
      implicit none

c ARGUMENTS
      character*(*) szRecord
      integer iOffset

c EXTERNAL FUNCTIONS
      integer iszeq

c INTERNAL VARIABLES
      integer ndx, i

c COMMON BLOCKS
#include "jobarc.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert job archive subsystem is up
      if (.not.bJAUp) then
         print *, '@ACES_JA_TRUNCATE: Assertion failed.'
         print *, '   szRecord = ',szRecord
         print *, '   bJAUp    = ',bJAUp
         i = 1
      end if
c   o the record name must be between 1 and 8 characters
      if ((len(szRecord).lt.1).or.(len(szRecord).gt.8)) then
         print *, '@ACES_JA_TRUNCATE: Assertion failed.'
         print *, '   szRecord = ',szRecord
         i = 1
      end if
c   o the record cannot be named 'OPENSLOT'
      if (len(szRecord).eq.8) then
      if (szRecord.eq.'OPENSLOT') then
         print *, '@ACES_JA_TRUNCATE: Assertion failed.'
         print *, '   szRecord = ',szRecord
         i = 1
      end if
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o get the index of the requested record
      ndx = iszeq(_MAX_JA_RECS,marker,1,szRecord)

c   o the record is not found or is at the last position
      if ((ndx.eq.0).or.(ndx.eq._MAX_JA_RECS)) return

c   o point to the truncation record
      ndx = ndx + iOffset

c   o the position is out of bounds
      if ((ndx.lt.1).or.(ndx.gt._MAX_JA_RECS)) return

c   o the record is already empty
      if (marker(ndx).eq.'OPENSLOT') return

c   o reset the remaining structures
      do i = ndx, _MAX_JA_RECS
         marker(i) = 'OPENSLOT'
      end do
      call izero(rloc(ndx), 1+_MAX_JA_RECS-ndx) 
      call izero(rsize(ndx),1+_MAX_JA_RECS-ndx) 

c   o mark JOBARC as modified
      bJAMod = .true.

      return
c     end subroutine aces_ja_truncate
      end

