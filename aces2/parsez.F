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

c Single-line parser of Z-matrix elements ONLY.

c INPUT
c    zline(80) : the Z-matrix line to parse

c OUTPUT
c    izl(2,7) : the beginning (1,*) and ending (2,*) indices of each
c               element (non-existent elements are at zline(0:0))
c               EXAMPLE:
c                  zline = '  X  1  R3  2  A2  3  D1   '
c                           123456789+123456789+1234567
c                  izl(1:2,1:7) = 3  6   9  13  16  20  23
c                                 3  6  10  13  17  20  24

#include "iachar.h"

      subroutine parsez(zline,izl)
      implicit none

#include "linelen.par"

      character*(linelen) zline
      integer izl(2,7)

      character*1 czTmp
      integer ndx, count
      logical find_char, not_done

c "parameters"
      character*1 achar, czTab, czSpace, czComment
      integer    max_cols
      parameter (max_cols = 80)

c ----------------------------------------------------------------------

      czTab     = achar(_IACHAR_TAB)
      czSpace   = achar(_IACHAR_SPACE)
      czComment = achar(_IACHAR_POUND)

      izl(1,1) = 0
      izl(2,1) = 0
      izl(1,2) = 0
      izl(2,2) = 0
      izl(1,3) = 0
      izl(2,3) = 0
      izl(1,4) = 0
      izl(2,4) = 0
      izl(1,5) = 0
      izl(2,5) = 0
      izl(1,6) = 0
      izl(2,6) = 0
      izl(1,7) = 0
      izl(2,7) = 0

c   o the zline pointer
      ndx = 1

c   o the izl pointer
      count = 1

c   o start looking for a char
      find_char = .true.

      not_done = .true.
      do while (not_done)

         czTmp = zline(ndx:ndx)
         if (find_char) then
            if ( (czTmp.ne.czSpace) .and.
     &           (czTmp.ne.czTab  )       ) then
               if (czTmp.eq.czComment) then
                  not_done = .false.
               else
                  izl(1,count) = ndx
                  find_char = .false.
               end if
            end if
         else
            if ( (czTmp.eq.czSpace).or.
     &           (czTmp.eq.czTab  )    ) then
               izl(2,count) = ndx - 1
               if (count.eq.7) then
                  return
               else
                  find_char = .true.
                  count = count + 1
               end if
            else
               if (czTmp.eq.czComment) then
                  izl(2,count) = ndx - 1
                  not_done = .false.
               end if
            end if
         end if

         if (not_done.and.(ndx.eq.max_cols)) then
            not_done = .false.
            if (.not.find_char) izl(2,count) = ndx
         else
            ndx = ndx + 1
         end if

c     end do while (not_done)
      end do

c      write(*,*) (izl(1,count),count=1,7)
c      write(*,*) (izl(2,count),count=1,7)

      return
      end

