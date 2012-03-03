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
      subroutine f_strtoken(source, delim, token)
c--------------------------------------------------------------------------
c   Returns the next token delimited by the character "delim".
c   If no delimiter is found, the NULL character is returned 
c   in "token".
c---------------------------------------------------------------------------
      implicit none
      character*(*) source, token
      character delim

      integer i, istart, istart2 

      character*256 src_save
      integer iptr, nsave 

      save src_save, iptr, nsave

      if (source(1:1) .eq. char(0)) then
c         print *,'Start from src_save : ',src_save(istart:nsave)
         istart = iptr
      else
         src_save = source
         istart = 1
         nsave  = len(source)
      endif

      token = ' '   ! blank-fill return argument.

c---------------------------------------------------------------------------
c   Skip past beginning whitespace.
c---------------------------------------------------------------------------

      do i = istart, nsave
         if (src_save(i:i) .ne. ' ' .and.
     *       src_save(i:i) .ne. char(0)) then
            istart2 = i
            go to 100
         endif
      enddo

c---------------------------------------------------------------------------
c   No non-blank, non-null characters exist.  Return a null.
c---------------------------------------------------------------------------

c      print *,'Return NULL token (#1)'
      token(1:1) = char(0)
      return
      
  100 continue
      do i = istart2, nsave
         if (src_save(i:i) .eq. delim) then

c---------------------------------------------------------------------------
c   We have a match.  Store the data in "token".
c---------------------------------------------------------------------------

            if (istart2 .eq. i) then
               token(1:1) = char(0)
c               print *,'Return NULL token (#2)'
            else
               token(1:i-istart2) = src_save(istart2:i-1) 
            endif

            iptr = i + 1   ! Next starting string.
            return
         endif
      enddo

c---------------------------------------------------------------------------
c   No more tokens.  Return a NULL.
c---------------------------------------------------------------------------

      token(1:1) = char(0) 
c      print *,'Return NULL token (#3)'
      return
      end


