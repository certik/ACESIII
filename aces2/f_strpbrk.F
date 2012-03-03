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

c This routine returns the first index in sz of any single char appearing in
c charset. WARNING: Since Fortran pads the end of all strings with spaces,
c be sure to pass only the relevant substring.

#define SZ_BUF      256
#define CHARSET_BUF 32

      integer function f_strpbrk(sz,charset)
      implicit none

      character*(*) sz, charset

      character*(SZ_BUF)      sz2
      character*(CHARSET_BUF) charset2

      character*1 achar
      integer f_strpbrk_core

      if (len(sz).ge.SZ_BUF) then
         print *, '@F_STRPBRK: The sz buffer is too small ',
     &            'to contain the input string.'
         print *, '                  Recompile with at least ',
     &            len(sz)+1, ' characters in the buffer.'
         print *, '                  (Currently ',SZ_BUF,
     &            ' characters.)'
         call c_exit(1)
      end if
      if (len(charset).ge.CHARSET_BUF) then
         print *, '@F_STRPBRK: The charset buffer is too small ',
     &            'to contain the input string.'
         print *, '                  Recompile with at least ',
     &            len(charset)+1, ' characters in the buffer.'
         print *, '                  (Currently ',CHARSET_BUF,
     &                                                   ' characters.)'
         call c_exit(1)
      end if

      sz2      = sz//achar(0)
      charset2 = charset//achar(0)

c      print *, 'sz is ',len(sz),' chars long'
c      print *, 'charset is ',len(charset),' chars long'

      f_strpbrk = f_strpbrk_core(sz2,charset2)

      return
      end

