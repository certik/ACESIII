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



c This routine removes a file (or directory) without having to open it and
c then close it with status='DELETE'.

c WARNING:
c Fortran blank-pads all strings. This routine will simply append the
c null character to szFile and pass it to f_remove_core. If there are trailing
c blanks in the file name, then f_remove_core will not remove it and the process
c will die.



      subroutine f_remove(szFile)
      implicit none

c ARGUMENTS
      character*(*) szFile

c EXTERNAL FUNCTIONS
      integer f_remove_core
      character*1 achar

c INTERNAL VARIABLES
      integer iLength, iTmp
      character*(256) sz

c ----------------------------------------------------------------------


      iLength = 1
      do while (szFile(iLength:iLength).ne.' '.and.
     &          iLength.le.len(szFile))
         iLength = iLength + 1
      end do
      iLength = iLength - 1
      if (iLength.eq.0) return

c ----------------------------------------------------------------------

      if (iLength.lt.256) then
         sz   = szFile(1:iLength)//achar(0)
         iTmp = f_remove_core(sz)
         if (iTmp.eq.0) return
         print *, '@F_REMOVE: The file "',szFile,
     &            '" could not be removed.'
         print *, '           error code = ',iTmp
      else
         print *, '@F_REMOVE: The sz buffer is too small ',
     &            'to contain the input string.'
         print *, '           Recompile with at least ',iLength+1,
     &            ' characters in the buffer.'
         print *, '           (Currently ',256,' characters.)'
      end if

      call c_exit(1)
      end

