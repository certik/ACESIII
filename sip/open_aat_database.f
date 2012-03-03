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
      subroutine open_aat_database(iunit, ierr_return)
c---------------------------------------------------------------------------
c   Opens the AAT_DATABASE file using the Fortran unit number specified in
c   "iunit".  If the file is opened successfully, this routine returns 
c   ierr_return = 0, otherwise an error has occurred in opening the file.
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'

      integer iunit, ierr_return
      integer i, n, lenenv, ierr
      character*256 fn
      character*120 envvar

c----------------------------------------------------------------------------
c   Strip any embedded nulls from the filename.
c----------------------------------------------------------------------------

      do i = 1, len(aat_database)
         if (aat_database(i:i) .eq. ' ' .or. 
     *       ichar(aat_database(i:i)) .eq. 0) then
            n = i-1
            go to 100
         else
            fn(i:i) = aat_database(i:i)
         endif
      enddo
 
c----------------------------------------------------------------------------
c   First, try to open the file using the copy in the current directory.
c----------------------------------------------------------------------------

  100 continue
      open (unit=iunit, file = fn, status = 'OLD',
     *      err=200, iostat = ierr)

c-----------------------------------------------------------------------------
c   The open was successful.
c-----------------------------------------------------------------------------

      ierr_return = 0
      return   
      
c----------------------------------------------------------------------------
c   Next, try to open the copy in the ACES_EXE_PATH directory.
c----------------------------------------------------------------------------
       
  200 continue
      call c_getenv('ACES_EXE_PATH'//char(0), envvar, lenenv, ierr)
      if (ierr .eq. 0) then
         fn = envvar(1:lenenv) // '/sio/' // aat_database(1:n)
         open (unit=iunit, file = fn, status = 'OLD',
     *      err=300, iostat = ierr)

c---------------------------------------------------------------------------
c   Successful open.
c---------------------------------------------------------------------------

         ierr_return = 0
         return

c---------------------------------------------------------------------------
c   The open failed.  Return an error.
c---------------------------------------------------------------------------

  300    continue
         ierr_return = ierr
      else

c---------------------------------------------------------------------------
c   The first attempt failed, and there is no env. variable coded.  
c   Return an error code. 
c---------------------------------------------------------------------------

         ierr_return = ierr
      endif

      return
      end
