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
      subroutine init_params(pfile)
c-------------------------------------------------------------------------
c   Reads parameters from a ZMAT file into an internal table.
c   The ZMAT file must be present in the run directory.
c-------------------------------------------------------------------------

      implicit none
      character*(*) pfile

      integer max_table 
      integer ntable
      integer iunit
      parameter (max_table = 1000)

      character*80 table(max_table)
      common /parameter_table/ntable, table

      integer i, ios, n
      integer str_trimlen
      character*80 data

      ntable = 0

c---------------------------------------------------------------------------
c   Open the parameter file.
c---------------------------------------------------------------------------

      iunit = 12
      open (file = pfile, unit = iunit, status = 'OLD',
     *      err=500, iostat = ios)

c---------------------------------------------------------------------------
c   Read each record down to the "*SIP" identifier.
c---------------------------------------------------------------------------

   50 continue
      read (iunit,fmt='(a80)', end = 150, err = 600) data
      if (data(1:4) .ne. '*SIP') go to 50

      go to 200
  150 continue
c      print *,'Error reading ZMAT: No *SIP record was found.'
c      call abort_job()
      go to 300 
  200 continue

c---------------------------------------------------------------------------
c   Read all parameter data until a completely blank line (or a line 
c   beginning with an "*" is found.
c---------------------------------------------------------------------------

      do i = 1, max_table
         read (iunit,'(A80)', end=100) table(i)
         
         n = str_trimlen(table)     ! find length w/o trailing blanks.
         if (n .eq. 0) go to 100    ! n=0 implies a blank line
         if (table(i)(1:1) .eq. '*') go to 100

         ntable = ntable + 1
      enddo

      print *,'ERROR: Maximum number of lines of parameters (',
     *         max_table,') has been exceeded.'
      call abort_job()

  100 continue

      if (table(ntable+1)(1:7) .eq. '*DROPAO') then

c---------------------------------------------------------------------------
c   DROPAO namelist: Read in the lines to the table.
c-------------------------------------------------------------------------

         data = ' '
         do while (data(1:1) .ne. '*')
            read (iunit,fmt='(a80)', end = 300, err = 600) data
            n = str_trimlen(data) 
            if (n .eq. 0) go to 300   ! blank line
            table(ntable) = data
            ntable = ntable + 1
         enddo
      endif

  300 continue
c--------------------------------------------------------------------------
c   Close the parameter file.
c--------------------------------------------------------------------------

      close(iunit)
      return

  500 continue
      print *,'Error: Cannot open parameter file.  I/O status = ',ios
      print *,'       Parameter file name = ',pfile
      call abort_job()

  600 continue
      print *,'I/O error while reading parameter: I/O status ',ios
      print *,'       Parameter file name = ',pfile
      call abort_job()
      return

      end

      subroutine rgetparam(keyword, instance, rvalue)
c------------------------------------------------------------------------
c   Returns the value corresponding to "keyword".  If the keyword is not
c   found, it returns nothing.
c------------------------------------------------------------------------

      implicit none

      integer max_table
      parameter (max_table = 1000)

      integer ntable
      character*80 table(max_table)

      common /parameter_table/ntable, table

      character*(*) keyword
      integer instance, icount
      double precision rvalue
      integer ivalue
      integer i, n
      integer str_trimlen

      character*80 cvalue

      n = str_trimlen(keyword)
      icount = 0
      do i = 1, ntable
         if (keyword(1:n) .eq. table(i)(1:n)) then
            if (table(i)(n+1:n+1) .eq. ' ' .or.
     *          table(i)(n+1:n+1) .eq. '=') then
               icount = icount + 1
               if (icount .eq. instance) then
                  call c_line_parse(table(i) // char(0), '=', 2, ivalue,
     *                        rvalue, cvalue)
                  return
               endif
            endif
         endif 
      enddo
      return
      end

      subroutine igetparam(keyword, instance, ivalue)
c------------------------------------------------------------------------
c   Returns the value corresponding to "keyword".  If the keyword is not
c   found, it returns nothing.
c------------------------------------------------------------------------

      implicit none

      integer max_table
      parameter (max_table = 1000)
      integer ntable
      character*80 table(max_table)

      common /parameter_table/ntable, table

      character*(*) keyword
      integer instance, icount
      double precision rvalue
      integer ivalue
      integer i, n
      integer str_trimlen

      character*80 cvalue

      n = str_trimlen(keyword)
      icount = 0
      do i = 1, ntable
         if (keyword(1:n) .eq. table(i)(1:n)) then
            if (table(i)(n+1:n+1) .eq. ' ' .or.
     *          table(i)(n+1:n+1) .eq. '=') then
               icount = icount + 1
               if (icount .eq. instance) then
                  call c_line_parse(table(i) // char(0), '=', 1, ivalue,
     *                        rvalue, cvalue)
                  return
               endif
            endif
         endif 
      enddo
      return
      end

      subroutine cgetparam(keyword, instance, cvalue)
c------------------------------------------------------------------------
c   Returns the value corresponding to "keyword".  If the keyword is not
c   found, it returns nothing.
c------------------------------------------------------------------------

      implicit none

      integer max_table
      parameter (max_table = 1000)
      integer ntable
      character*80 table(max_table)

      common /parameter_table/ntable, table

      character*(*) keyword
      integer instance, icount
      double precision rvalue
      integer ivalue
      integer i, n, nn, j, indx
      integer str_trimlen

      character*(*) cvalue

      n = str_trimlen(keyword)
      icount = 0
      do i = 1, ntable
         if (keyword(1:n) .eq. table(i)(1:n)) then
            if (table(i)(n+1:n+1) .eq. ' ' .or.
     *          table(i)(n+1:n+1) .eq. '=') then
               icount = icount + 1
               if (icount .eq. instance) then
                  call c_line_parse(table(i) // char(0), '=', 3, ivalue,
     *                        rvalue, cvalue)
                  return
               endif
            endif
         endif 
      enddo
      return
      end

