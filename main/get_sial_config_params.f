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
      subroutine get_sial_config_params(sial_program)
c--------------------------------------------------------------------------
c   Scans for any program-specific parameters on the SIAL config file.
c--------------------------------------------------------------------------
      implicit none
      include 'sial_config_params.h'

      character*(*) sial_program
      character*80 line, keyword, cval
      character*256 config_file, envvar

      integer iunit, ival, lenenv, ierr
      integer i, n, nn, ios, str_trimlen
      double precision fval
 
c--------------------------------------------------------------------------
c   Default values.
c--------------------------------------------------------------------------

      no_servers = .false.
      ignore_dropmo = .false.
      use_2der_integrals = .false.
      vvvi_stack = .false.
 
c--------------------------------------------------------------------------
c   Open the SIAL config file.
c--------------------------------------------------------------------------

      config_file = './sial_config'   ! first try.

      iunit = 12
      open (unit=iunit, file = config_file, status = 'OLD',
     *      err=500, iostat = ios)

      print *,'Reading SIAL config file ',config_file

c--------------------------------------------------------------------------
c   Search file for a match to sial_program.
c--------------------------------------------------------------------------

   50 continue
      nn = str_trimlen(sial_program)
  100 continue
      read (iunit,'(A80)', end = 200) line
      if (line(1:nn) .ne. sial_program(1:nn)) go to 100

c---------------------------------------------------------------------------
c   Read parameters line.
c---------------------------------------------------------------------------

  120 continue
      read (iunit,'(A80)', end = 200) line

c---------------------------------------------------------------------------
c   Find the '=' sign.
c---------------------------------------------------------------------------

      cval = ' '
      call c_line_parse(line // char(0), '=',3,ival, fval, cval)
      if (cval(1:1) .eq. ' ') go to 200   ! no more params found.

c---------------------------------------------------------------------------
c   Find the keyword.
c---------------------------------------------------------------------------

      n = str_trimlen(line)
      do i = 1, n
         if (line(i:i) .eq. '=') then
            keyword = line(1:i-1)
            go to 150
         endif
      enddo
      
  150 continue
      n = str_trimlen(keyword)
  
      if (keyword(1:n) .eq. 'NO_SERVERS') then
         call c_line_parse(line//char(0),'=',3,ival, fval, cval)
         if (cval(1:3) .eq. 'YES') no_servers = .true.
      else if (keyword(1:n) .eq. 'IGNORE_DROPMO') then
         call c_line_parse(line//char(0),'=',3,ival, fval, cval)
         if (cval(1:3) .eq. 'YES') ignore_dropmo = .true.
      else if (keyword(1:n) .eq. 'SECOND_DERIVATIVE_INTEGRALS') then
         call c_line_parse(line//char(0),'=',3,ival, fval, cval)
         if (cval(1:3) .eq. 'YES') use_2der_integrals = .true.
      else if (keyword(1:n) .eq. 'VVVI_STACK') then
         call c_line_parse(line//char(0),'=',3,ival, fval, cval)
         if (cval(1:3) .eq. 'YES') vvvi_stack = .true.
      else
         print *,'Warning: Invalid keyword ',keyword(1:n),' found on ',
     *     'config file for SIAL program ',sial_program
         call abort_job() 
      endif

      go to 120   ! try reading another param line for this program.

  200 continue
      close (iunit)
      return

  500 continue
      
c--------------------------------------------------------------------------
c   Check for existence of $ACES_EXE_PATH/sial_config
c--------------------------------------------------------------------------

      call c_getenv('ACES_EXE_PATH'//char(0), envvar, lenenv, ierr)
      if (ierr .eq. 0) then
         config_file = envvar(1:lenenv) // '/sio/sial_config'
         open (unit=iunit, file = config_file, status = 'OLD',
     *      err=600, iostat = ios)

         go to 50   ! successfully opened.
      endif 

  600 continue
      return
      end
