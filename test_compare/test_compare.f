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
      program test_compare
c-------------------------------------------------------------------------
c   Program to compare test results with that found in a standardized
c   database of results.
c-------------------------------------------------------------------------

      implicit none
      include 'mpif.h'

      integer ios, lenenv, ierr
      integer argc, f_iargc
      integer i, n, nw, ntitle, nkey
      integer str_trimlen 
      integer me, mpierr
      double precision jobarc_vals(1000)
      double precision rvalue, rmserr

      character*256 db_fname, envvar
      character*20 test_name
      character*80 title, line, key

c--------------------------------------------------------------------------
c   Initialize mpi.
c--------------------------------------------------------------------------

      call mpi_init(mpierr)
      call mpi_comm_rank(mpi_comm_world, me, mpierr)
      if (me .ne. 0) go to 500

c--------------------------------------------------------------------------
c   Initialize ACES2 file system, open JOBARC.
c--------------------------------------------------------------------------

      call init_machine_types()
      call aces_init_rte()
      call aces_ja_init()

c--------------------------------------------------------------------------
c   Get ACES_EXE_PATH environment variable.  Test database should be in 
c   this directory.
c--------------------------------------------------------------------------

      call c_getenv('ACES_EXE_PATH'//char(0), envvar, lenenv, ierr)
      if (ierr .ne. 0) then
         print *,'Error: ACES_EXE_PATH env. variable was not found.'
         stop
      endif

      db_fname = envvar(1:lenenv) // '/' // 'test_results'

c--------------------------------------------------------------------------
c   Open database of test results.
c--------------------------------------------------------------------------

      open (unit=20,file=db_fname,status='OLD',err=1000,iostat=ios)

c---------------------------------------------------------------------------
c   Get command-line argument.  This tells us which test we are doing.
c---------------------------------------------------------------------------

      argc = f_iargc()
      call f_getarg(1, test_name)
      n = str_trimlen(test_name)

c----------------------------------------------------------------------------
c   Search the database for this test.
c----------------------------------------------------------------------------

   10 continue
      read (unit=20, fmt=50, end=100 ) line
      if (line(1:n) .ne. test_name(1:n)) go to 10

c----------------------------------------------------------------------------
c   Test was found.  Read the next line. 
c----------------------------------------------------------------------------

      ntitle = str_trimlen(line) 
      title = line (1:ntitle)
      read (unit=20,fmt=50) line
   50 format(a80)

      do i = 1, 80
         if (line(i:i) .eq. ' ') then
            nkey = i-1
            read (line(i+1:80),'(i3)') nw
            go to 70
         else
            key(i:i) = line(i:i)
         endif
      enddo
  70  continue 

      call dgetrec(1,'JOBARC', key(1:nkey), nw, jobarc_vals)

c----------------------------------------------------------------------------
c   Calculate the RMS error between the standardized values and those 
c   found on JOBARC.
c----------------------------------------------------------------------------

      rmserr = 0.d0
      do i = 1, nw
         read (unit=20,fmt=9000) rvalue
         rmserr = rmserr + (rvalue-jobarc_vals(i)) * 
     *                     (rvalue-jobarc_vals(i))  
      enddo 

      rmserr = dsqrt(rmserr) / dble(nw)

c---------------------------------------------------------------------------
c   Print the test result.
c---------------------------------------------------------------------------

      print *,'TEST ',title(1:ntitle),' RMS error ',rmserr       
  100 continue
      call aces_ja_fin() 

  500 continue
      call mpi_barrier(mpi_comm_world, mpierr)
      call mpi_finalize(mpierr)
      stop

 1000 continue
      print *,'Error: Cannot open test database: status = ',ios
      print *,'Test database filename ',db_fname
      stop
 9000 format(g20.12)
      end
