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
      subroutine joda_exit(istat, string)
      implicit none
      include 'mpif.h'
      include 'parallel_info.h'
      
#ifdef ACES64
      integer*8 istat
#else
      integer istat
#endif
      integer ierr
      character*(*) string

c---------------------------------------------------------------------------
c   Perform xaces3 cleanup so we can stop cleanly if a joda-based routine
c   decides to execute a "STOP".
c---------------------------------------------------------------------------

      call c_flush_stdout()  ! flush joda print to stdout

c---------------------------------------------------------------------------
c   Signal the other processors that we are through.
c---------------------------------------------------------------------------

      call mpi_bcast(1, 1, mpi_integer, 0,
     *               mpi_comm_world, ierr)

c---------------------------------------------------------------------------
c   Remove the restart file.
c---------------------------------------------------------------------------

      print *,'Joda: stop with istat ',istat
      call prt_time(string)  
      close (restart_file)
      call f_deletefile('aces3.restart' // char(0), -1)

c---------------------------------------------------------------------------
c   Finalize MPI.
c---------------------------------------------------------------------------

      call mpi_finalize(ierr) 
      stop
      end
