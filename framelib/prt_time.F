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
      subroutine prt_time(msg)
c-----------------------------------------------------------------------
c   Prints the current time and date, with a text message.
c-----------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'dbugcom.h'

      character*(*) msg
      integer me, ierr
      character*8 dt
      character*10 tm

#ifdef HP
      character*24 ctime
      integer system_time
      integer time
      external time
      external ctime
#endif

      if (.not. dbg) return
      call mpi_comm_rank(mpi_comm_world, me, ierr)

#ifdef HP
      system_time = time()
      print *,'Task ',me,' ',ctime(system_time),' ',msg 
#else
      call date_and_time(dt, tm)

      print *,'TASK ',me,tm(1:2),':',tm(3:4),':',tm(5:6),' ',
     *     msg
#endif
      call c_flush_stdout()
      return
      end
