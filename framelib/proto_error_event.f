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
      subroutine proto_error_event(msg)
c--------------------------------------------------------------------------
c   Prints a message to stdout, and sends the master task an error
c   event, tirggering job termination.
c--------------------------------------------------------------------------
      implicit none

      include 'mpif.h'
      include 'proto_events.h'

      character*(*) msg

      integer master, me, ierr
      integer pst_get_master
      integer event 
      save event

      call mpi_comm_rank(mpi_comm_world, me, ierr)
      print *,'ERROR on task ',me,': ',msg
      call c_flush_stdout()
 
      master = pst_get_master()
      event = mstr_error_event
      call mpi_send(event, 1, mpi_integer,
     *              master, master_event_tag, mpi_comm_world, ierr)
 
      return
      end
