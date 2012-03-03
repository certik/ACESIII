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
      subroutine sip_barrier(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'dbugcom.h'

      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)

      integer status(MPI_STATUS_SIZE)
      integer request
      integer msg(5)
      integer i
      integer ierr
      integer company_comm, pst_get_company_comm
      integer company_rank, pst_get_company_rank
      integer nmessages, nmessages_global
      integer query_message_counter
      integer bmaster
      integer get_barrier_master
    
      if (dbg) then
         print *, 'Task ',me,
     *   ' is waiting at SIP barrier:  line number ',
     *   current_line
         call c_flush_stdout()
      endif

      if (my_company_size .eq. 1) return
      company_comm = pst_get_company_comm(me)
      company_rank = pst_get_company_rank(me)

      bmaster = get_barrier_master(company_rank)
      if (company_rank .ne. bmaster) then

c---------------------------------------------------------------------------
c   Send a barrier message to the master.
c---------------------------------------------------------------------------

         msg(1) = -7   ! barrier signal is -7, send with tag 3456.
         call mpi_isend(msg, 5, MPI_INTEGER, bmaster,
     *                  3456, company_comm, request, ierr)
      endif

c---------------------------------------------------------------------------
c   Subroutine exec_server_threads blocks until there is no thread on this
c   processor that is actively engaged in communication.
c---------------------------------------------------------------------------

      call exec_thread_server(1)
      if (dbg) print *,'Task ',me,' After exec_thread_server line ',
     *              current_line

c---------------------------------------------------------------------------
c   Wait for the request to be satisfied (should be done at this point).
c---------------------------------------------------------------------------

      if (company_rank .ne. bmaster) then
         call mpi_wait(request, status, ierr)
         if (dbg) print *,'Task ',me,
     *                    ' SIP_BARRIER: After mpi_wait, line ',
     *                     current_line
      endif

  100 continue
 
c---------------------------------------------------------------------------
c   Get the collective number of sent messages (i. e. GETs/PUTs) and the
c   number of received messages.
c---------------------------------------------------------------------------

      nmessages =  query_message_counter()   ! number on this processor
      call mpi_reduce(nmessages, nmessages_global, 1, MPI_INTEGER, 
     *                MPI_SUM, 0, company_comm, ierr)
      call mpi_bcast(nmessages_global, 1, MPI_INTEGER, 0,
     *               company_comm, ierr)
      if (dbg) print *,'Task ',me,' After exec_thread_server line ',
     *              current_line ,' nmessages, nmessages_global ',
     *              nmessages, nmessages_global
      if (nmessages_global .ne. 0) then
         call exec_thread_server(2)   ! checks for more messages
         go to 100
      endif

c---------------------------------------------------------------------------
c   Reset all persistent blocks  of all arrays which have had a "PUT' since
c   the last barrier.
c   This is done to insure data integrity.
c---------------------------------------------------------------------------

      do i = 1, narray_table
         if (array_table(c_put_flag,i) .ne. 0)  then
            call free_persistent_blocks(i, array_table, narray_table,
     *               index_table, nindex_table, block_map_table)
            array_table(c_put_flag,i) = 0
         endif
      enddo
     
      if (dbg) then
         print *,'Task ',me,' After free loop at line ',
     *      current_line
         call c_flush_stdout()
      endif 
      
c---------------------------------------------------------------------------
c   Reset the message counter in the thread server code.
c---------------------------------------------------------------------------

      call clear_message_counter()

c---------------------------------------------------------------------------
c   The final barrier guarantees that all processors have achieved a
c   quiet state with their threads.
c---------------------------------------------------------------------------

       call mpi_barrier(company_comm, ierr)
       if (dbg) then
          print *,'Task ',me,' Passed sip_barrier'
          call c_flush_stdout()
       endif

      return
      end
