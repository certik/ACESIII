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
      subroutine distribute_server_info(array_table, narray_table,
     *               index_table, nindex_table,
     *               segment_table, nsegment_table,
     *               block_map_table, nblock_map_table,
     *               iocompany_id, niocompany, scr, nscr,
     *               stack_blocksizes, nstacks )
c---------------------------------------------------------------------------
c   Builds the server table for each server, then sends it to the proper 
c   server.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'server_struct.h'
      include 'mpif.h'
      include 'dbugcom.h'

      integer narray_table, nindex_table, nsegment_table, 
     *        nblock_map_table, iocompany_id, niocompany, nscr
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer block_map_table(lblock_map_entry,nblock_map_table)
      integer scr(lserver_table_entry,*)
      integer nstacks
      integer stack_blocksizes(nstacks)

      integer i, j, array, nind, segment
      integer start, rank, iblock_map, nwsend, ierr, iserver, nprocs
      integer nblocks
      integer next_entry 
      integer server_rank(niocompany), request(niocompany)
      integer statuses(MPI_STATUS_SIZE,niocompany)
      integer index(mx_array_index)
      integer pst_get_company

      if (niocompany .eq. 0) return   ! no servers

c---------------------------------------------------------------------------
c   Determine the ranks of each server.
c---------------------------------------------------------------------------

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      next_entry = 0
      do i = 1, nprocs
         if (pst_get_company(i-1) .eq. iocompany_id) then
            next_entry = next_entry + 1
            server_rank(next_entry) = i-1
         endif
      enddo

      if (dbg) print *,'Server ranks: ',(server_rank(i),i=1,niocompany)

      next_entry = 0

      print *,'Number of block_map entries ',nblock_map_table

      do iserver = 1, niocompany

c---------------------------------------------------------------------------
c   Build the table for this server.
c---------------------------------------------------------------------------

         rank = server_rank(iserver)
         start = next_entry + 1

c---------------------------------------------------------------------------
c   Find each served array in the array table.
c---------------------------------------------------------------------------

         do array = 1, narray_table
            if (array_table(c_array_type,array) .eq. served_array) then
               iblock_map = array_table(c_block_map,array)
               nind       = array_table(c_nindex,array)
               nblocks    = array_table(c_numblks,array)
                
               do i = 1, nind
                  index(i) = array_table(c_index_array1+i-1,array)
               enddo

c--------------------------------------------------------------------------
c   Scan the block_map_table for the blocks needed by this server.
c--------------------------------------------------------------------------

               do i = 1, nblocks
                  if (block_map_table(c_processor,iblock_map+i-1) .eq.
     *                  rank) then

c--------------------------------------------------------------------------
c   Add an entry to the server table for this rank.
c--------------------------------------------------------------------------

                     next_entry = next_entry + 1
                     if (next_entry*lserver_table_entry .gt. nscr) then
                        print *,'Error: Not enough scratch space ',
     *                     'in distribute_server_info. nscr ',nscr,
     *                     ' but need at least ',
     *                     next_entry*lserver_table_entry
                        call abort_job()
                     endif

                     do j = 1, lserver_table_entry
                        scr(j,next_entry) = 0
                     enddo

                     scr(c_server_array,next_entry) = array
                     scr(c_server_nind,next_entry) = nind

                     do j = 1, nind
                        segment = block_map_table(c_block_map_seg+j-1,
     *                                            iblock_map+i-1)
                        call get_index_segment(index(j), segment,
     *                             segment_table, nsegment_table, 
     *                             index_table, nindex_table, 
     *                             scr(c_server_bsegs+j-1,next_entry),
     *                             scr(c_server_esegs+j-1,next_entry))
                     enddo
                  endif
               enddo    ! block loop
            endif
         enddo          ! array loop

c---------------------------------------------------------------------------
c   Send this table to its server.
c---------------------------------------------------------------------------

         nwsend = (next_entry-start+1)*lserver_table_entry
         call mpi_isend(scr(1,start), nwsend, MPI_INTEGER, rank,
     *                  sip_server_data_message, mpi_comm_world,
     *                  request(iserver), ierr)
      enddo             ! server loop

c--------------------------------------------------------------------------
c   Wait for the requests to be received by the servers.
c-------------------------------------------------------------------------

      call mpi_waitall(niocompany, request, statuses,
     *                 ierr)

c---------------------------------------------------------------------------
c   Send the blkmgr array of stacksizes to each server.
c---------------------------------------------------------------------------

      do iserver = 1, niocompany
         rank = server_rank(iserver)
         call mpi_isend(stack_blocksizes, nstacks, 
     *                  MPI_INTEGER, rank,
     *                  sip_server_data_message, mpi_comm_world,
     *                  request(iserver), ierr)
      enddo

c--------------------------------------------------------------------------
c   Wait for the requests to be received by the servers.
c-------------------------------------------------------------------------

      call mpi_waitall(niocompany, request, statuses,
     *                 ierr)
      return
      end
