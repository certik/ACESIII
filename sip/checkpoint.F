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
      subroutine checkpoint(array_table, narray_table,
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
      include 'int_gen_parms.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'

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
      integer msg(5)
      integer i
      integer ierr
      integer company_comm, pst_get_company_comm
      integer array, array_type

      integer pst_get_company
      integer niocompany
      integer request(10000)
      integer statuses(MPI_STATUS_SIZE,10000)

      character*40 ckpt_msg

      if (.not. mpi_io_support) return

      array = op(c_result_array)
      if (array .ne. 0) then
         array_type = array_table(c_array_type, array)
         if (array_type .ne. served_array .and. 
     *       array_type .ne. distributed_array) then
            print *,
     *       'Error: checkpoint instruction requires a served or ',
     *       'distributed array'
            call abort_job()
         endif
      else
         array_type = 0   ! no arg for checkpoint instruction
      endif
 
      if (me .eq. 0) then

c-------------------------------------------------------------------------
c   Send a checkpoint message to each I/O server.
c-------------------------------------------------------------------------

         if (array_type .eq. served_array) then
            niocompany = 0
            msg(1) = sip_server_checkpoint_signal
            msg(2) = array
            msg(4) = current_line   ! in the "tag" field of the msg.
            do i = 1, nprocs
               if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Send the checkpoint message.
c---------------------------------------------------------------------------

                  niocompany = niocompany + 1
                  call mpi_isend(msg, 4,
     *               MPI_INTEGER, i-1,
     *               sip_server_message, mpi_comm_world,
     *               request(niocompany), ierr)
               endif
            enddo

c---------------------------------------------------------------------------
c   Wait until all servers have received their checkpoint signal.
c---------------------------------------------------------------------------

            call mpi_waitall(niocompany, request,
     *                    statuses, ierr)
            call prt_time('All server ckpt signals were received.')

            niocompany = 0
            do i = 1, nprocs
               if (pst_get_company(i-1) .eq. io_company_id) then

c---------------------------------------------------------------------------
c   Proc i-1 is a server.  Receive the "checkpoint complete" message.
c   Each server will send a 1-word indicator message using the
c   sip_server_checkpoint_signal as the message tag.
c---------------------------------------------------------------------------

                  niocompany = niocompany + 1
                  call mpi_irecv(msg, 1,
     *               MPI_INTEGER, i-1,
     *               sip_server_checkpoint_signal, mpi_comm_world,
     *               request(niocompany), ierr)
               endif
            enddo

c-----------------------------------------------------------------------
c   Wait for all the acknowledgement messages.
c------------------------------------------------------------------------

            call mpi_waitall(niocompany, request, statuses,
     *                    ierr)
         endif   ! array_type .eq. served_array
      endif   ! me .eq. 0


      if (array_type .eq. distributed_array) then
         call checkpoint_distributed_array(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
      endif

c--------------------------------------------------------------------------
c   Workers must wait at this barrier untile all have completed their task.
c--------------------------------------------------------------------------

      company_comm = pst_get_company_comm(me) 
      call mpi_barrier(company_comm, ierr)
      if (dbg) then
         print *,'Task ',me,' Passed checkpoint worker barrier'
         call c_flush_stdout()
      endif

      return
      end

      subroutine checkpoint_distributed_array(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
c---------------------------------------------------------------------------
c   Performs a checkpoint of a distributed array in parallel using MPI-IO.
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'int_gen_parms.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

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

      integer array
      integer company_comm, pst_get_company_comm
      integer fh, mode, ierr
      integer nblk, iblock, iblock_map
      integer ind(mx_array_index)
      integer seg(mx_array_index)
      integer i, k, size , nreq
      integer*8 offset

      integer mx_mpi_req
      parameter (mx_mpi_req = 10000)
      integer request(mx_mpi_req)
      integer statuses(MPI_STATUS_SIZE,mx_mpi_req) 
      integer blkndx, stack, nindex, val1, val2
      integer ientry
      integer*8 filelen
      integer*8 indblk, get_block_index

      logical new_array

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      array = op(c_result_array)
      nindex = array_table(c_nindex,array)
#ifdef ALTIX
      dptr = dshptr
#endif

c---------------------------------------------------------------------------
c   Determine if this array has previously been checkpointed.
c---------------------------------------------------------------------------

      new_array = .false.
      do i = 1, nckpt_arrays
         if (ckpt_arrays(i) .eq. array) then
            ientry = i
            go to 100
         endif
      enddo

c---------------------------------------------------------------------------
c   Add array to table of checkpointed arrays.
c---------------------------------------------------------------------------

      nckpt_arrays          = nckpt_arrays + 1
      ientry                = nckpt_arrays
      ckpt_arrays(ientry)   = array
      ckpt_diskaddr(ientry) = 0
      call determine_array_size(array, array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      disk_sizes(ientry))

      new_array             = .true.
 
  100 continue

c---------------------------------------------------------------------------
c   Open the checkpoint data file.
c---------------------------------------------------------------------------

      company_comm = pst_get_company_comm(me)
      mode = MPI_MODE_CREATE + MPI_MODE_WRONLY 
      call mpi_file_open(company_comm, 'CKPT_DIST', mode,
     *                   MPI_INFO_NULL, fh, ierr)
      offset = 0
      call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION,
     *                       MPI_DOUBLE_PRECISION, 'native', 
     *                       MPI_INFO_NULL, ierr)

      call mpi_file_get_size(fh, filelen, ierr)
      filelen = filelen / bytes_per_double
      call mpi_barrier(company_comm, ierr)

c---------------------------------------------------------------------------
c   Determine beginning array file offset to use for this array.
c---------------------------------------------------------------------------

      if (new_array) then

c---------------------------------------------------------------------------
c   Can we find free space to hold the array?
c---------------------------------------------------------------------------

         do i = 1, nfree_space
            if (free_space_sizes(i) .ge. disk_sizes(ientry)) then
               ckpt_diskaddr(ientry) = free_space(i)
               free_space(i) = -1
               free_space_sizes(i) = -1
               go to 150
            endif
         enddo
 
c---------------------------------------------------------------------------
c   No free space available is large enough.  Add to end of file.
c---------------------------------------------------------------------------

         ckpt_diskaddr(ientry) = filelen   
  150    continue
         offset = ckpt_diskaddr(ientry)  
      else
         
c---------------------------------------------------------------------------
c   Check for a free space entry large enough to hold the array.
c---------------------------------------------------------------------------

         do i = 1, nfree_space
            if (free_space_sizes(i) .ge. disk_sizes(ientry)) then
               offset = free_space(i)

c---------------------------------------------------------------------------
c   Add existing data entry to list of free space candidates.
c---------------------------------------------------------------------------

               nfree_space_candidates = nfree_space_candidates + 1
               free_space_candidates(nfree_space_candidates) = 
     *                          ckpt_diskaddr(ientry)
               candidate_sizes(nfree_space_candidates) = 
     *                          disk_sizes(ientry)

c---------------------------------------------------------------------------
c   Mark existing free space entry as used.
c---------------------------------------------------------------------------

               free_space(i)       = -1
               free_space_sizes(i) = -1
               ckpt_diskaddr(ientry) = offset
               go to 200
            endif
         enddo

c---------------------------------------------------------------------------
c   No free space entry was found that is large enough.  Add to the end of 
c   file, and use the existing copy of the array as a new free space 
c   candidate.  When the checkpoint is committed, the candidates become
c   actual free space entries.
c---------------------------------------------------------------------------

         nfree_space_candidates = nfree_space_candidates + 1
         if (nfree_space .gt. mx_ckpt_arrays) then
            print *,'Error: Checkpoint has run out of room ',
     *          'for free space entries.'
            call abort_job()
         endif
 
         offset                        = filelen
         free_space_candidates(nfree_space_candidates) = 
     *                ckpt_diskaddr(ientry)
         candidate_sizes(nfree_space_candidates) = disk_sizes(ientry)

         ckpt_diskaddr(ientry)         = offset
  200    continue
      endif

c---------------------------------------------------------------------------
c   Loop over the blocks of the array, writing them into the appropriate slot 
c   within the file.
c----------------------------------------------------------------------------

      nblk = array_table(c_numblks, array)
      iblock_map = array_table(c_block_map, array)
      nreq = 0
      do iblock = 1, nblk

c--------------------------------------------------------------------------
c   Determine block size in words.
c--------------------------------------------------------------------------

         size = 1
         do k = 1, nindex
            ind(k) = array_table(c_index_array1+k-1,array)
            seg(k) = block_map_table(c_block_map_seg+k-1,
     *                       iblock_map+iblock-1)
            call get_index_segment(ind(k), seg(k), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            size = size * (val2-val1+1)
         enddo

         blkndx = block_map_table(c_bmap_blkndx,iblock_map+iblock-1)
         if (block_map_table(c_processor,iblock_map+iblock-1) .eq.
     *          my_company_rank) then

c----------------------------------------------------------------------------
c   The block resides on this processor.   Write it to disk.
c----------------------------------------------------------------------------

            stack  = array_table(c_array_stack,array)
            indblk = get_block_index(array, iblock, stack, blkndx, 
     *                               x, .true.)
 
            if (nreq .eq. mx_mpi_req) then
               call mpi_waitall(nreq, request, statuses, ierr)
               nreq = 0
            endif

            nreq = nreq + 1
            call mpi_file_iwrite_at(fh, offset, x(indblk), size, 
     *                              MPI_DOUBLE_PRECISION, request(nreq),
     *                              ierr)
         endif   

         offset = offset + size
      enddo    ! iblock
    
c----------------------------------------------------------------------------
c   Wait on MPI requests to complete.
c---------------------------------------------------------------------------

      if (nreq .gt. 0) call mpi_waitall(nreq, request, statuses, ierr)

c----------------------------------------------------------------------------
c   Sync the data to the file, close file.
c----------------------------------------------------------------------------

      call mpi_file_sync(fh, ierr)
      call mpi_file_close(fh, ierr)

      return
      end

      subroutine determine_array_size(array, array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      size)
c---------------------------------------------------------------------------
c   Determines the size of an array (in words).
c   Only works for distributed and served arrays.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h' 

      integer array
      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer*8 size

      integer k, nindex, nblk, iblock_map
      integer iblock, ind, seg
      integer blksize, val1, val2

      nindex     = array_table(c_nindex,array)
      nblk       = array_table(c_numblks, array)
      iblock_map = array_table(c_block_map, array)

      size = 0
      do iblock = 1, nblk

c--------------------------------------------------------------------------
c   Determine block size in words.
c--------------------------------------------------------------------------

         blksize = 1
         do k = 1, nindex
            ind = array_table(c_index_array1+k-1,array)
            seg = block_map_table(c_block_map_seg+k-1,
     *                       iblock_map+iblock-1)
            call get_index_segment(ind, seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1, val2)
            blksize = blksize * (val2-val1+1)
         enddo

         size = size + blksize 
      enddo   ! iblock

      return
      end
     
