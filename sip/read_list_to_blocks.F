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
      subroutine read_list_to_blocks(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table,  op)
c---------------------------------------------------------------------------
c   Wrapper routine for read_blocks_to_list instruction.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
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

      if (mpi_io_support) then
         call read_list_to_blocks_mpi_io(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table,  op)
      else
         call read_list_to_blocks_no_mpi_io(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table,  op)
      endif
     
      return
      end

      subroutine read_list_to_blocks_mpi_io(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table,  op)
c--------------------------------------------------------------------------
c   Reads block data from a list file and distributes them to their owners.
c--------------------------------------------------------------------------

      implicit none
      include 'aces3_version.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h' 
      include 'block_index.h'
      include 'dbugcom.h'
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

      integer i, k, ierr 
      integer iarray, array, type, ind, seg, iblock, iblock_map
      integer company, comm, company_comm
      integer niocompany , nbarray
      integer pst_get_company, pst_get_company_comm
      integer msg(len_sip_server_message)
      integer*8 msg8(len_sip_server_message/2)
      equivalence (msg, msg8)
      integer mx_mpi_req
      parameter (mx_mpi_req = 10000)
      integer request(mx_mpi_req)
      integer statuses(MPI_STATUS_SIZE, mx_mpi_req)
      integer mode, fh, fhndx
      integer*8 offset, size, data_file_loc, index_file_loc
      integer*8 bdata_table(2,narray_list)
      integer*8 indblk, get_block_index
      integer stack, blkndx, val1, val2, nindex
      integer status(MPI_STATUS_SIZE)

      logical served_data, file_exist
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif

      if (aces3_version .eq. '3.0.3') then

c---------------------------------------------------------------------------
c   Check for presence of BLOCK_INDEX_FORMAT file.  If it exists, then
c   the BLOCK_INDEX and BLOCKDATA files were created by the new algorithm
c   for write_blocks_to_list that uses MPI_IO, and we must read the files
c   using the new algorithm as well.  If the file does not exist, then the
c   BLOCK_INDEX and BLOCKDATA were created with the previous version, and
c   during this switchover, we will revert back to the non-MPI_IO version.
c---------------------------------------------------------------------------

         inquire (file='BLOCK_INDEX_FORMAT', exist=file_exist)
         if (.not. file_exist) then
            if (me .eq. 0) print *,'*** Restarting using non-parallel',
     *            'I/O algorithm (old format BLOCKDATA file.'
            call read_list_to_blocks_no_mpi_io(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table,  op)
            return
         endif
      endif
 
      data_file_loc  = 0
      index_file_loc = 0
      served_data = .false.
      do iarray = 1, narray_list
         array = array_list(iarray)
         type  = array_table(c_array_type, array)
         if (type .ne. served_array .and.
     *       type .ne. distributed_array) then
            print *,'Error: READ_LIST_TO_BLOCKS array ',
     *          array,' should be either served or distributed.'
            call abort_job()
         endif  

         if (type .eq. served_array) served_data = .true. 

         call determine_array_size(array, array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      size)
         nbarray = array_table(c_numblks, array)

         bdata_table(1,iarray) = index_file_loc
         bdata_table(2,iarray) = data_file_loc
         data_file_loc  = data_file_loc + size
         index_file_loc = index_file_loc + nbarray * lblk_index_entry
      enddo

      company_comm = pst_get_company_comm(me)
      if (served_data) then

c-------------------------------------------------------------------------
c   Send a message to each server indicating that the BLOCKDATA file should
c   be opened.
c--------------------------------------------------------------------------

         comm = mpi_comm_world

c--------------------------------------------------------------------------
c   Send a blocks_to_list "open" message (i. e. array = 0) to each server.
c--------------------------------------------------------------------------

         if (me .eq. 0) then
            msg(1) = sip_server_list_to_blocks
            msg(2) = 0         ! signal to open the files.
            msg(3) = 0         ! nindex
            niocompany = 0
            do i = 1, nprocs
               company = pst_get_company(i-1)
               if (company .eq. io_company_id) then
                  niocompany = niocompany + 1
                  call mpi_isend(msg, len_sip_server_message,
     *                     mpi_integer,
     *                     i-1, sip_server_message,
     *                     mpi_comm_world,
     *                     request(niocompany), ierr)
               endif
            enddo
         endif   ! me .eq. 0
      else
         comm = company_comm
      endif

      mode = MPI_MODE_RDONLY
      call mpi_file_open(comm, 'BLOCKDATA', mode,
     *                   MPI_INFO_NULL, fh, ierr)
      offset = 0
      call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION,
     *                       MPI_DOUBLE_PRECISION, 'native',
     *                       MPI_INFO_NULL, ierr)

      call mpi_file_open(comm, 'BLOCK_INDEX', mode,
     *                   MPI_INFO_NULL, fhndx, ierr)
      offset = 0
      call mpi_file_set_view(fhndx, offset, MPI_INTEGER,
     *                       MPI_INTEGER, 'native',
     *                       MPI_INFO_NULL, ierr)

c-----------------------------------------------------------------------------
c   Wait for completion of the server messages.
c-----------------------------------------------------------------------------

      if (me .eq. 0 .and. served_data) then
         call mpi_waitall(niocompany, request, statuses, ierr)
      endif

c---------------------------------------------------------------------------
c   Process each array.
c---------------------------------------------------------------------------

      do iarray = 1, narray_list
         array      = array_list(iarray)
         type       = array_table(c_array_type,array)
         nindex     = array_table(c_nindex, array)
         iblock_map = array_table(c_block_map, array)
         index_file_loc = bdata_table(1,iarray)
         data_file_loc  = bdata_table(2,iarray)
         if (type .eq. served_array) then

c---------------------------------------------------------------------------
c   Send a list_to_blocks message to each server.
c---------------------------------------------------------------------------

            msg(1) = sip_server_list_to_blocks
            msg(2) = array
            msg(3) = 0         ! nindex
            msg8(2) = index_file_loc
            msg8(3) = data_file_loc
            msg8(4) = array_table(c_numblks,array)
            niocompany = 0
            do i = 1, nprocs
               company = pst_get_company(i-1)
               if (company .eq. io_company_id) then
                  if (me .eq. 0) then
                     niocompany = niocompany + 1
                     call mpi_isend(msg, len_sip_server_message,
     *                  mpi_integer,
     *                  i-1, sip_server_message,
     *                  mpi_comm_world,
     *                  request(niocompany), ierr)
                  endif
               endif
            enddo

            if (me .eq. 0) then

c--------------------------------------------------------------------------
c   Wait for completion of server messages.
c--------------------------------------------------------------------------

               call mpi_waitall(niocompany, request,
     *                          statuses, ierr)
            endif   ! me .eq. 0 
         else

c---------------------------------------------------------------------------
c   Process the blocks of a distributed array.
c---------------------------------------------------------------------------

            offset = data_file_loc
            nbarray = array_table(c_numblks, array)
            do iblock = 1, nbarray

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

               blkndx = block_map_table(c_bmap_blkndx,
     *                                  iblock_map+iblock-1)
               if (block_map_table(c_processor,iblock_map+iblock-1)
     *                .eq. my_company_rank) then

c----------------------------------------------------------------------------
c   The block resides on this processor.   Read it from disk.
c----------------------------------------------------------------------------

                  stack  = array_table(c_array_stack,array)
                  indblk = get_block_index(array, iblock, stack,
     *                               blkndx, x, .true.)

                  call mpi_file_read_at(fh, offset, x(indblk),
     *                     blksize, MPI_DOUBLE_PRECISION, status,
     *                     ierr)
               endif

               offset = offset + blksize
            enddo    ! iblock
         endif   
      enddo   ! iarray

c-------------------------------------------------------------------------
c   Send servers the signal to close their files.
c-------------------------------------------------------------------------

      if (me .eq. 0 .and. served_data) then
         msg(1) = sip_server_list_to_blocks 
         msg(2) = list_to_blocks_done

         niocompany = 0
         do i = 1, nprocs
            company = pst_get_company(i-1)
            if (company .eq. io_company_id) then
               niocompany = niocompany + 1
               call mpi_isend(msg, len_sip_server_message,
     *                  mpi_integer,
     *                  i-1, sip_server_message,
     *                  mpi_comm_world,
     *                  request(niocompany), ierr)
            endif
         enddo

c--------------------------------------------------------------------------
c   Wait for completion of server messages.
c--------------------------------------------------------------------------

         call mpi_waitall(niocompany, request, statuses, ierr)
      endif   

c-------------------------------------------------------------------------
c   Close the BLOCKDATA file.
c-------------------------------------------------------------------------

      call mpi_file_close(fh, ierr)
      call mpi_file_close(fhndx, ierr)
      call mpi_barrier(company_comm, ierr)

      return
      end
