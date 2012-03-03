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
      subroutine write_blocks_to_list(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c---------------------------------------------------------------------------
c   Wrapper routine for write_blcoks_to_list instruction.
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
         call write_blocks_to_list_mpi_io(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
      else
         call write_blocks_to_list_no_mpi_io(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, op)
      endif

      return
      end

      subroutine write_blocks_to_list_mpi_io(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Dumps all blocks of an array to a list.
c--------------------------------------------------------------------------

      implicit none
      include 'aces3_version.h'
      include 'interpreter.h'
      include 'int_gen_parms.h'
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

      integer*8 index_file_loc, data_file_loc, offset
      integer*8 size
      integer i, company, niocompany, nbarray
      integer array, type, iarray
      integer fh, fhndx, mode, ierr
      integer*8 bdata_table(2,narray_list)
      integer comm, company_comm, pst_get_company_comm, pst_get_company
      integer msg(len_sip_server_message)
      integer*8 msg8(len_sip_server_message/2)
      equivalence (msg, msg8)

      integer*8 indblk, get_block_index
      integer blkndx, val1, val2, stack, iblock_map, iblock, seg, ind, k
      integer nindex, nreq
      integer bseg(mx_array_index)
      integer eseg(mx_array_index)

      integer mx_mpi_req
      parameter (mx_mpi_req = 10000)
      integer request(mx_mpi_req)
      integer statuses(MPI_STATUS_SIZE, mx_mpi_req)
      integer status(MPI_STATUS_SIZE)

      integer block_index_entry(lblk_index_entry)

      logical served_data

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      if (aces3_version .eq. '3.0.3') then
    
c--------------------------------------------------------------------------
c   Write BLOCK_INDEX_FORMAT file, used to check for proper format on
c   restart.
c--------------------------------------------------------------------------

         open (unit=29, file='BLOCK_INDEX_FORMAT')
         write (29,*) aces3_version
         close(29)
      endif

c---------------------------------------------------------------------------
c   Pre-calculate each array's position in the BLOCKDATA file.
c---------------------------------------------------------------------------

#ifdef ALTIX
      dptr = dshptr
#endif
      data_file_loc  = 0
      index_file_loc = 0
      served_data = .false.
      company_comm = pst_get_company_comm(me)

      do iarray = 1, narray_list
         array = array_list(iarray)
         type  = array_table(c_array_type,array)
         if (type .ne. distributed_array .and.
     *       type .ne. served_array) then
            print *,'Error: BLOCKS_TO_LIST array must be either ',
     *        'a distributed or served array.'
            print *,'Array ',array,' has type = ',type
            call abort_job()
         endif

         if (type .eq. served_array) served_data = .true.
         nbarray = array_table(c_numblks,array) 
         call determine_array_size(array, array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      size)

         bdata_table(1, iarray) = index_file_loc
         bdata_table(2, iarray) = data_file_loc
         index_file_loc  = index_file_loc + nbarray * lblk_index_entry
         data_file_loc  = data_file_loc + size
      enddo

c-----------------------------------------------------------------------------
c   Open the files.
c-----------------------------------------------------------------------------

      if (served_data) then
         comm = mpi_comm_world

c--------------------------------------------------------------------------
c   Send a blocks_to_list "open" message (i. e. array = 0) to each server.
c--------------------------------------------------------------------------

         if (me .eq. 0) then
            msg(1) = sip_server_blocks_to_list
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
         comm = pst_get_company_comm(me)
      endif

      mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
      call mpi_file_open(comm, 'BLOCKDATA2', mode,
     *                   MPI_INFO_NULL, fh, ierr)
      offset = 0
      call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION,
     *                       MPI_DOUBLE_PRECISION, 'native',
     *                       MPI_INFO_NULL, ierr)

      call mpi_file_open(comm, 'BLOCK_INDEX2', mode,
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
         array          = array_list(iarray)
         type           = array_table(c_array_type,array)
         nindex         = array_table(c_nindex,array)
         iblock_map     = array_table(c_block_map, array)
         index_file_loc = bdata_table(1,iarray)
         data_file_loc  = bdata_table(2,iarray)
 
         if (type .eq. served_array) then
            
c---------------------------------------------------------------------------
c   Send a blocks_to_list message to each server.
c---------------------------------------------------------------------------

            msg(1) = sip_server_blocks_to_list
            msg(2) = array
            msg(3) = 0         ! nindex
            niocompany = 0
            do i = 1, nprocs
               company = pst_get_company(i-1)
               if (company .eq. io_company_id) then

                  msg8(2) = index_file_loc
                  msg8(3) = data_file_loc

                  call determine_server_contribution(array, i-1,
     *                      array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      nbarray, size)

                  if (me .eq. 0) then 
                     niocompany = niocompany + 1
                     call mpi_isend(msg, len_sip_server_message, 
     *                  mpi_integer,
     *                  i-1, sip_server_message,
     *                  mpi_comm_world,
     *                  request(niocompany), ierr)
                  endif

                  data_file_loc  = data_file_loc + size
                  index_file_loc = index_file_loc + 
     *                             nbarray * lblk_index_entry
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

            nreq = 0
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
                  bseg(k) = val1
                  eseg(k) = val2
                  blksize = blksize * (val2-val1+1)
               enddo

               blkndx = block_map_table(c_bmap_blkndx,
     *                                  iblock_map+iblock-1)
               if (block_map_table(c_processor,iblock_map+iblock-1) 
     *                .eq. my_company_rank) then

c----------------------------------------------------------------------------
c   The block resides on this processor.   Write it to disk.
c----------------------------------------------------------------------------

                  stack  = array_table(c_array_stack,array)
                  indblk = get_block_index(array, iblock, stack, 
     *                               blkndx, x, .true.)

                  if (nreq .eq. mx_mpi_req) then
                     call mpi_waitall(nreq, request, statuses, ierr)
                     nreq = 0
                  endif

                  nreq = nreq + 1
                  call mpi_file_iwrite_at(fh, offset, x(indblk),
     *                     blksize, MPI_DOUBLE_PRECISION, request(nreq),
     *                     ierr)

c---------------------------------------------------------------------------
c   Create the BLOCK_INDEX entry for this block.
c---------------------------------------------------------------------------
    
                  block_index_entry(c_blk_index_array) = array
                  block_index_entry(c_blk_index_blkno) = iblock
                  block_index_entry(c_blk_index_size)  = blksize
                  block_index_entry(c_blk_index_nind)  = nindex
                  do k = 1, nindex
                     block_index_entry(c_blk_index_bsegs+k-1) = bseg(k)
                     block_index_entry(c_blk_index_esegs+k-1) = eseg(k)
                  enddo
 
c---------------------------------------------------------------------------
c   Write the BLOCK_INDEX entry.
c---------------------------------------------------------------------------

                  call mpi_file_write_at(fhndx, index_file_loc, 
     *                  block_index_entry, lblk_index_entry,
     *                  MPI_INTEGER, status,
     *                  ierr)
               endif

               offset = offset + blksize
               index_file_loc = index_file_loc + lblk_index_entry
            enddo    ! iblock

c--------------------------------------------------------------------------
c   Wait for completion of any dist. data blocks.
c---------------------------------------------------------------------------

            if (nreq .gt. 0) call mpi_waitall(nreq, request, 
     *                                        statuses, ierr)
         endif
      enddo

c---------------------------------------------------------------------------
c   Send servers the signal to close their files.
c---------------------------------------------------------------------------

      if (me .eq. 0 .and. served_data) then
         msg(1) = sip_server_blocks_to_list
         msg(2) = blocks_list_done

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

c--------------------------------------------------------------------------
c   Sync the data and close the files.
c   Wait until all in the worker's company have finished processing.
c--------------------------------------------------------------------------

      call mpi_barrier(company_comm, ierr)
      call mpi_file_sync(fh, ierr)
      call mpi_file_close(fh, ierr)
      call mpi_file_sync(fhndx, ierr)
      call mpi_file_close(fhndx, ierr)

c--------------------------------------------------------------------------
c   Rename the files to their actual values.
c--------------------------------------------------------------------------

      if (me .eq. 0) then
         call f_renamefile('BLOCKDATA2'//char(0), 'BLOCKDATA'//char(0))
         call f_renamefile('BLOCK_INDEX2'//char(0), 
     *                      'BLOCK_INDEX'//char(0))
         
      endif
      return
      end

      subroutine determine_server_contribution(array, proc,
     *                      array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      nb, size)
c---------------------------------------------------------------------------
c   Calculates the number of blocks and number of words contained in all
c   the blocks of an array residing on server "proc".
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'

      integer array, proc
      integer narray_table, nindex_table, nsegment_table,
     *        nblock_map_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nb
      integer*8 size

      integer nindex, nblk, iblock_map, iblock, i,k, blksize
      integer ind, seg, val1, val2
      nb = 0
      size = 0

      nindex     = array_table(c_nindex,array)
      nblk       = array_table(c_numblks, array)
      iblock_map = array_table(c_block_map, array)

      do iblock = 1, nblk
         if (block_map_table(c_processor,iblock_map+iblock-1) 
     *                                    .eq. proc) then
            nb = nb + 1

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
         endif
      enddo

      return
      end

      double precision function block_energy(x, n)
      integer n
      double precision x(1)

      integer i
      double precision sum

      sum = 0.
      do i = 1, n
         sum = sum + x(i) * x(i)
      enddo

      block_energy = sum
      return
      end

