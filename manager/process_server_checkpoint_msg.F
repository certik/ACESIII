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
      subroutine process_server_checkpoint_msg(node, server_table,
     *                                   nserver_table)
c---------------------------------------------------------------------------
c   This subroutine handles the processing and manages state transitions
c   of a server_checkpoint_msg message.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'server_ckpt_data.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'dbugcom.h'
      include 'server_stat.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)
      integer node
      integer i, j, n, istart, nsearch, iblock
      integer*8 indblk, get_index_from_base
      integer*8 ind, imsg
      integer msgbuffer, state, ierr
      integer memloc, diskloc, size, ifile, nxt
      integer request, status(MPI_STATUS_SIZE)
      integer array, nindex, tag, f_form_msg_tag 
      integer msglen
      integer msg(100)

      integer nblk_file, header_entry, next
      integer*8 loc, index_file_diskloc, data_file_diskloc 
      integer*8 ndx_file_len, array_size
      integer str_trimlen
      logical created_file, match, flag

      integer index_file_entry_size
      parameter (index_file_entry_size = 3+mx_array_index*2)
      integer*8 index_file_entry(index_file_entry_size)

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

#ifdef ALTIX
      dptr = dshptr
#endif
      state = server_msg(c_msg_state,node)
      if (state .eq. begin_state) then
         array = server_msg(c_msg_array,node)

c---------------------------------------------------------------------------
c   Locate the server_table entries for this array.
c---------------------------------------------------------------------------

         do i = 1, nserved_arrays
            if (served_array_table(i) .eq. array) then
               istart = served_array_entry(i)
               nsearch = served_numblocks(i)
               go to 50
            endif
         enddo

c---------------------------------------------------------------------------
c   Array not found on this server.   We have nothing to do.
c---------------------------------------------------------------------------

	 server_msg(c_msg_state,node) = null_state
	 return
 
   50    continue
         indblk = get_index_from_base(base_mem_addr, x, 2)

c---------------------------------------------------------------------------
c   Calculate the index of the message buffer.
c--------------------------------------------------------------------------

         msgbuffer = server_msg(c_msg_msgbuffer,node)
         imsg = indblk + (msgbuffer-1)* server_mem_blocksize

c---------------------------------------------------------------------------
c   Open the checkpoint files.
c----------------------------------------------------------------------------

         n = str_trimlen(ckpt_ndx_filename)
         inquire (file = ckpt_ndx_filename(1:n), 
     *            exist=flag)
         if (flag) then
            call f_openfile(ckpt_ndx_filename(1:n) // char(0), 
     *                   ckpt_ndx_unit)
            created_file = .false.
         else
            call f_creatfile(ckpt_ndx_filename(1:n) // char(0), 
     *                        ckpt_ndx_unit)
            if (ckpt_ndx_unit .eq. -1) then
               print *,'Error: Server ',me,' cannot create checkpoint ',
     *                 'index file.'
               print *,'ckpt_ndx_filename = ',ckpt_ndx_filename
               call server_abort_job(server_table, nserver_table)
            else
               created_file = .true.
            endif
         endif 

         n = str_trimlen(ckpt_dat_filename)
         inquire (file = ckpt_dat_filename(1:n),
     *            exist=flag)
         if (flag) then
            call f_openfile(ckpt_dat_filename(1:n) // char(0), 
     *                   ckpt_dat_unit)
         else
            call f_creatfile(ckpt_dat_filename(1:n) // char(0),
     *                        ckpt_dat_unit)
            if (ckpt_dat_unit .eq. -1) then
               print *,'Error: Server ',me,' cannot create checkpoint ',
     *                 'data file.'
               print *,'ckpt_dat_filename = ',ckpt_dat_filename
               call server_abort_job(server_table, nserver_table)
            endif
         endif

c---------------------------------------------------------------------------
c   Determine the array size (in words) of all blocks of this array on the
c   current server.
c---------------------------------------------------------------------------

         array_size = 0
         do iblock = istart, istart+nsearch-1 
            array_size = array_size + 
     *           server_table(c_server_size,iblock)
         enddo

c---------------------------------------------------------------------------
c   Initialize the file header.
c---------------------------------------------------------------------------

         if (created_file) then
 
c---------------------------------------------------------------------------
c   File was opened for the first time.  Clear the file header.
c---------------------------------------------------------------------------

            do j = 1, mx_ckpt_arrays
               do i = 1, 4 
                  ndx_file_header(i,j) = 0
               enddo 

               data_file_sizes(j) = 0
            enddo 

            nblk_file          = 0
            header_entry       = 1
            data_file_diskloc  = 0
            ndx_file_header(1,header_entry) = array
            ndx_file_header(2,header_entry) = nsearch   ! num. of blocks
            ndx_file_header(3,header_entry) = -1
            ndx_file_header(4,header_entry) = -1
            data_file_sizes(header_entry)   = array_size
            data_file_location(header_entry) = -1

            commit_flag = .false.

c--------------------------------------------------------------------------
c   Initialize free space data structures.
c--------------------------------------------------------------------------

            nfree_space_ndx = 0
            nfree_space_dat = 0
            nfree_ndx_can = 0
            nfree_dat_can = 0
 
         else   ! created_file

c---------------------------------------------------------------------------
c   For each entry in the file header:
c   Word 1 = array
c   Word 2 = number of blocks checkpointed by this server.
c   Word 3 = disk pointer to index file data.
c   Word 4 = Number of array indices.
c---------------------------------------------------------------------------

            do i = 1, mx_ckpt_arrays
               if (ndx_file_header(1,i) .eq. array) then
                  header_entry = i   ! array was ckptd previously.
                  go to 200
               endif
 
               if (ndx_file_header(1,i) .eq. 0) then

c---------------------------------------------------------------------------
c   Array was never checkpointed.   Create a new header entry.
c---------------------------------------------------------------------------

                  header_entry = i

                  ndx_file_header(1,i) = array
                  ndx_file_header(2,i) = nsearch
                  ndx_file_header(3,i) = -1
                  ndx_file_header(4,i) = -1
                  data_file_location(i) = -1
                  data_file_sizes(i)   = array_size
                  go to 200
               endif 
            enddo
  200       continue

         endif   ! created_file

c---------------------------------------------------------------------------
c   Find locations on disk for the index file entries and data file entries.
c---------------------------------------------------------------------------

        call find_ckpt_disk_location(array, header_entry, 
     *                               array_size, nsearch) 
        index_file_diskloc = ndx_file_header(3,header_entry)
        data_file_diskloc  = data_file_location(header_entry) 

c---------------------------------------------------------------------------
c   Process each block of the array found in the server_table.
c---------------------------------------------------------------------------

         do iblock = istart, istart+nsearch-1
            if (server_table(c_server_array,iblock) .eq. array) then
               if (server_table(c_server_iblk,iblock) .eq. 
     *                               99999999) then
                  print *,'No message for iblock ',iblock
                  print *,(server_table(i,iblock),
     *                     i=1,lserver_table_entry)
                  call server_abort_job(server_table, nserver_table)
               endif

               if (server_table(c_server_memloc,iblock) .ne. 0) then
                  
c---------------------------------------------------------------------------
c   The data for this block is already in memory.
c---------------------------------------------------------------------------
   
                  j = server_table(c_server_memloc,iblock)
                  ind = indblk + (j-1) * server_mem_blocksize
                  size    = server_table(c_server_size,iblock)
               else if (server_table(c_server_diskloc,iblock) .gt. 
     *                      0) then 

c---------------------------------------------------------------------------
c   Read the data from disk into the node's message buffer.
c---------------------------------------------------------------------------

                  diskloc = server_table(c_server_diskloc,iblock)
                  size    = server_table(c_server_size,iblock)
                  ifile   = server_table(c_server_file, iblock)
                  call f_restoreram(server_unit(ifile), diskloc,
     *                           server_blocksizes(ifile),
     *                           x(imsg), size)
                  ind = imsg
               else
                  print *,'Server ',me,': Error in write_blocks_to_list'
                  print *,'   Array ',array,
     *                    ' erroneous server_table entry'
                  print *,(server_table(i,iblock),i=1,
     *                      lserver_table_entry)
                  call server_abort_job(server_table, nserver_table)
               endif
            endif

c-------------------------------------------------------------------------
c   Create the index file entry for this block.
c-------------------------------------------------------------------------

            index_file_entry(1) = server_table(c_server_iblk,iblock)
            index_file_entry(2) = server_table(c_server_size,iblock)
            index_file_entry(3) = data_file_diskloc
            next = 4
            do i = 1, mx_array_index
               index_file_entry(next) = 
     *                 server_table(c_server_bsegs+i-1,iblock)
               next = next + 1
            enddo

            do i = 1, mx_array_index
               index_file_entry(next) = 
     *                 server_table(c_server_esegs+i-1,iblock)
               next = next + 1
            enddo
 
c-------------------------------------------------------------------------
c   Write the index file entry for the block.
c-------------------------------------------------------------------------

            call f_write_disk(ckpt_ndx_unit, index_file_diskloc, 
     *                        index_file_entry, 
     *                        index_file_entry_size)
            index_file_diskloc = index_file_diskloc + 
     *                           index_file_entry_size

c-------------------------------------------------------------------------
c   Write the data for this block to the checkpoint file.
c-------------------------------------------------------------------------

            call f_write_disk(ckpt_dat_unit, data_file_diskloc,
     *                        x(ind), size)
            data_file_diskloc = data_file_diskloc + size

c--------------------------------------------------------------------------
c   Word 4 of file header entry contains the number of indices in the array.
c--------------------------------------------------------------------------

            ndx_file_header(4,header_entry) = 
     *                      server_table(c_server_nind,iblock) 
         enddo   ! iblock

 
c--------------------------------------------------------------------------
c   Close the checkpoint files.
c--------------------------------------------------------------------------

         call f_close_file(ckpt_ndx_unit)
         call f_close_file(ckpt_dat_unit)

c--------------------------------------------------------------------------
c   Send an acknowledgement to the master.  This tells him that the 
c   server has finished his checkpoint operation.
c   The message tag for the ack. msg is server_checkpoint_msg.
c--------------------------------------------------------------------------

            call mpi_isend(msg, 1,
     *            mpi_integer, 0, server_checkpoint_msg,
     *            mpi_comm_world,  server_msg(c_msg_request,node), 
     *            ierr)
            server_msg(c_msg_state,node) = wait_for_send_state
            return
      else if (state .eq. wait_for_send_state) then

c-------------------------------------------------------------------------
c   Test for completion of the send.
c-------------------------------------------------------------------------

         if (server_msg(c_msg_request,node) .eq.
     *                               MPI_REQUEST_NULL) then
            flag = .true.
         else
            call mpi_test(server_msg(c_msg_request,node), flag,
     *                 status, ierr)
         endif

         if (flag) then
               server_msg(c_msg_state,node) = null_state   ! done
         else
            return
         endif
      endif   ! state
      return
      end

      subroutine find_ckpt_disk_location(array, ientry, array_size, 
     *                                   array_nblks)
c----------------------------------------------------------------------------
c   Determines the locations on the index file and data file for an 
c   array.  Stores the resulting disk pointers in the ndx_file_header 
c   data structure, and manages the free space data structures.
c----------------------------------------------------------------------------
      implicit none
      include 'server_ckpt_data.h'
      include 'machine_types.h'
      include 'parallel_info.h'

      integer array, array_nblks
      integer*8 array_size

      integer i, j, ientry
      integer*8 index_file_diskloc, data_file_diskloc

      index_file_diskloc = ndx_file_header(3,ientry)
      data_file_diskloc  = data_file_location(ientry)

c---------------------------------------------------------------------------
c   Do we already have data in the ndx_file_header for this array?  If so 
c   (i. e. data_file_diskloc .gt. -1), we must move the existing entry to 
c   the list of free space candidates, as the data is old (from a previous
c   checkpoint).
c---------------------------------------------------------------------------

      if (data_file_diskloc .gt. -1) then
         nfree_ndx_can = nfree_ndx_can + 1
         if (nfree_ndx_can .gt. mx_ckpt_arrays) then
            print *,'Error: Out of free space candidates ',
     *                    'for checkpoint ndx file'
            call server_abort_job()
         endif

         free_space_ndx_candidate(1,nfree_ndx_can) =
     *                                 ndx_file_header(3,ientry)
         free_space_ndx_candidate(2,nfree_ndx_can) =
     *                          ndx_file_header(2,ientry)

         nfree_dat_can = nfree_dat_can + 1
         if (nfree_dat_can .gt. mx_ckpt_arrays) then
            print *,'Error: Out of free space candidates ',
     *              'for checkpoint dat file'
            call server_abort_job()
         endif

         free_space_dat_candidate(1,nfree_dat_can) =
     *                           ndx_file_header(4,ientry)
         free_space_dat_candidate(2,nfree_dat_can) =
     *                           data_file_sizes(ientry)
      endif

c----------------------------------------------------------------------------
c   Is there a block of free space available that is large enough?
c----------------------------------------------------------------------------

      do i = 1, nfree_space_ndx
         if (free_space_ndx(1,i) .ge. 0 .and.
     *       free_space_ndx(2,i) .ge. array_nblks) then

c-----------------------------------------------------------------------------
c   Use this block of free space for this array.
c-----------------------------------------------------------------------------

            ndx_file_header(2,ientry) = array_nblks
            ndx_file_header(3,ientry) = free_space_ndx(1,i)
            free_space_ndx(1,i) = -1
            free_space_ndx(2,i) = -1
            go to 1000 
         endif
      enddo

c----------------------------------------------------------------------------
c   Either there was no free space available, or there were no blocks 
c   large enough to hold this data.  Add to the end of file.
c----------------------------------------------------------------------------

      call get_filelen(ckpt_ndx_unit, index_file_diskloc)
      index_file_diskloc = index_file_diskloc / 8
      if (index_file_diskloc .eq. 0) 
     *      index_file_diskloc = index_file_header_size
      ndx_file_header(2,ientry) = array_nblks
      ndx_file_header(3,ientry) = index_file_diskloc
 
1000  continue

c----------------------------------------------------------------------------
c   Determine the data_file_diskloc.
c   Is there a block of free space available that is large enough?
c----------------------------------------------------------------------------

      do i = 1, nfree_space_dat
         if (free_space_dat(1,i) .ge. 0 .and.
     *       free_space_dat(2,i) .ge. array_size) then

c-----------------------------------------------------------------------------
c   Use this block of free space for this array.
c-----------------------------------------------------------------------------

            data_file_sizes(ientry) = array_size
            data_file_location(ientry) = free_space_dat(1,i)
            free_space_dat(1,i) = -1
            free_space_dat(2,i) = -1
            go to 2000 
         endif
      enddo

c----------------------------------------------------------------------------
c   Either there was no free space available, or there were no blocks 
c   large enough to hold this data.  Add to the end of file.
c----------------------------------------------------------------------------

      call get_filelen(ckpt_dat_unit, data_file_diskloc)
      data_file_diskloc = data_file_diskloc / bytes_per_double
      data_file_sizes(ientry) = array_size
      data_file_location(ientry) = data_file_diskloc
 
 2000 continue

      return
      end

