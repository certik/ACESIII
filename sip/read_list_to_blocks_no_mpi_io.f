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
      subroutine read_list_to_blocks_no_mpi_io(array_table,
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table,  op)
c--------------------------------------------------------------------------
c   Reads block data from a list file and distributes them to their owners.
c   This version does not use MPI_IO.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h' 
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

      integer i, j, k, n
      integer array, index(mx_array_index), nindex, ierr
      integer block, seg(mx_array_index)
      integer blkndx, get_block_number
      integer saveseg(mx_array_index)
      integer*8 indblk, get_block_data_index
      integer stack
      integer f_form_msg_tag
      integer datatag
      integer type, kreq, ireq, nreq
      
      integer desc_msg(2)

      integer my_comm_rank, comm
      integer home, iblk
      integer val1(mx_array_index), val2(mx_array_index)
      double precision x(1)
    
      integer nb, proc, map, company_comm
      integer tag
      integer pst_get_company_comm, pst_get_company_rank
      integer status(MPI_STATUS_SIZE)
      integer iarray

      company_comm    = pst_get_company_comm(me)
      proc  = pst_get_company_rank(me)
      if (proc .eq. 0) then
         call master_list_to_blocks(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, op)
         return
      endif

  100 continue

c---------------------------------------------------------------------------
c   Receive a block_list_array_tag message.  This contains the data for the 
c   next array to be processed.
c---------------------------------------------------------------------------

      if (dbg) print *,'Task ',me,
     *             ' LIST_TO_BLOCKS: wait for array data '

      call mpi_recv(array, 1, mpi_integer, 0,
     *                     block_list_array_tag,
     *                     company_comm, status, ierr)
      if (dbg) print *,'Task ',me,'   Received array message: array = ',
     *          array

c---------------------------------------------------------------------------
c   Check the array data and determine its type.
c---------------------------------------------------------------------------

      if (array .eq. -1) then          ! this is a "quit" signal
         if (dbg) print *,'LIST_TO_BLOCKS: task ',me,
     *               ' completed processing'
         narray_list = 0
         return 
      else if (array .lt. 1 .or. array .gt. narray_table) then  ! out of range
         print *,'Task ',me,' LIST_TO_BLOCKS: received array ',array,
     *      ' which is out of range.'
         call abort_job()
      else
         type = array_table(c_array_type,array)
         if (dbg) print *,'Task ',me,' LIST_TO_BLOCKS: begin array ',
     *         array,' type ',type 
      endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      if (type .eq. served_array) go to 100   ! master handles this.

      nb = array_table(c_numblks, array)
      map = array_table(c_block_map, array)

      nindex = array_table(c_nindex, array)

c---------------------------------------------------------------------------
c   Pick up the array's indices, and save the current segments of each.
c---------------------------------------------------------------------------

      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1,array)
         saveseg(i) = index_table(c_current_seg,index(i))
      enddo

  200 continue

c-------------------------------------------------------------------------
c   Receive a block_descriptor message.  If a descriptor contains a -1,
c   that indicates processing for this array is complete.
c-------------------------------------------------------------------------

      call mpi_recv(desc_msg, 2, MPI_INTEGER,
     *                    0, block_list_descriptor_tag,
     *                    company_comm, status,ierr)

c-------------------------------------------------------------------------
c   Unpack the descriptor message.
c-------------------------------------------------------------------------

         iblk    = desc_msg(1)
         datatag = desc_msg(2)

         if (datatag .eq. -1) then    ! thru with this array
            if (dbg) print *,'Task ',me,' LIST_TO_BLOCKS recvd ',
     *                  ' array complete msg for array ',array
            go to 100
         endif

c-------------------------------------------------------------------------
c   Make sure this is the block's "home" processor.
c-------------------------------------------------------------------------

         home = block_map_table(c_processor, map+iblk-1)
         if (home .ne. proc) then
            print *,'Error: array ',array,' block ',iblk,
     *          ' resides on ',
     *          'proc ',home,' but is being processed on ',me
            call abort_job()
         endif

c-------------------------------------------------------------------------
c   Find the block's memory location.
c-------------------------------------------------------------------------

         blkndx = get_block_number(array, iblk)
         stack  = array_table(c_array_stack,array)
         indblk = get_block_data_index(array, iblk, stack, blkndx, x)
         call get_actual_blocksize(array, iblk, blkndx,
     *              array_table, narray_table,
     *              index_table, nindex_table,
     *              segment_table, nsegment_table, n)

c--------------------------------------------------------------------------
c   Post a recv for the block.
c--------------------------------------------------------------------------

         call mpi_recv(x(indblk), n, mpi_double_precision, 0,
     *                     datatag,
     *                     company_comm, status, ierr)
      go to 200   ! process another block.

      return
      end

      subroutine master_list_to_blocks(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, op)
c--------------------------------------------------------------------------
c   Reads the blocks of an array from a list file and distributes them to 
c   their owners.
c   This routine handles the master's role.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'saved_data.h'
      include 'blockdata_list.h'
      include 'machine_types.h'
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

      integer i, j, k, n
      integer array, index(mx_array_index), nindex, ierr
      integer block, seg(mx_array_index)
      integer dblock, dnindex, dseg(mx_array_index)
      integer*8 iscr
      integer*8 indblk
      integer*8 get_block_data_index
      integer stack, blkndx
      integer allocate_scratch_block
      integer type
      
      integer my_comm_rank, comm
      integer home, iblk
      integer val1(mx_array_index), val2(mx_array_index)
      integer dval1(mx_array_index), dval2(mx_array_index)
      integer index_type(mx_array_index)
      integer dindex_type(mx_array_index)
      double precision x(1)
    
      integer nb, proc, map, company_comm
      integer pst_get_company_comm, pst_get_company_rank
      integer list_unit, data_unit, ios
      integer datatag, handle
      integer icall, nrec

      integer request, drequest
      integer msg(len_sip_server_message)
      integer f_form_msg_tag
      integer darray, next_array
      integer tag, nxt, mgr
      integer status(MPI_STATUS_SIZE)
      integer iarray
      logical is_file_open
      logical dist_flag
      integer*8 data_ptr
      integer block_map_lookup

      double precision esum, block_energy

      company_comm    = pst_get_company_comm(me)
      proc  = pst_get_company_rank(me)
      if (proc .ne. 0) then
         print *,'ERROR: master_block_to_list called on non-master'
         call abort_job()
         return
      endif

c---------------------------------------------------------------------------
c   Allocate a scratch block.
c---------------------------------------------------------------------------
   
      ierr = allocate_scratch_block(x, iscr, handle,
     *                                array_table, narray_table,
     *                              index_table, nindex_table,
     *                              block_map_table)
      if (ierr .ne. 0) then
         print *,'Error: Cannot allocate a scratch block'
         call array_block_summary(array_table, narray_table)
         call dump_block_ids()
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Open the files.
c--------------------------------------------------------------------------

      list_unit = 29
      inquire (unit = list_unit, opened=is_file_open)
      if (is_file_open) then
         rewind list_unit
      else
         open (unit = list_unit, file='BLOCK_INDEX',err = 1000, 
     *      iostat=ios, access='SEQUENTIAL', form='UNFORMATTED',
     *      status='OLD')
      endif

      call f_openfile('BLOCKDATA' // char(0), data_unit)
      data_ptr = 0

c---------------------------------------------------------------------------
c   Find the indices of the array.
c---------------------------------------------------------------------------
       
      dist_flag = .false.

      do iarray = 1, narray_list
      esum = 0.
      array = array_list(iarray)
      type  = array_table(c_array_type,array)
      if (dbg) then
         print *,'Task ',me,' LIST_TO_BLOCKS: Reading array ',array,
     *             ' type ',type
         call prt_time('Worker time')
      endif
      if (type .ne. distributed_array .and.
     *    type .ne. served_array) then
         print *,'Error: LIST_TO_BLOCKS array must be either ',
     *           'a distributed or served array.'
         print *,'Array ',array,' has type = ',type
         call abort_job()
      endif

      nb = array_table(c_numblks, array)
      map = array_table(c_block_map, array)
      nindex = array_table(c_nindex, array)

      if (type .eq. served_array) then

c--------------------------------------------------------------------------
c   Set the "prepared" flag in the array_table, indicating that a PREPARE
c   has been executed on this array.
c--------------------------------------------------------------------------

         array_table(c_prepare_flag, array) = 1
      else

c--------------------------------------------------------------------------
c   Send the array signal to each processor in the company.
c--------------------------------------------------------------------------

         do i = 1, my_company_size
            if (i-1 .ne. my_company_rank) then
               call mpi_send(array, 1, MPI_INTEGER, i-1, 
     *                     block_list_array_tag, company_comm,
     *                     status, ierr)
               if (dbg) print *,'Task ',me,' Sent array message to ',
     *               i-1,' array = ',array
            endif
         enddo
      endif

c---------------------------------------------------------------------------
c   Pick up the array's indices, and save the current segments of each.
c---------------------------------------------------------------------------

      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1,array)
         index_type(i) = index_table(c_index_type, index(i))
      enddo

c-------------------------------------------------------------------------
c   Main processing loop.
c-------------------------------------------------------------------------

      request = MPI_REQUEST_NULL
      drequest = MPI_REQUEST_NULL

      do iblk = 1, nb

c-------------------------------------------------------------------------
c   Read a descriptor record from the BLOCK_INDEX file.
c-------------------------------------------------------------------------

         read (list_unit) darray, dblock, dnindex, n,
     *                       (dindex_type(i),i=1,nindex),
     *                       (dseg(i),i=1,nindex),
     *                       (dval1(i),i=1,nindex),
     *                       (dval2(i),i=1,nindex)

c--------------------------------------------------------------------------
c   Make sure the indices match with respect to number and type.
c--------------------------------------------------------------------------

         ierr = 0
         if (dnindex .ne. nindex) ierr = 1
         do i = 1, nindex
            if (dindex_type(i) .ne. index_type(i)) ierr = 1
         enddo 

c---------------------------------------------------------------------------
c   Make sure the segment values match the segment values on disk.
c---------------------------------------------------------------------------

         do i = 1, nindex
            seg(i) = block_map_table(c_block_map_seg+i-1,map+dblock-1)
            if (seg(i) .ne. dseg(i)) ierr = 1

            call get_index_segment(index(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
            if (dval1(i) .ne. val1(i)) ierr = 1
            if (dval2(i) .ne. val2(i)) ierr = 1
         enddo

         if (ierr .eq. 1) then
            print *,'Error in LIST_TO_BLOCKS: ',
     *               'Index data does not match data ',
     *              'for current job.'
            print *,'nindex ',nindex,' dnindex ',dnindex
            print *,'index_type: ',(index_type(i),i=1,nindex),
     *               ' dindex_type: ',(dindex_type(i),i=1,nindex)
            print *,'seg: ',(seg(i),i=1,nindex),' dseg ',
     *              (dseg(i),i=1,nindex)
            print *,'val1: ',(val1(i),i=1,nindex),' dval1: ',
     *                (dval1(i),i=1,nindex)
            print *,'val2: ',(val2(i),i=1,nindex),' dval2: ',
     *                (dval2(i),i=1,nindex)
            call abort_job()
         endif

c--------------------------------------------------------------------------
c   Wait for a previous descriptor message to complete.
c--------------------------------------------------------------------------

         if (request .ne. MPI_REQUEST_NULL) 
     *      call mpi_wait(request, status, ierr)

         if (type .eq. distributed_array) then
            home = block_map_table(c_processor, map+dblock-1)
            if (home .ne. proc) then

c--------------------------------------------------------------------------
c   Send a descriptor message to the block's owner.
c--------------------------------------------------------------------------

               datatag = f_form_msg_tag() 

               msg(1) = dblock
               msg(2) = datatag
               call mpi_isend(msg, 2, mpi_integer, home,
     *                     block_list_descriptor_tag,
     *                     company_comm, request, ierr)
            endif
         else

c---------------------------------------------------------------------------
c   Post a PREPARE descriptor to the I/O server.
c---------------------------------------------------------------------------

            mgr = block_map_table(c_processor, map+dblock-1)

            msg(1) = sip_server_prepare
            msg(2) = array
            msg(3) = nindex
            tag    = f_form_msg_tag()
            msg(4) = tag
            msg(5) = current_line
            msg(6) = op(c_server_stat_key)
            nxt = 7 
            do i = 1, mx_array_index
               if (i .le. nindex) then
                  msg(nxt) = index(i)
               else
                  msg(nxt) = 0
               endif
               nxt = nxt + 1
            enddo

            n = 1
            do i = 1, nindex
               call get_index_segment(index(i), seg(i), 
     *                  segment_table, nsegment_table, index_table,
     *                  nindex_table, msg(nxt), msg(nxt+1))
               n = n * (msg(nxt+1)-msg(nxt)+1)
               nxt = nxt + 2
            enddo

            do i = nindex+1, mx_array_index
               msg(nxt)   = 0
               msg(nxt+1) = 0
               nxt = nxt + 2
            enddo

            msg(7) = block_map_lookup(seg, nindex, array,
     *                         array_table(1,array),
     *                         index_table, nindex_table)

            call mpi_isend(msg, len_sip_server_message, mpi_integer,
     *                 mgr, sip_server_message,
     *                 mpi_comm_world, request, ierr)
         endif

c---------------------------------------------------------------------------
c   Read the data block into the scratch area.
c---------------------------------------------------------------------------

         if (drequest .ne. MPI_REQUEST_NULL) 
     *       call mpi_wait(drequest, status, ierr)

c--------------------------------------------------------------------------
c   Send the data in scratch to its owner.
c--------------------------------------------------------------------------

         if (type .eq. distributed_array) then
            dist_flag = .true.
            if (home .eq. proc) then

c--------------------------------------------------------------------------
c   The block is owned by the master.  Simply read it into its block,
c   which should already exist due to a previous CREATE instruction.
c--------------------------------------------------------------------------

               stack = array_table(c_array_stack,array)
               blkndx = block_map_table(c_bmap_blkndx,map+dblock-1)
               indblk = get_block_data_index(array, dblock, stack, 
     *                                       blkndx, x)
c               read (data_unit) (x(indblk+i-1), i=1,n) 
               call f_read_disk(data_unit, data_ptr, x(indblk), n)
               data_ptr = data_ptr + n 

c               call check_orbital_ind(x(indblk),index_type, val1, 
c     *                           val2, nindex) 
               if (dbg) esum = esum + block_energy(x(indblk), n)
            else

c--------------------------------------------------------------------------
c   Disributed array, not owned by master.
c--------------------------------------------------------------------------

c               read (data_unit) (x(iscr+i-1),i=1,n)
               call f_read_disk(data_unit, data_ptr, x(iscr), n)
               data_ptr = data_ptr + n

c               call check_orbital_ind(x(iscr),index_type, val1,
c     *                           val2, nindex)
               if (dbg) esum = esum + block_energy(x(iscr), n)
               call mpi_isend(x(iscr), n, mpi_double_precision, home,
     *                     datatag,
     *                     company_comm, drequest, ierr)
            endif   
         else

c--------------------------------------------------------------------------
c   Served array.
c--------------------------------------------------------------------------

            call f_read_disk(data_unit, data_ptr, x(iscr), n)
            data_ptr = data_ptr + n 
c            call check_orbital_ind(x(iscr),index_type, val1,
c     *                           val2, nindex)
            if (dbg) esum = esum + block_energy(x(iscr), n)
            call mpi_isend(x(iscr), n, mpi_double_precision,
     *                     mgr, tag,
     *                     mpi_comm_world, drequest, ierr)
         endif
      enddo

c--------------------------------------------------------------------------
c   Do the final waits if necessary.
c--------------------------------------------------------------------------

      call mpi_wait(request, status, ierr)
      call mpi_wait(drequest, status, ierr)

      if (dbg) print *,'LIST_TO_BLOCKS array ',array,' at line ',
     *    current_line,' esum ',esum

c---------------------------------------------------------------------------
c   Send a descriptor message with a -1 for the tag as a signal that we are
c   through with the current array.
c---------------------------------------------------------------------------

         if (type .ne. served_array) then
            msg(1) = -1
            msg(2) = -1
            do i = 1, my_company_size
               if (i-1 .ne. me) then
                  call mpi_send(msg, 2, MPI_INTEGER, i-1,
     *                    block_list_descriptor_tag,
     *                    company_comm, status,ierr)
               endif
            enddo
         endif

      enddo   ! iarray

c--------------------------------------------------------------------------
c   Send an array signal with a -1 to each non-server process.
c   This functions as a quit signal.
c--------------------------------------------------------------------------

      msg(1) = -1
      do i = 1, my_company_size
         if (i-1 .ne. me) then
            call mpi_send(msg, 1, MPI_INTEGER, i-1,
     *                    block_list_array_tag,
     *                    company_comm, status,ierr)
            if (dbg) print *,'Task ',me,' Sent quit signal to ',
     *                         i-1
         endif
      enddo

c---------------------------------------------------------------------------
c   Close the files, release the scratch block.
c---------------------------------------------------------------------------

      close(list_unit)
      call f_close_file(data_unit)
      call free_scratch_block(handle) 
      narray_list = 0
      call prt_time('END OF LIST_TO_BLOCKS')
      
      return

 1000 continue
      print *,'Cannot open BLOCK_INDEX or BLOCKDATA files.'
      print *,'I/O status = ',ios
      call abort_job()
      end

      subroutine check_orbital_ind(x,index_type, val1, val2, nindex)
      implicit none
      include 'interpreter.h'
      include 'dropmo.h'
      integer*8 x(*)

      integer nindex
      integer index_type(nindex), val1(nindex), val2(nindex)

      integer i, j, k, next, nzero
      integer orb(4), match(4)
      integer*8 data, save_data
      integer*8 unpackit(4), shift1, shift2, shift3

      logical mismatch, printit

      next = 0

      do i = 1, nindex
         orb(i) = val1(i)
      enddo
 
      shift1 = 2**16

  100 continue
      next = next + 1
      data = x(next)
      save_data = data

      do i = 1, 4
         unpackit(i) = -1
      enddo

      mismatch = .false.
      nzero = 0
      do i = 1, 4
         match(i) = orb(i)
         if (index_type(i) .eq. moindex .or. 
     *       index_type(i) .eq. moaindex) then
            do k = 1, ndropmo_a
               if (orb(i) .eq. modrop_a(k)) then
                  match(i) = 0
                  nzero = nzero + 1
               endif
            enddo
         else if (index_type(i) .eq. mobindex) then
            do k = 1, ndropmo_b
               if (orb(i) .eq. modrop_b(k)) then
                  match(i) = 0
                  nzero = nzero + 1
               endif
            enddo
         endif
      enddo

c      print 101,(orb(k),k=1,4),(match(k),k=1,4),data,nzero
c  101 format(' orb ',4(i4,1x),' match ',4(i4,1x),z16,' nzero ',i4)
 
      if (nzero .eq. 0) then
         unpackit(1) = data - data / shift1 * shift1
         data        = data / shift1
         unpackit(2) = data - data / shift1 * shift1
         data        = data / shift1
         unpackit(3) = data - data / shift1 * shift1
         unpackit(4) = data / shift1
      
         do i = 1, 4
            if (unpackit(i) .ne. match(i) ) mismatch = .true.
         enddo 
      else
         if (save_data .ne. 0) mismatch = .true.
      endif

      if (mismatch) then
         print *,'Error: data mismatch at sample ',(orb(k),
     *           k=1,4)
         print 102,(unpackit(k),k=1,4),(match(k),k=1,4),save_data
  102 format(' data is ',4(i4,1x),' should be ',4(i4,1x),' raw data ',
     *        z16)
      endif

      j = 1
  200 continue
      orb(j) = orb(j) + 1
      if (orb(j) .gt. val2(j)) then
         orb(j) = val1(j)
         j = j + 1
         if (j .gt. nindex) return
         go to 200
      else
        go to 100
      endif

      return
      end
