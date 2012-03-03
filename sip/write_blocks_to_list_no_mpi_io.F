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
      subroutine write_blocks_to_list_no_mpi_io(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, 
     *                      address_table, op)
c--------------------------------------------------------------------------
c   Dumps all blocks of an array to a list.  This version is used when 
c   a system does not support MPI_IO.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
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
      integer saveseg(mx_array_index)
      integer find_current_block
      integer blkndx
      integer*8 indblk, get_block_data_index
      integer stack
      integer f_form_msg_tag
      integer type, next, nreq
      
      integer mx_blk
      parameter (mx_blk = 200)

      integer request(mx_blk)
      integer msg(2,mx_blk)
      integer statuses(MPI_STATUS_SIZE,mx_blk)

      integer my_comm_rank, comm
      integer home, iblk
      integer val1(mx_array_index), val2(mx_array_index)
      double precision x(1)
    
      integer nb, proc, map, company_comm
      integer tag
      integer pst_get_company_comm, pst_get_company_rank
      integer iarray

      company_comm    = pst_get_company_comm(me)
      proc  = pst_get_company_rank(me)
      if (proc .eq. 0) then
         call master_blocks_to_list(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, op)
         return
      endif

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      do iarray = 1, narray_list
      array = array_list(iarray)
      type  = array_table(c_array_type,array)
      if (type .eq. served_array) go to 900 ! master handles all served arrays
      if (type .ne. distributed_array) then
         print *,'Error: BLOCKS_TO_LIST array must be either ',
     *     'a distributed or served array.'
         print *,'Array ',array,' has type = ',type
         call abort_job()
      endif

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

c-------------------------------------------------------------------------
c   Loop over all blocks in the block_map_table, processing only the
c   blocks with this processor as "home".
c-------------------------------------------------------------------------

      nreq = 0
      next = 0
      do iblk = 1, nb
         home = block_map_table(c_processor, map+iblk-1)
         if (home .eq. proc) then 

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

            n = 1
            do i = 1, nindex
               seg(i) = block_map_table(c_block_map_seg+i-1,map+iblk-1)
               call get_index_segment(index(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))

c--------------------------------------------------------------------------
c   Set the current_segment of the appropriate index to this segment.
c--------------------------------------------------------------------------

                index_table(c_current_seg,index(i)) = seg(i)
                n = n * (val2(i)-val1(i)+1)
            enddo

c---------------------------------------------------------------------------
c   Get array data address.
c---------------------------------------------------------------------------

            block = find_current_block(array, array_table(1,array),
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)
            if (block .lt. 0) then
               print *,'Error in BLOCKS_TO_LIST: Array ',array,
     *            ' Block ',block,' not found.'
               call abort_job()
            endif

            stack  = array_table(c_array_stack,array)
            indblk = get_block_data_index(array, block, stack,
     *                                    blkndx, x)

c--------------------------------------------------------------------------
c   Send a descriptor message to the master.
c--------------------------------------------------------------------------
   
            tag = f_form_msg_tag() 
            next = next + 1
            if (nreq .gt. mx_blk) then
               call mpi_waitall(nreq, request, statuses, ierr)
               next = 1
            endif

            msg(1,next) = iblk
            msg(2,next) = tag
            nreq = nreq + 1
            call mpi_isend(msg(1,next), 2, MPI_INTEGER, 0, 
     *                  block_list_descriptor_tag,
     *                  company_comm, request(nreq), ierr)

c---------------------------------------------------------------------------
c   Send the data to the master process.
c---------------------------------------------------------------------------

            nreq = nreq + 1
            call mpi_isend(x(indblk), n, MPI_DOUBLE_PRECISION, 0,
     *                  tag, company_comm, request(nreq), ierr)
         endif 
      enddo

c---------------------------------------------------------------------------
c   Restore the current segments to their saved values.
c---------------------------------------------------------------------------

      do i = 1, nindex
         index_table(c_current_seg,index(i)) = saveseg(i)
      enddo

c----------------------------------------------------------------------------
c   Wait for all messages to complete.
c----------------------------------------------------------------------------

      call mpi_waitall(nreq, request, statuses, ierr)
  900 continue
      enddo   ! iarray

      narray_list = 0
      return
      end

      subroutine master_blocks_to_list(array_table, narray_table, 
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, op)
c--------------------------------------------------------------------------
c   Dumps all blocks of an array to a list.
c   This routine handles the master's role.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blockdata_list.h'
      include 'machine_types.h'
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

      integer i, j, k, n
      integer array, index(mx_array_index), nindex, ierr
      integer seg(mx_array_index)
      integer*8 iscr
      integer*8 indblk
      integer*8 get_block_data_index
      integer allocate_scratch_block
      integer type
      
      integer my_comm_rank, comm
      integer home, iblk
      integer val1(mx_array_index), val2(mx_array_index)
      integer index_type(mx_array_index)
      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif
    
      integer nblk, proc, map, company_comm
      integer pst_get_company_comm, pst_get_company_rank
      integer list_unit, data_unit, ios
      logical file_open
      integer irecv, nrecv
      integer datatag, nind, handle

      integer mx_blk
      parameter (mx_blk = 25000)
      integer request(mx_blk)
      integer tag(mx_blk)
      integer desc_msg(2,mx_blk)
      integer stack, blkndx
      integer status(MPI_STATUS_SIZE)
      integer iarray
      integer*8 data_ptr

      logical served_flag
      double precision esum, block_energy

#ifdef ALTIX
      dptr = dshptr
#endif

      company_comm    = pst_get_company_comm(me)
      proc  = pst_get_company_rank(me)
      if (proc .ne. 0) then
         print *,'ERROR: master_block_to_list called on non-master'
         call abort_job()
         return
      endif

c---------------------------------------------------------------------------
c   Verify array type info.
c---------------------------------------------------------------------------

      do iarray = 1, narray_list
         array = array_list(iarray)
         type  = array_table(c_array_type,array)
         if (type .ne. distributed_array .and.
     *       type .ne. served_array) then
            print *,'Error: BLOCKS_TO_LIST array must be ',
     *              'either a distributed or served array.'
            print *,'Array ',array,' has type = ',type
            call abort_job()
         endif
      enddo

c---------------------------------------------------------------------------
c   Allocate a scratch block.
c---------------------------------------------------------------------------
   
      ierr = allocate_scratch_block(x, iscr, handle, array_table,
     *                              narray_table, index_table, 
     *                              nindex_table, block_map_table)
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
      inquire (unit=list_unit, opened=file_open)
      if (.not. file_open) then
         open (unit = list_unit, file='BLOCK_INDEX2',err = 1000,
     *      iostat=ios, access='SEQUENTIAL', form='UNFORMATTED')
      endif

      data_unit = -1
      inquire (file='BLOCKDATA2', exist = file_open)
      if (file_open) call f_deletefile('BLOCKDATA2'//char(0), 
     *                                 data_unit)
      call f_creatfile('BLOCKDATA2' // char(0), data_unit)
      data_ptr = 0
 
c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      served_flag = .false.

      do iarray = 1, narray_list
      esum = 0.

      array = array_list(iarray)
      type  = array_table(c_array_type,array)

      if (type .eq. served_array) then
         call master_blocks_to_list_served(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, op, x(iscr),
     *                      list_unit, data_unit, array, data_ptr,
     *                      esum)
         go to 900    ! process next array
      endif

      nblk = array_table(c_numblks, array)
      map = array_table(c_block_map, array)
      nindex = array_table(c_nindex, array)

      if (dbg) print *,'Task ',me,' BLOCKS_TO_LIST: array ',
     *     array,' type ',
     *     type,' nb ',nblk
      call prt_time('Worker time')
      if (nblk .gt. mx_blk) then
         print *,'Error: Array ',array,' contains ',nblk,' blocks'
         print *,'Current max. for lists is ',mx_blk
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Pick up the array's indices, and save the current segments of each.
c---------------------------------------------------------------------------

      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1,array)
         index_type(i) = index_table(c_index_type, index(i))
      enddo

c-------------------------------------------------------------------------
c   Loop over all blocks in the block_map_table, processing only the
c   blocks with this processor as "home".
c-------------------------------------------------------------------------

      nrecv = 0
      do iblk = 1, nblk
         home = block_map_table(c_processor, map+iblk-1)

c--------------------------------------------------------------------------
c   Post a descriptor message recv from the block's owner.
c--------------------------------------------------------------------------

         if (home .ne. proc) then
            nrecv = nrecv + 1
            call mpi_irecv(desc_msg(1,nrecv), 2, mpi_integer, home,
     *                     block_list_descriptor_tag,
     *                     company_comm, request(nrecv), ierr)
         endif
      enddo

c---------------------------------------------------------------------------
c   Process the data messages.
c---------------------------------------------------------------------------

      do irecv = 1, nrecv
         call mpi_wait(request(irecv), status, ierr)

         iblk    = desc_msg(1,irecv)
         datatag = desc_msg(2,irecv)
         home    = block_map_table(c_processor, map+iblk-1)

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

         n = 1
         do i = 1, nindex
            seg(i) = block_map_table(c_block_map_seg+i-1,map+iblk-1)
            call get_index_segment(index(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
            n = n * (val2(i)-val1(i)+1)
         enddo

c---------------------------------------------------------------------------
c   Write an index record for this block.
c---------------------------------------------------------------------------

         write (list_unit) array, iblk, nindex, n,
     *                       (index_type(i),i=1,nindex),
     *                       (seg(i),i=1,nindex),
     *                       (val1(i),i=1,nindex), 
     *                       (val2(i),i=1,nindex)

c----------------------------------------------------------------------------
c   Receive the data into the scratch block.
c----------------------------------------------------------------------------

         call mpi_recv(x(iscr), n, mpi_double_precision, home,
     *                     datatag,
     *                     company_comm, status, ierr)

c----------------------------------------------------------------------------
c   Write the data record for this block.
c----------------------------------------------------------------------------

c         write (data_unit) (x(iscr+i-1),i=1,n)
c          call encode_orbital_ind(x(iscr), index_type, val1, val2, 
c     *                            nindex)
          call f_write_disk(data_unit, data_ptr, x(iscr), n)
          data_ptr = data_ptr + n
          if (dbg) esum = esum + block_energy(x(iscr), n)
      enddo

c--------------------------------------------------------------------------
c   Write all the master's blocks.
c--------------------------------------------------------------------------

      do iblk = 1, nblk
         home = block_map_table(c_processor, map+iblk-1)
         if (home .eq. proc) then

c-------------------------------------------------------------------------
c   Get segment ranges.
c-------------------------------------------------------------------------

            n = 1
            do i = 1, nindex
               seg(i) = block_map_table(c_block_map_seg+i-1,map+iblk-1)
               call get_index_segment(index(i), seg(i), segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, val1(i), val2(i))
               n = n * (val2(i)-val1(i)+1)
            enddo

c---------------------------------------------------------------------------
c   Write an index record to the index file for this block.
c---------------------------------------------------------------------------

            write (list_unit) array, iblk, nindex, n,
     *                       (index_type(i),i=1,nindex),
     *                       (seg(i),i=1,nindex),
     *                       (val1(i),i=1,nindex), 
     *                       (val2(i),i=1,nindex)

c--------------------------------------------------------------------------
c   Write the data block.
c--------------------------------------------------------------------------

            stack = array_table(c_array_stack,array)
            blkndx = block_map_table(c_bmap_blkndx,map+iblk-1)
            indblk = get_block_data_index(array, iblk, 
     *                                    stack, blkndx, x)
c            write (data_unit) (x(indblk+i-1),i=1,n)
c            call encode_orbital_ind(x(indblk), index_type, val1,  
c     *                              val2, nindex)
            call f_write_disk(data_unit, data_ptr, x(indblk), n) 
            data_ptr = data_ptr + n
            if (dbg) esum = esum + block_energy(x(indblk), n)
         endif
      enddo

  900 continue
         if (dbg) then
            print *,'BLOCKS_TO_LIST: array ',array,' at line ',
     *          current_line,' esum = ',esum
         endif
      enddo   ! iarray

c---------------------------------------------------------------------------
c   Release the scratch block.
c---------------------------------------------------------------------------

      call free_scratch_block(handle) 
      narray_list = 0

c---------------------------------------------------------------------------
c   Close the files.
c---------------------------------------------------------------------------

      close(list_unit)
      call f_close_file(data_unit)

c---------------------------------------------------------------------------
c   Rename BLOCK_INDEX2 --> BLOCK_INDEX, BLOCKDATA2 --> BLOCKDATA.
c---------------------------------------------------------------------------

      call f_renamefile('BLOCK_INDEX2'//char(0),'BLOCK_INDEX'//char(0))
      call f_renamefile('BLOCKDATA2'//char(0),'BLOCKDATA'//char(0))
      call prt_time('END OF BLOCKS_TO_LIST')
      return

 1000 continue
      print *,'Cannot open BLOCK_INDEX or BLOCKDATA files.'
      print *,'I/O status = ',ios
      call abort_job()
      end

      subroutine master_blocks_to_list_served(array_table, 
     *                      narray_table, index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, op, scr,
     *                      list_unit, data_unit, array, data_ptr,
     *                      esum)
c--------------------------------------------------------------------------
c   Dumps all blocks of an array to a list.
c   This routine handles the master's role.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'int_gen_parms.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
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
      double precision scr(*)

      integer i, j, k, n
      integer array, index(mx_array_index), nindex, ierr
      integer msg_array, msg_nindex
      integer msg_tag
      integer seg(mx_array_index)
      integer type
      
      integer my_comm_rank
      integer msg_iblk
      integer val1(mx_array_index), val2(mx_array_index)
      integer msg_val1(mx_array_index), msg_val2(mx_array_index)
      integer index_type(mx_array_index)
    
      integer proc, map
      integer list_unit, data_unit
      integer request(nprocs)
      integer server_request(nprocs)
      integer server_id(nprocs)
      integer nxt, match
      integer msg(len_sip_server_message),
     *        desc_msg(len_sip_server_message,nprocs)
      integer f_form_msg_tag
      integer status(MPI_STATUS_SIZE)
      integer company, pst_get_company

      integer*8 data_ptr
      integer nblocks, nwpblk, niocompany
      integer nb_processed, irecv, source
      double precision esum, block_energy
      integer statuses(MPI_STATUS_SIZE, nprocs)
      logical invalid_data

c---------------------------------------------------------------------------
c   Find the indices of the array block.
c---------------------------------------------------------------------------
       
      type  = array_table(c_array_type,array)
      if (type .ne. served_array) then
         print *,'Error: MASTER_BLOCKS_TO_LIST_SERVED must be ',
     *           'called with a served array.'
         print *,'Array ',array,' has type = ',type
         call abort_job()
      endif

      nb = array_table(c_numblks, array)
      map = array_table(c_block_map, array)
      nindex = array_table(c_nindex, array)
      if (dbg) then
        print *,'BLOCKS_TO_LIST array ',array,' type ',
     *     type,' nb ',
     *     nb
         call prt_time('Worker time')
      endif

c---------------------------------------------------------------------------
c   Send a "sip_server_blocks_to_list" message to each server, requesting
c   that the servers transfer the blocks of the array.
c---------------------------------------------------------------------------

      msg(1) = sip_server_blocks_to_list
      msg(2) = array 
      msg(3) = 0         ! nindex
      niocompany = 0
      do i = 1, nprocs
         company = pst_get_company(i-1) 
         if (company .eq. io_company_id) then
            niocompany = niocompany + 1
            call mpi_isend(msg, len_sip_server_message, mpi_integer,
     *                  i-1, sip_server_message,
     *                  mpi_comm_world, 
     *                  server_request(niocompany), ierr) 
            server_id(niocompany) = i-1
         endif 
      enddo

      call mpi_waitall(niocompany, server_request, statuses, ierr)

c---------------------------------------------------------------------------
c   Post an initial recv for each server.
c---------------------------------------------------------------------------

      do i = 1, niocompany
         call mpi_irecv(desc_msg(1,i), len_sip_server_message, 
     *                     mpi_integer, server_id(i),
     *                     sip_server_blocks_to_list,
     *                     mpi_comm_world, request(i), ierr)
      enddo 
      
c---------------------------------------------------------------------------
c   Pick up the array's indices and index_type's.
c---------------------------------------------------------------------------

      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1,array)
         index_type(i) = index_table(c_index_type, index(i))
      enddo

c----------------------------------------------------------------------------
c   Process the blocks of this array.
c----------------------------------------------------------------------------

      nb_processed = 0
  100 continue
      
c----------------------------------------------------------------------------
c   Wait for a recv of the desc_msg.
c----------------------------------------------------------------------------

      call mpi_waitany(niocompany, request, irecv, status, ierr)

c----------------------------------------------------------------------------
c   Decode the descriptor message and write the data to BLOCK_INDEX.
c----------------------------------------------------------------------------

         msg_array = desc_msg(2,irecv)
         msg_nindex = desc_msg(3,irecv) 
         msg_tag    = desc_msg(4,irecv)
         msg_iblk   = desc_msg(7,irecv) 

         nxt = 7 + mx_array_index

         n = 1
         do i = 1, mx_array_index
            if (i .le. msg_nindex) then
               msg_val1(i) = desc_msg(nxt,irecv)
               msg_val2(i) = desc_msg(nxt+1,irecv)
               n = n * (msg_val2(i) - msg_val1(i) + 1)
            endif
            nxt = nxt + 2
         enddo

c---------------------------------------------------------------------------
c   Validate the data in the message.
c---------------------------------------------------------------------------

         invalid_data = .false. 
         if (array .ne. msg_array .or. nindex .ne. msg_nindex) then
            print *,'Error in write_blocks_to_list: Invalid data in msg'
            print *,'array ',array,' msg_array ',msg_array
            print *,'nindex ',nindex,' msg_nindex ',msg_nindex
            do i = 1, nindex
               print *,'val1, val2 ',val1(i),val2(i),
     *                 ' msg_val1, msg_val2 ',msg_val1(i),msg_val2(i)
            enddo
            call abort_job()
         endif
     
         do j = 1, nindex
            seg(j) = 
     *         block_map_table(c_block_map_seg+j-1,msg_iblk)
            call get_index_segment(index(j), seg(j), segment_table,
     *                          nsegment_table, index_table,
     *                          nindex_table, val1(j), val2(j))
         enddo

         do j = 1, nindex
            if (val1(j) .ne. msg_val1(j) .or. 
     *          val2(j) .ne. msg_val2(j)) then
               print *,'Error: seg values do not match.'
               print *,'index ',(index(i),i=1,nindex) 
               print *,'seg ',(seg(i),i=1,nindex) 
               print *,'val1 ',(val1(i),i=1,nindex) 
               print *,'msg_val1 ',(msg_val1(i),i=1,nindex) 
               print *,'val2 ',(val2(i),i=1,nindex) 
               print *,'msg_val2 ',(msg_val2(i),i=1,nindex) 
            endif
         enddo

c--------------------------------------------------------------------------
c   Look up the block in the blk_map_table, verify the blocksize.
c--------------------------------------------------------------------------

         source = block_map_table(c_processor, msg_iblk)
         if (source .ne. server_id(irecv)) then
            print *,'Error: write_blocks_to_list'
            print *,'Block ',msg_iblk,' recvd from wrong source.'
            print *,'Msg source ',server_id(irecv),' should be ',
     *                source
            call abort_job()
         endif
 
c---------------------------------------------------------------------------
c   Write an index record to the index file for this block.
c---------------------------------------------------------------------------

         write (list_unit) array, msg_iblk-map+1, nindex, n,
     *                       (index_type(i),i=1,nindex),
     *                       (seg(i),i=1,nindex),
     *                       (val1(i),i=1,nindex),
     *                       (val2(i),i=1,nindex)

c--------------------------------------------------------------------------
c   Now recv the data block, write it to disk.
c--------------------------------------------------------------------------

      call mpi_recv(scr, n, mpi_double_precision,
     *               server_id(irecv), msg_tag, mpi_comm_world,
     *               status, ierr)

      call f_write_disk(data_unit, data_ptr, scr, n)
      data_ptr = data_ptr + n
      if (dbg) esum = esum + block_energy(scr, n)

c---------------------------------------------------------------------------
c   Post a new recv for this server.
c---------------------------------------------------------------------------

      call mpi_irecv(desc_msg(1,irecv), len_sip_server_message,
     *                     mpi_integer, server_id(irecv),
     *                     sip_server_blocks_to_list,
     *                     mpi_comm_world, request(irecv), ierr)
      nb_processed = nb_processed + 1
      if (nb_processed .lt. nb) go to 100

c--------------------------------------------------------------------------
c   Free up the posted recvs of the desc_msg's.
c--------------------------------------------------------------------------

      do i = 1, niocompany
         if ( request(i) .ne. MPI_REQUEST_NULL) then
            call mpi_cancel(request(i), ierr)
            call mpi_wait(request(i), status, ierr)
         endif
      enddo

      return
      end

      subroutine encode_orbital_ind(x, index_type, val1, val2, nindex)
      implicit none
      include 'interpreter.h'
      include 'dropmo.h'

      integer nindex
      integer index_type(nindex), val1(nindex), val2(nindex)
      integer*8 x(*)

      integer i, j, k
      integer*8 orb(4), neworb(4), next
      integer*8 packit, shift1, shift2, shift3

      if (nindex .ne. 4) return
      do i = 1, 4
         orb(i) = val1(i)
      enddo

      shift1 = 2**16
      shift2 = shift1 * shift1
      shift3 = shift2 * shift1

      next = 0
  100 continue
      do i = 1, 4
         neworb(i) = orb(i)
         if (index_type(i) .eq. moaindex .or.
     *       index_type(i) .eq. moindex) then
            neworb(i) = moamap(orb(i))
         else if (index_type(i) .eq. mobindex) then
            neworb(i) = mobmap(orb(i))
         endif
      enddo

      packit = neworb(1)
      packit = packit + neworb(2) * shift1
      packit = packit + neworb(3) * shift2 
      packit = packit + neworb(4) * shift3
      next = next + 1
      x(next) = packit
c      print 101,(orb(i),i=1,4),x(next),packit 
c  101 format(' orb ',4(i4,1x),' x(next) ',z16,' packit ',z16)

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
