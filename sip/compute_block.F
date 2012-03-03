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
      subroutine compute_block(op, array_table, narray_table,
     *                         index_table, nindex_table, 
     *                         block_map_table, nblock_map_table,
     *                         segment_table, nsegment_table,
     *                         scalar_table, nscalar_table,
     *                         address_table, 
     *                         debugit, validate, 
     *                         flopcount, comm, comm_timer, 
     *                         instruction_timer)
c--------------------------------------------------------------------------
c   Compute the block corresponding to the result array of the "op" 
c   entry. 
c
c   op             	Operation table entry.
c   array_table		Table of array information.
c   narray_table	Number of arrays in table.
c   debug		Logical flag for debug print
c   validate            Logical flag for block validation.
c   flopcount           Number of floating point operations performed.
c   send_count          Returned number of blocks sent (0 or 1).
c   comm                MPI Communicator
c   comm_timer          Timer handle to accumulate wait time on blocks.
c   instruction_timer   Timer handle for current instruction timer.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'blkmgr.h'
      include 'machine_types.h'
      include 'dbugcom.h'
      include 'checkpoint_data.h'

      integer narray_table
      integer nindex_table, index_table(lindex_table_entry,nindex_table)
      integer op(*), array_table(larray_table_entry,narray_table)
      integer nblock_map_table
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer nsegment_table
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer nscalar_table
      double precision scalar_table(nscalar_table)
      integer*8 address_table(narray_table)
      integer comm
      integer comm_timer, instruction_timer
      logical debugit, debug, validate, create_block
      logical partial_create

      integer find_current_block
      integer find_current_block_map
      integer i, j, k, opcode, block, blkndx, result
      integer op1_blkndx, op2_blkndx
      integer nindex
      integer op1, op2
      integer op1_block, op2_block
      integer ncind, cind(mx_array_index)
      integer flag, flag1, flag2
      integer result_type, op1_type, op2_type
      integer flopcount
      integer ind(mx_array_index), segsave(mx_array_index), 
     *        nb, iblock, l,ind_input(mx_array_index)
      integer iblock_map
      integer get_block_number
      integer request
      integer get_block_request
      integer status(MPI_STATUS_SIZE)
      integer ierr, allocate_block
      integer result_flag
      integer size
      integer dummy
      integer stack
      
      integer push_do, pop_do
      logical direct_flag
      logical prefetch_flag

      opcode = op(c_opcode)
      if (opcode .eq. do_op .or.
     *    opcode .eq. enddo_op .or.
     *    opcode .eq. pardo_op .or.
     *    opcode .eq. endpardo_op .or.
     *    opcode .eq. go_to_op .or.
     *    opcode .eq. jz_op .or.
     *    opcode .eq. exit_op .or.
     *    opcode .eq. cycle_op .or.
     *    opcode .eq. where_op) return

      result = op(c_result_array)
      result_type = array_table(c_array_type, result)
      if (result .lt. 0) then
         print *,'Error: Invalid array handle in operation'
         print *,'Array = ',result
         print *,'OP: ',(op(i),i=1,loptable_entry)
         call abort_job()
      endif

      create_block = .false.
      debug        = debugit
      flopcount    = 0
      
      if (opcode .eq. reindex_op) then

c---------------------------------------------------------------------------
c   Reindex operation. 
c---------------------------------------------------------------------------

         nindex = array_table(c_nindex,result)
         do i = 1, nindex
            cind(i) = op(c_ind1+i-1)
            if (cind(i) .le. 0 .or.
     *          cind(i) .gt. nindex_table) then
               print *,'Task ',me,' Error: Invalid index in operation'
               print *,'Operation: ',(op(j),j=1,loptable_entry)
               call abort_job()
            endif
         enddo

         if (nindex .lt. mx_array_index) then
            do i = nindex+1,mx_array_index
               cind(i) = 0
            enddo
         endif
 
         call set_effective_indices(array_table(1,result), cind)
         return
      endif

c---------------------------------------------------------------------------
c   Index arithmetic and logic operations
c---------------------------------------------------------------------------

      if ((opcode .ge. sp_add_op .and. 
     *    opcode .le. sp_ldindex_op) .or.
     *    opcode .eq. sp_ldi_sym_op) then
         call handle_sp_op(op, index_table, nindex_table)
         return
      endif

c--------------------------------------------------------------------------
c   Scalar logic and arithmetic instructions.
c--------------------------------------------------------------------------

      if (opcode .ge. fl_add_op .and.
     *    opcode .le. fl_load_value_op) then
         call handle_fsp_op(op, array_table, narray_table,
     *                      scalar_table, nscalar_table)
         return
      endif

c---------------------------------------------------------------------------
c   If a distributed array is being used in any instruction except "create"
c   (or reindex_op), check to make sure it has already been created.
c---------------------------------------------------------------------------

      call check_create_flag(op, array_table, narray_table)

c--------------------------------------------------------------------------
c   See if the instruction's result array block is available.
c--------------------------------------------------------------------------

      direct_flag = .false.
      if (result .gt. 0 .and. result_type .ne. static_array) then
         block = find_current_block(result, array_table(1,result),
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           block_map_table, blkndx)
         if (block .gt. 0) then

c--------------------------------------------------------------------------
c   Check the block's communication request to make sure it is
c   not engaged in communication.  If it is, we must wait until
c   the communication is complete before proceeding with the calculation.
c
c   The only exceptions are the GET and REQUEST instructions (get_op, 
c   request_op), which are checked internally to their handler subroutines.
c--------------------------------------------------------------------------

            if (opcode .ne. get_op .and.
     *          opcode .ne. request_op) then
               request = get_block_request(result, block, blkndx)
               if (request .ne. mpi_request_null) 
     *             call wait_on_block(result, block, blkndx,
     *                result_type, request,
     *                instruction_timer, comm_timer)
            endif   
         else
            direct_flag = .true.   ! block does not exist, no need to attempt
                                   ! further lookups.
         endif
      endif

      if (opcode .eq. user_sub_op) then

c---------------------------------------------------------------------------
c   Call user-written code with standard args.
c---------------------------------------------------------------------------

         call exec_user_sub(op(c_user_sub), array_table, narray_table,
     *               index_table, nindex_table, segment_table, 
     *               nsegment_table, block_map_table, nblock_map_table,
     *               scalar_table, nscalar_table, address_table,
     *               op)
         return
      endif

      if (opcode .eq. put_op .or.
     *    opcode .eq. put_replace_op) then

c------------------------------------------------------------------------
c   Accumulate operation.
c------------------------------------------------------------------------

         call global_accumulate(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table, address_table,
     *                      comm, op, comm_timer, instruction_timer)
         return
      endif

      if (opcode .eq. collective_sum_op) then

c------------------------------------------------------------------------
c   Collective sum operation (i. e., scalar accumulate on all procs).
c------------------------------------------------------------------------

         call collective_sum(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      scalar_table, nscalar_table,
     *                      comm, op)
         return
      endif

c---------------------------------------------------------------------------
c   Create operation.   Create the blocks of a distributed array.
c---------------------------------------------------------------------------

      if (opcode .eq. create_op) then
         if (result_type .ne. distributed_array) then
            print *,'Error: Create operation called for array ',
     *               result,'.  Array is not distributed.'
            print *,'Array type is ',result_type
            call abort_job()
         endif
         
         if (array_table(c_create_flag, result) .ne. 0) then
            print *,'Task ',me,' Error: Attempt to create an array ',
     *              'that has already been created.'
            print *,'Array = ',result,' Array_table entry = ',
     *              (array_table(j,result),j=1,larray_table_entry)
            call abort_job()
         endif

c----------------------------------------------------------------------------
c   Save the segment states of this array.
c----------------------------------------------------------------------------

         nindex = array_table(c_nindex,result)
         do k = 1, nindex
            ind(k) = array_table(c_index_original+k-1,result)
            segsave(k) = index_table(c_current_seg,ind(k))
            ind_input(k) = array_table(c_index_array1+k-1,result)
            array_table(c_index_array1+k-1,result) =
     *          array_table(c_index_original+k-1,result)
         enddo

c---------------------------------------------------------------------------
c   Warn user if the array was declared with repeated indices.
c---------------------------------------------------------------------------

         do k = 1, nindex
            do l = 1, nindex
               if (k .ne. l .and.
     *             ind(k) .eq. ind(l)) then
                  print *,'*** WARNING *** Array number ',result,
     *               ' was declared with repeated indices.  This could',
     *               ' result in the creation of only the diagonal ',
     *               'blocks of the array, causing surprising results.'
                  go to 100
               endif
            enddo
         enddo

c--------------------------------------------------------------------------
c   Partial create, or full create?
c---------------------------------------------------------------------------

         partial_create = .false.
         do k = 1, nindex
            if (op(c_ind1+k-1) .eq. wildcard_indicator) 
     *            partial_create = .true.
         enddo
         if (dbg .and. partial_create) 
     *      print *,'PARTIAL CREATE at line ',
     *         current_line,' op ',(op(l),l=1,loptable_entry)
 
  100    continue
         nb = array_table(c_numblks, result)
         iblock_map = array_table(c_block_map, result)

         do 200 k = 1, nb
            if (block_map_table(c_processor,iblock_map+k-1) .eq. 
     *          my_company_rank) then

               if (partial_create) then
                  do l = 1, nindex
                     if (op(c_ind1+l-1) .ne. wildcard_indicator) then

c-----------------------------------------------------------------------------
c   Check block_map table's index values against the one's saved in ind.
c   If all such indices match, we fall through, and the block will be 
c   created.
c-----------------------------------------------------------------------------

                        if (index_table(c_current_seg,ind(l)) .ne.
     *                      block_map_table(c_block_map_seg+l-1,
     *                               iblock_map+k-1)) go to 200
                     endif
                  enddo
               endif
 
c----------------------------------------------------------------------------
c   Set up the "current segments" of this block and calculate the blocksize.
c----------------------------------------------------------------------------

               do l = 1, nindex
                  index_table(c_current_seg,ind(l)) =
     *               block_map_table(c_block_map_seg+l-1,
     *                               iblock_map+k-1)
               enddo

               if (dbg) then
                  if (partial_create) print *,'Creating block ',k,
     *            ' indices ',
     *            (index_table(c_current_seg,ind(l)),l=1,nindex)
               endif

               call determine_current_block_size(ind, 
     *                    nindex, index_table, nindex_table,
     *                    segment_table, nsegment_table, size)

c----------------------------------------------------------------------------
c   Attempt to set up the new block.
c----------------------------------------------------------------------------

               block = get_block_number(result,k)
               if (block .le. 0) then

c----------------------------------------------------------------------------
c   Now allocate the block.
c----------------------------------------------------------------------------

                  ierr = allocate_block(result, k, size,
     *                            array_table, narray_table,
     *                            index_table, nindex_table,
     *                            block_map_table)
                  if (ierr .le. 0) then
                     print *,'Error: During create operation.'
                     print *,'Cannot allocate a block for ',
     *                  'array',result,' block number ',k,
     *                  ' on processor ',me
                     call array_block_summary(array_table, 
     *                                        narray_table)
                     call dump_block_ids()
                     call abort_job()
                  else
                     blkndx = ierr
                  endif

c--------------------------------------------------------------------------
c   Store the block's index and segment data in the blkmgr table.
c---------------------------------------------------------------------------

                  call blkmgr_insert_block_in_list(
     *                      array_table(c_block_list,result),dummy,
     *                      blkndx, c_block_list_ptr, .false.)
                  call set_block_indices(result, k, blkndx,
     *                                   array_table(1,result))
                  call set_block_segments(result, k, blkndx, 
     *                        index_table, nindex_table)
               else

c-------------------------------------------------------------------------
c   The block is left over from previous operations. 
c   Since the block was not allocated through the allocate_block routine
c   it must be busied.
c-------------------------------------------------------------------------

                 blkndx = block
                  call set_block_busy_flag(result, k, blkndx, 1)

                 call blkmgr_insert_block_in_list(
     *                     array_table(c_block_list,result), dummy,
     *                     blkndx, c_block_list_ptr, .false.) 
               endif  

c-------------------------------------------------------------------------
c   Set the block_computed_flag so the data will be preserved until 
c   a delete occurs.
c-------------------------------------------------------------------------

               call set_block_computed_flag(result, k, blkndx, 1)
               call set_block_created_flag(result,k, blkndx, 1)

c---------------------------------------------------------------------------
c   Store the blkndx in the block_map_table entry.
c---------------------------------------------------------------------------

               block_map_table(c_bmap_blkndx,iblock_map+k-1) = 
     *                              blkndx

c---------------------------------------------------------------------------
c   Clear the block.  Only the actual extent of the block's data is zeroed.
c---------------------------------------------------------------------------

               stack = array_table(c_array_stack,result)
               call clear_block(result, k, stack, blkndx, size)
            endif
  200    enddo 

c---------------------------------------------------------------------------

c   Set the create flag in the array_table entry.
c---------------------------------------------------------------------------

         array_table(c_create_flag, result) = 1

c----------------------------------------------------------------------------
c   Restore the saved segment values.
c----------------------------------------------------------------------------

         do k = 1, nindex
            index_table(c_current_seg, ind(k)) = segsave(k)
            array_table(c_index_array1+k-1,result) = ind_input(k)
         enddo
    
c---------------------------------------------------------------------------
c   Add the create data to the checkpoint data common block.
c---------------------------------------------------------------------------

         if (.not. restart_job) then
            do i = 1, nactive_create_table
               if (active_create_table(i) .eq. result .and.
     *          active_create_op(i) .eq. current_op) go to 300
            enddo

            nactive_create_table   = nactive_create_table + 1
            active_create_table(nactive_create_table) = result
            active_create_op(nactive_create_table)    = current_op
         endif
  300    continue

         return
      endif

c---------------------------------------------------------------------------
c   Delete operation.  Delete the blocks of a distributed array.
c---------------------------------------------------------------------------

      if (opcode .eq. delete_op) then
         if (result_type .ne. distributed_array) then
            print *,'Error: Delete operation called for array ',
     *               result,'.  Array is not distributed.'
            call abort_job()
         endif

c         call scrub_blocks(.false., array_table, narray_table, 
c     *                     index_table, nindex_table, 
c     *                     block_map_table, ierr)

         nb = array_table(c_numblks, result)
         iblock_map = array_table(c_block_map, result)
         do k = 1, nb
c            iblock = get_block_number(result, k)
            iblock = block_map_table(c_bmap_blkndx,iblock_map+k-1)
            if (iblock .gt. 0) then

c---------------------------------------------------------------------------
c   Check the block request to determine if there is an outstanding
c   GET or PUT still going on.
c---------------------------------------------------------------------------

               request = get_block_request(result, k, iblock)
               if (request .ne. MPI_REQUEST_NULL) 
     *            call wait_on_block(result, k, iblock,
     *                result_type, 
     *                request, instruction_timer, comm_timer)

               call clear_block_created_flag(result, k, iblock)
               call_marker = 20202
               call free_block(result, k, iblock, array_table,
     *                    narray_table, index_table, nindex_table,
     *                    block_map_table)
            endif   
         enddo 

c---------------------------------------------------------------------------
c   Clear the create flag in the array_table entry.
c---------------------------------------------------------------------------

         array_table(c_create_flag, result) = 0

c---------------------------------------------------------------------------
c   Remove the last create data for this array from the checkpoint_data 
c   common block.
c---------------------------------------------------------------------------

         do i = nactive_create_table, 1, -1
            if (active_create_table(i) .eq. result) then

c--------------------------------------------------------------------------
c   Shift remaining table data over this entry.
c--------------------------------------------------------------------------

               do j = i+1, nactive_create_table
                  active_create_table(j-1) = active_create_table(j)
                  active_create_op(j-1)    = active_create_op(j)
               enddo

               nactive_create_table = nactive_create_table - 1
               return
            endif
         enddo

         return
      endif   ! delete_op

c--------------------------------------------------------------------------
c   allocate_op: Allocate the (non-wildcard) blocks of a local array.
c--------------------------------------------------------------------------

      if (opcode .eq. allocate_op) then
         call allocate_instruction(array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op)
         return
      endif

c--------------------------------------------------------------------------
c   deallocate_op: Deallocate the (non-wildcard) blocks of a local array.
c--------------------------------------------------------------------------

      if (opcode .eq. deallocate_op) then
         call deallocate_instruction(array_table, narray_table,
     *                      index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      block_map_table, op, 
     *                      comm_timer, instruction_timer)
         return
      endif

c--------------------------------------------------------------------------
c   get_op: Fetch a block of a distributed array into a static copy.
c--------------------------------------------------------------------------

      if (opcode .eq. get_op) then
         call fetch_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, direct_flag)
         return
      endif

c---------------------------------------------------------------------------
c   request_op: Request a block of a "served" array.
c---------------------------------------------------------------------------

      if (opcode .eq. request_op) then
         call request_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer, 
     *                      direct_flag)

#ifdef BLUEGENE
         call get_prefetch_flag(prefetch_flag)
         if (prefetch_flag) then
            call pre_request_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer,
     *                      direct_flag)
         endif
#endif
         return
      endif

c----------------------------------------------------------------------------
c   Partial request (prequest).
c----------------------------------------------------------------------------

      if (opcode .eq. prequest_op) then
         call prequest_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer,
     *                      direct_flag)
         return
      endif

c----------------------------------------------------------------------------
c   prepare_op: Sends a static block of a "served" array to a server.
c----------------------------------------------------------------------------

      if (opcode .eq. prepare_op .or.
     *    opcode .eq. prepare_increment_op) then
         call prepare_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer)
         return
      endif

c----------------------------------------------------------------------------
c   Destroy: Frees blocks of a served array on the server side, allows reuse
c            of the array's disk space for other blocks.
c----------------------------------------------------------------------------

      if (opcode .eq. destroy_op) then
         call destroy_instruction(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, instruction_timer, comm_timer)
         return
      endif

c----------------------------------------------------------------------------
c   compute_integrals_op: Compute a block of integrals locally.
c----------------------------------------------------------------------------

      if (opcode .eq. compute_integrals_op) then
         call local_integral_block(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table,
     *                      op, direct_flag)
         return
      endif

      if (result_type .eq. static_array .or.
     *    result_type .eq. scalar_value) then
         
c--------------------------------------------------------------------------
c   The result array is "static".  Do not look for a block
c   since it is not managed by blkmgr.
c--------------------------------------------------------------------------
 
         block = -1
         go to 1000
      endif

c--------------------------------------------------------------------------
c    Must we do the computation of the result block?
c--------------------------------------------------------------------------

      if (block .le. 0) then

c----------------------------------------------------------------------------
c   Block is not present.  Create a new block on the fly.
c----------------------------------------------------------------------------

         call create_current_block(result,array_table,
     *                 narray_table, index_table,
     *                 nindex_table, segment_table, nsegment_table,
     *                 block_map_table, nblock_map_table, op, 
     *                 .true., .true., block, ierr)
         create_block = .true.
         if (ierr .le. 0) then
            print *,'Error in compute_block: Cannot create block for ',
     *            'array ',result
            call dump_block_ids()
            call abort_job()
         else
            blkndx = ierr
         endif
      endif

c---------------------------------------------------------------------------
c   Compute the block if one or more of the operands has been updated.
c---------------------------------------------------------------------------

 1000    continue
         op1    = op(c_op1_array)
         op2    = op(c_op2_array)
         op1_type = array_table(c_array_type,op1)
         op2_type = array_table(c_array_type,op2)

c---------------------------------------------------------------------------
c   Unpack the array_table data for the operands.
c---------------------------------------------------------------------------

         if (debug) print *,'op1, op2 = ',op1,op2
         if (op1_type .eq. static_array .or.
     *       op1_type .eq. scalar_value) then
            op1_block = -1 
            flag1     = 0
         else
            op1_block = find_current_block(op1, array_table(1,op1),
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           block_map_table, op1_blkndx)
            if (op1_block .le. 0) then
               print *,'Task ',me,' Error: line ',current_line,
     *                 ' op1 array = ',
     *                op1,' block = ',op1_block
               print *,'Cannot compute result unless blocks exist.'
               print *,'optable entry: ',(op(j),j=1,loptable_entry)
               print *,'Array table entry for ',op1,' is ',
     *           (array_table(j,op1),j=1,larray_table_entry)
               call dump_array_list(op1, array_table, narray_table,
     *              'DUMP OF ARRAY BLOCKS LIST: ')
               call array_block_summary(array_table, narray_table)
               print *,'Current indices:'
               do i = 1, nindex_table
                  print *,'   index ',i,' current segment ',
     *              index_table(c_current_seg,i)
               enddo
               call dump_block_ids()
               call abort_job()
            endif

            call get_block_computed_flag(op1, op1_block, op1_blkndx,
     *                                   flag1)

c--------------------------------------------------------------------------
c   Check for possible block communication in progress, wait if necessary.
c--------------------------------------------------------------------------

            request = get_block_request(op1, op1_block, op1_blkndx)
            if (request .ne. mpi_request_null) 
     *          call wait_on_block(op1, op1_block, op1_blkndx,
     *              op1_type, request,
     *             instruction_timer, comm_timer)
         endif
         if (debug) print *,'op1_block,flag1 = ',op1_block, flag1
        
         if (op2_type .eq. static_array .or.
     *       op2_type .eq. scalar_value .or.
     *       opcode .eq. assignment_op) then
            op2_block = -1
            flag2     = 0
         else
            op2_block = find_current_block(op2, array_table(1,op2),
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           block_map_table, op2_blkndx)
            if (op2_block .le. 0) then
               print *,'Error: op2 array = ',op2,' block = ',op2_block
               print *,'Cannot compute result unless both blocks exist.'
               print *,'optable entry: ',(op(j),j=1,loptable_entry)
               call array_block_summary(array_table, narray_table)
               print *,'Current indices:'
               do i = 1, nindex_table
                  print *,'   index ',i,' current segment ',
     *              index_table(c_current_seg,i)
               enddo
               call dump_block_ids()
               call abort_job()
            endif

            call get_block_computed_flag(op2, op2_block, op2_blkndx,
     *                                    flag2)

c--------------------------------------------------------------------------
c   Check for possible block communication in progress, wait if necessary.
c--------------------------------------------------------------------------

            request = get_block_request(op2, op2_block, op2_blkndx)
            if (request .ne. mpi_request_null) 
     *          call wait_on_block(op2, op2_block, op2_blkndx,
     *             op2_type, request,
     *             instruction_timer, comm_timer)
         endif

         if (debug) print *,'op2_block,flag2 = ',op2_block, flag2

c-----------------------------------------------------------------------
c   If neither of the operands were computed on this pass, there is no
c   need to compute the block. 
c-----------------------------------------------------------------------

         if (.not. create_block) then
            if (flag1 .eq. 0 .and. flag2 .eq. 0) then

c--------------------------------------------------------------------------
c   Check for a scalar-scalar --> scalar operation
c--------------------------------------------------------------------------

               if (result_type .eq. scalar_value .and.
     *             array_table(c_array_type,op1) .eq. scalar_value .and.
     *             array_table(c_array_type,op2) .eq. scalar_value) then
                  call handle_sss_op(op, array_table, narray_table,
     *                      scalar_table, nscalar_table)
                  return
               else if (result_type .eq. scalar_value .and.
     *                  array_table(c_array_type,op1) .eq. 
     *                                       scalar_value .and.
     *                  opcode .eq. assignment_op) then
                  call handle_sss_op(op, array_table, narray_table,
     *                      scalar_table, nscalar_table)
                  return
               endif
            endif
         endif

         if (opcode .eq. contraction_op) then

c---------------------------------------------------------------------------
c   Contraction.
c---------------------------------------------------------------------------

            ncind = 0
            do i = 1, mx_array_index
               cind(i) = op(c_ind1+i-1)
               if (cind(i) .ne. 0) ncind = ncind + 1
            enddo
       
c            print *,'Contraction: ',(op(k),k=1,loptable_entry) 
            call contract_blocks(op1, op1_block, op1_blkndx,
     *                           op2, op2_block, op2_blkndx,
     *                           result, block, blkndx, cind, ncind,
     *                           array_table, narray_table,
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           scalar_table, nscalar_table,
     *                           address_table, block_map_table,
     *                           flopcount)
         else if (opcode .eq. tensor_op) then
            call tensor_contraction(op1, op1_block, op1_blkndx,
     *                           op2, op2_block, op2_blkndx,
     *                           result, block, blkndx,
     *                           array_table, narray_table,
     *                           index_table, nindex_table,
     *                           segment_table, nsegment_table,
     *                           scalar_table, nscalar_table,
     *                           address_table, flopcount)
         else if (opcode .eq. sum_op .or.
     *            opcode .eq. subtract_op) then

c---------------------------------------------------------------------------
c   Summation.
c---------------------------------------------------------------------------

            call sum_blocks(op1, op1_block, op1_blkndx,
     *                      op2, op2_block, op2_blkndx,
     *                      result, block, blkndx, opcode, array_table,
     *                      narray_table, index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
         else if (opcode .eq. assignment_op) then
            call assign_block(op1, op1_block, op1_blkndx, 
     *                      result, block, blkndx,
     *                         array_table,
     *                      narray_table, index_table, nindex_table,
     *                      segment_table, nsegment_table,
     *                      scalar_table, nscalar_table,
     *                      address_table, flopcount)
         else
            print *,'Error: Invalid opcode in compute_block.' 
            print *,'       optable entry is ',
     *                (op(i),i=1,loptable_entry)
            call abort_job()
         endif

c--------------------------------------------------------------------------
c   Set the "block_computed_flag" with the blkmgr.
c--------------------------------------------------------------------------

      if (result_type .ne. scalar_value .and.
     *    result_type .ne. static_array) then
         call get_block_computed_flag(result, block, blkndx,
     *                                result_flag)
         if (result_flag .eq. 0) then
            call get_block_created_flag(result, block, blkndx, 
     *                                  result_flag)
            if (result_flag .eq. 0) then
               call set_opblock(result, block, blkndx, op)
               call set_block_computed_flag(result, block, blkndx, 1)
            endif
         endif
      endif

      return
      end
      
      subroutine check_create_flag(op, array_table, narray_table)
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'trace.h'
      integer narray_table
      integer op(loptable_entry)
      integer array_table(larray_table_entry, narray_table)

      integer i, ierr
      integer result, result_type
      integer op1, op1_type, op2, op2_type

      if (op(c_opcode) .eq. create_op .or.
     *    op(c_opcode) .eq. reindex_op .or.
     *    (op(c_opcode) .ge. sp_add_op .and. 
     *     op(c_opcode) .le. sp_ldindex_op)) return

      ierr = 0
      result = op(c_result_array)
      if (result .gt. 0) then
         result_type = array_table(c_array_type, result)

         if (result_type .eq. distributed_array) then
            if (array_table(c_create_flag,result) .eq. 0) then
              print *,'Task ',me,
     *              ' Error: Attempt to use a distributed ',
     *              'array before it has been created.'
              print *,'Line ',current_line,' array ',result
              print *,'Operation = ',(op(i),i=1,loptable_entry)
              print *,'Result array = ',result,' has not been created.'
              ierr = 1
            endif
         endif
      endif

      op1 = op(c_op1_array)
      if (op1 .gt. 0) then
         op1_type = array_table(c_array_type,op1)
         if (op1_type .eq. distributed_array) then
            if (array_table(c_create_flag,op1) .eq. 0) then
               print *,'Task ',me,
     *               ' Error: Attempt to use a distributed ',
     *              'array before it has been created.'
               print *,'Line ',current_line,' array ',op1
               print *,'Operation = ',(op(i),i=1,loptable_entry)
               print *,'Operand 1 array = ',op1,' has not been created.'
               ierr = 1
            endif
         endif
      endif

      op2 = op(c_op2_array)
      if (op2 .gt. 0) then
         op2_type = array_table(c_array_type,op2)
         if (op2_type .eq. distributed_array) then
            if (array_table(c_create_flag,op2) .eq. 0) then
               print *,'Task ',me,
     *                 ' Error: Attempt to use a distributed ',
     *              'array before it has been created.'
               print *,'Line ',current_line,' array ',op2
               print *,'Operation = ',(op(i),i=1,loptable_entry)
               print *,'Operand 1 array = ',op2,' has not been created.'
               ierr = 1
            endif
         endif
      endif

      if (ierr .eq. 1) then
         print *,'Error at line number ',current_line
         call c_flush_stdout()
         call abort_job()
      endif 
      
      return
      end

      subroutine array_block_summary(array_table, narray_table)
c-------------------------------------------------------------------------
c   Prints a summary of the blocks in use for each array in the array_table.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      integer get_block_array_handle

      integer i, j, nblks
      integer ierr

      do i = 1, nblkmgr_stacks
         print *,'Stack ',i,' start block ',stack_start(i),
     *      ' end block ',
     *      stack_start(i)+stack_size(i)-1,' blocksize ',
     *      stack_blocksize(i),' stack pointer ',free_stack_ptr(i)
      enddo

      print *,'ARRAY BLOCK SUMMARY FOR TASK ',me
      print *,'Array    Type     Total blks    Blks in use '
      print *,'-----    -----    ----------    -----------'

      do i = 1, narray_table
         if (array_table(c_array_type,i) .ne. scalar_value) then
            nblks = 0
            do j = 1, blkmgr_blocks
               if (get_block_array_handle(j) .eq. i) 
     *            nblks = nblks + 1
            enddo

            if (nblks .gt. 0) then
               print 100,i,array_table(c_array_type,i),
     *                array_table(c_numblks,i),nblks
            endif
         endif
      enddo

  100 format(i5,6x,i4,6x,i5,5x,i6)

      return
      end      

      subroutine dump_array_list(array, array_table, narray_table,
     *                           msg)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array
      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      character*(*) msg
      integer next, i, x, nind, flag
      integer seg(mx_array_index)
      integer blkndx
      integer get_block_array_handle

      next = array_table(c_block_list,array)
      nind = array_table(c_nindex, array)
      print *,msg
      print *,'Blocks for array ',array
  100 continue
      if (next .eq. 0) return
      blkndx = next
      x = get_block_array_handle(blkndx)
      call get_block_segments(blkndx, seg)
      call get_blk_header(flag, blkndx,c_block_flags)
      print 200,x,blkndx,flag, (seg(i),i=1,nind)
  200 format(' array, blkndx ',2i5,' flags ',i8,' segs ',6(i3,1x))
      call get_blk_header(next, blkndx,c_block_list_ptr) 
      go to 100

      return
      end
