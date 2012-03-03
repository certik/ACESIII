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
      subroutine allocate_instruction(array_table, narray_table,
     *                      index_table,
     *                      nindex_table, segment_table, nsegment_table,
     *                      block_map_table, nblock_map_table, op)
c---------------------------------------------------------------------------
c   Runtime implementation code for the ALLOCATE instruction.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'checkpoint_data.h'

      integer narray_table, nindex_table, nsegment_table
      integer array_table(larray_table_entry,narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer segment_table(lsegment_table_entry,nsegment_table)
      integer nblock_map_table
      integer block_map_table(lblock_map_entry, nblock_map_table)
      integer op(loptable_entry)

      integer array, type, size, nind
      integer i, iblock, nblock, nwild
      integer blkndx
      integer ierr, block, dummy
      integer allocate_block
      integer ind(mx_array_index)
      integer ind_save(mx_array_index)
      integer seg(mx_array_index)
      integer nseg(mx_array_index)
      integer bseg(mx_array_index)
      integer eseg(mx_array_index)
      integer iwild(mx_array_index)
      integer iseg(mx_array_index)
      integer stack

      array = op(c_result_array)
      type  = array_table(c_array_type,array)
      nind = array_table(c_nindex, array)
      if (type .ne. local_array) then
         print *,'Error: Allocate instruction requires a local_array'
         print *,'Array, type = ',array,type
         call abort_job()
      endif

c      print *,'Task ',me,' ALLOC at line ',current_line,' op ',
c     *   (op(c_ind1+i-1),i=1,mx_array_index)

c----------------------------------------------------------------------------
c   Switch the array's indices to the indices used in its SIAL reference,
c   which are the indices passed in on the current instruction.
c   Wildcard indices use the indices used in the original array declaration.
c----------------------------------------------------------------------------

      do i = 1, nind
         ind_save(i) = array_table(c_index_array1+i-1,array) 
         if (op(c_ind1+i-1) .eq. wildcard_indicator) then ! wild card index
            array_table(c_index_array1+i-1,array) =
     *             array_table(c_index_original+i-1,array)
         else
            array_table(c_index_array1+i-1,array) = op(c_ind1+i-1)
         endif
      enddo

c----------------------------------------------------------------------------
c   Determine the indices to use.
c----------------------------------------------------------------------------

      nblock = 1
      nwild  = 0
      do i = 1, nind
         ind(i) = array_table(c_index_array1+i-1,array)
         if (ind(i) .lt. 1 .or. ind(i) .gt. nindex_table) then
            print *,'Task ',me,' Error in ALLOCATE at line ',
     *           current_line
            print *,'Invalid index: ',ind(i)
            call abort_job()
         endif

         seg(i) = index_table(c_current_seg,ind(i))
         nseg(i) = index_table(c_nsegments, ind(i))
         bseg(i) = index_table(c_bseg, ind(i))
         eseg(i) = index_table(c_eseg, ind(i))

         if (op(c_ind1+i-1) .eq. wildcard_indicator) then
            nblock = nblock*nseg(i)
            nwild = nwild + 1
            iwild(nwild) = i
            iseg(iwild(nwild)) = bseg(iwild(nwild)) 
         endif
      enddo

c      print *,'Task ',me,' iwild ',(iwild(i),i=1,nwild)
c      print *,'Task ',me,' iseg ',(iseg(i),i=1,mx_array_index)

c VFL skip to allow nwild = 0 
c     if (nwild .eq. 0) then
c        print *,'Error: No wildcard index in allocate at line ',
c    *          current_line
c        call abort_job()
c     endif
c ENDVFL skip to allow nwild = 0 

      do iblock = 1, nblock
         
c---------------------------------------------------------------------------
c   Set the "current segment" field in the index table, for each wildcard
c   index.
c--------------------------------------------------------------------------

         if (nwild .gt. 0) then 
         do i = 1, nwild
            index_table(c_current_seg,ind(iwild(i))) = iseg(iwild(i))
         enddo
         endif 

c-------------------------------------------------------------------------
c   Determine the blocksize to be allocated based on the index and segment 
c   data.
c-------------------------------------------------------------------------

c         print 10101,me,iblock,(ind(i),
c     *          index_table(c_current_seg,ind(i)),i=1,nind)
c10101 format(' Task ',i2,' iblock ',i3,' ind, curseg: ',
c     *      4(i3,1x,i6,1x))
         call determine_current_block_size(ind, nind,
     *        index_table, nindex_table,
     *        segment_table, nsegment_table, size)

c----------------------------------------------------------------------------
c   Release any blocks that are eligible to be freed up.
c----------------------------------------------------------------------------

c         call release_blocks(array_table, narray_table,
c     *                       index_table, nindex_table,
c     *                       block_map_table, nblock_map_table)

c----------------------------------------------------------------------------
c   Attempt to set up the new block.
c----------------------------------------------------------------------------

         ierr = allocate_block(array, iblock, size, array_table, 
     *                         narray_table, index_table, 
     *                         nindex_table, block_map_table)
         if (ierr .le. 0) then
            print *,'Error: During allocate operation.'
            print *,'Cannot allocate a block for ',
     *            'array',array,' block number ',iblock,
     *            ' on processor ',me
            call array_block_summary(array_table,
     *                               narray_table)
            call dump_block_ids()
            call abort_job()
         else
            blkndx = ierr
	    endif

         call blkmgr_insert_block_in_list(
     *                 array_table(c_block_list,array),
     *                 dummy, blkndx, c_block_list_ptr, .false.)
         call set_block_indices(array, iblock, blkndx,
     *                          array_table(1,array))
         call set_block_segments(array, iblock, blkndx, index_table,
     *                      nindex_table)

c-------------------------------------------------------------------------
c   Set the block_computed_flag so the data will be preserved until
c   a delete occurs.
c-------------------------------------------------------------------------

         call set_block_computed_flag(array, iblock, blkndx, 1)
         call set_block_created_flag(array, iblock, blkndx, 1)

c---------------------------------------------------------------------------
c   Clear the block.  Only the actual extent of the block's data is zeroed.
c---------------------------------------------------------------------------

         stack = array_table(c_array_stack,array)
         call clear_block(array, iblock, stack, blkndx, size)

c--------------------------------------------------------------------------
c   Increment the "wildcard" segments.
c--------------------------------------------------------------------------

         if (nwild .gt. 0) then 
         iseg(iwild(1)) = iseg(iwild(1)) + 1
         do i = 2, nwild
            if (iseg(iwild(i-1)) .le. eseg(iwild(i-1))) go to 100
            iseg(iwild(i-1)) = bseg(iwild(i-1))
            iseg(iwild(i))   = iseg(iwild(i)) + 1
         enddo
  100    continue
         endif 

      enddo

c-----------------------------------------------------------------------------
c   Restore the "current segment" fields of the index_table.
c   Restore the array's indices to those at the beginning of the instruction.
c-----------------------------------------------------------------------------

      do i = 1, nind
         index_table(c_current_seg,ind(i)) = seg(i)
c         if (ind_save(i) .gt. 0)
c     *       array_table(c_index_array1+i-1,array) = ind_save(i)
      enddo

c-------------------------------------------------------------------------
c   Remember this allocate for the checkpoint data.
c-------------------------------------------------------------------------

      if (.not. restart_job) then
         do i = 1, nactive_allocate_table
            if (active_allocate_table(i) .eq. array .and. 
     *          active_allocate_op(i) .eq. current_op) go to 200
         enddo      

         nactive_allocate_table   = nactive_allocate_table + 1
         active_allocate_table(nactive_allocate_table) = array
         active_allocate_op(nactive_allocate_table)    = current_op 
      endif
  200 continue

      return
      end

      
