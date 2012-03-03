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
      subroutine blkmgr_init(x, nx, nblocks_total,
     *                       blk_header,
     *                       free_stack,
     *                       blk_addr,
     *                       blk_in_use,
     *                       blksize, stack_blocks,
     *                       nstacks, comm)
c------------------------------------------------------------------------
c   Initialization of block manager subsystem.
c   x         = array to use for block management.
c   nx        = total number of double precision words in x
c   blksize   = array of blocksizes to use in each stack of blocks.
c   stack_blocks = array of blocks used for each stack.
c   nstacks   = number of stacks
c   comm      = company communicator.
c-------------------------------------------------------------------------  

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'machine_types.h'
      include 'dbugcom.h'
 
      integer nblocks_total
      integer blk_header(nblocks_total,lblk_header)
      integer free_stack(nblocks_total) 
      integer blk_in_use(nblocks_total) 
      integer*8 blk_addr(nblocks_total)

      integer nstacks
      integer i, j, nx, blksize(nstacks), stack_blocks(nstacks)
      double precision x(*)
      integer comm
      integer*8  count_word, nword_stack
      integer*8 c_loc64, ixx
      integer count_block
      integer size, ierr
      integer nblocks_stack, wsum
      integer*8 indblk
      integer stk_weights(10)
      integer*8 addr
      integer k
      integer iii

      blkmgr_comm   = comm
      blkmgr_next_block = 0
      call mpi_comm_rank(comm, blkmgr_comm_rank, ierr)

      ixx = 1
      addr_blk_header = c_loc64(blk_header, ixx, intsize) 
      addr_free_stack = c_loc64(free_stack, ixx, intsize)
      addr_blk_addr   = c_loc64(blk_addr,   ixx, 8)
      addr_blk_in_use = c_loc64(blk_in_use, ixx, intsize)

c---------------------------------------------------------------------------
c   Set up the stack blocks.
c---------------------------------------------------------------------------

      count_block = 0
      count_word  = 0
      do i = 1, nstacks
         stack_size(i) = stack_blocks(i)
         count_block = count_block + stack_size(i)
         count_word = count_word + stack_size(i) * 
     *                  (blksize(i) + lblock_id_data)
      enddo
      
c--------------------------------------------------------------------------
c   Inform the user about the actual stack sizes.
c--------------------------------------------------------------------------

      if (dbg) then
         do i = 1, nstacks
            print *,'Task ',blkmgr_comm_rank,' Stack ',i,' contains ',
     *        stack_size(i),' blocks of ',blksize(i)+lblock_id_data,
     *        ' words per block.'
         enddo
      endif

      blkmgr_blocks = count_block
      if (blkmgr_blocks .gt. nblocks_total) then
         print *,'Error: Total number of blocks (',blkmgr_blocks,
     *     ') exceeds the blkmgr limit of ',nblocks_total
         call abort_job()
      endif

      if (count_word .gt. nx) then
         print *,'Error: Total number of words (',count_word,
     *       ') exceeds amount allowed for blkmgr data (',
     *       nx,')'
         call abort_job()
      endif

      do j = 1, lblk_header
      do i = 1, blkmgr_blocks
        call put_blk_header(0, i, j)
      enddo
      enddo

c-----------------------------------------------------------------------
c   Set up the block stacks.
c-----------------------------------------------------------------------

      nblkmgr_stacks = nstacks
      count_block = 1
      count_word  = 1
      
      do i = 1, nblkmgr_stacks
         stack_start(i)     = count_block
         free_stack_ptr(i)  = stack_start(i) - 1
         stack_base_addr(i) = c_loc64(x, count_word, bytes_per_double)
         stack_blocksize(i) = blksize(i) + lblock_id_data
         stack_datasize(i)  = blksize(i)
         
c--------------------------------------------------------------------------
c   Each block is pushed onto a free stack.
c   Also, calculate the block address.
c--------------------------------------------------------------------------

         do j = 1, stack_size(i)
            call push_free_stack(i, count_block) 
            call put_blk_header(i, count_block,c_block_stack)
            addr = c_loc64(x, count_word, bytes_per_double)
            call put_blk_addr(addr, count_block)
            count_block = count_block + 1
            count_word  = count_word + stack_blocksize(i)
         enddo
      enddo

      last_blk_addr = blk_addr(count_block-1)
      first_blk_addr = blk_addr(1)

c--------------------------------------------------------------------------
c   Adjust block and word counters.
c--------------------------------------------------------------------------

      count_block = count_block - 1
      count_word  = count_word - 1
      
      nblk_in_use = 0
      do i = 1, blkmgr_blocks
         blk_in_use(i) = 0
      enddo

      do i = 1, blkmgr_blocks
         call put_blk_header(mpi_request_null, i, c_block_request)
      enddo
 
      comm_list_head = 0
      comm_list_tail = 0
      persistent_list_head = 0
      persistent_list_tail = 0

      blkmgr_init_flag = .true.
      if (dbg) then
         print *,'Task ',blkmgr_comm_rank,
     *        ' BLKMGR has been initialized with ',
     *        count_word,' double precision ',
     *        'words, or ',count_block,' blocks.'
         print *,'Task ',blkmgr_comm_rank,
     *     ' stack base_addr: ',(stack_base_addr(i),i=1,nblkmgr_stacks)
      endif
      return
      end

      subroutine put_blk_in_use(blkndx)
c----------------------------------------------------------------------------
c   Adds a block to the "blk_in_use" table.
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer blkndx

      integer*8 get_index_from_base, ixx
      integer blk_in_use(1)
#ifdef ALTIX
      pointer (iptr, blk_in_use)
#else
      common blk_in_use
#endif

#ifdef ALTIX 
      iptr = ishptr
#endif

      ixx = get_index_from_base(addr_blk_in_use, blk_in_use, 1)
      if (nblk_in_use .ge. blkmgr_blocks) then
         print *,'Error: Attempting to set too many blocks in use.'
         call abort_job()
      endif

      nblk_in_use = nblk_in_use + 1
      blk_in_use(ixx + nblk_in_use - 1) = blkndx
      return
      end

      subroutine remove_block_from_use(blkndx)
c----------------------------------------------------------------------------
c   Removes an entry from the "blk_in_use" table. 
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer blkndx

      integer i
      integer*8 get_index_from_base, ixx
      integer blk_in_use(1)
#ifdef ALTIX
      pointer (iptr, blk_in_use)
#else 
      common blk_in_use
#endif
      
#ifdef ALTIX
      iptr = ishptr
#endif
      ixx = get_index_from_base(addr_blk_in_use, blk_in_use, 1)
      do i = 1, nblk_in_use
         if (blk_in_use(ixx+i-1) .eq. blkndx) then
            
c--------------------------------------------------------------------------
c   Remove entry by swapping entry "nblk_in_use" over it.
c--------------------------------------------------------------------------

            blk_in_use(ixx+i-1) = blk_in_use(ixx+nblk_in_use-1)
            nblk_in_use = nblk_in_use - 1
            return
         endif
      enddo

      print *,'Error: Attempted to remove ',blkndx,' from use, ',
     *          'but it was not registered as in use already.'
      call abort_job()
      return
      end

      subroutine put_blk_header(data, blkndx, col)
c----------------------------------------------------------------------------
c   Stores a word in the blk_header table at row "blkndx", cloumn "col".
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer data, blkndx, col

      integer*8 get_index_from_base, ixx
      integer blk_header(1)
#ifdef ALTIX
      pointer (iptr, blk_header)
#else
      common blk_header
#endif

#ifdef ALTIX
      iptr = ishptr  
#endif
      ixx = get_index_from_base(addr_blk_header, blk_header, 1)
      blk_header(ixx + (col-1)* blkmgr_blocks + blkndx - 1) = data 
      return
      end

      subroutine get_blk_header(data, blkndx, col)
c----------------------------------------------------------------------------
c   Retrieves a word in the blk_header table from row "blkndx", column "col".
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer data, blkndx, col

      integer*8 get_index_from_base, ixx
      integer blk_header(1)
#ifdef ALTIX
      pointer (iptr, blk_header) 
#else
      common blk_header
#endif

#ifdef ALTIX
      iptr = ishptr
#endif

      ixx = get_index_from_base(addr_blk_header, blk_header, 1)
      data = blk_header(ixx + (col-1)* blkmgr_blocks + blkndx - 1)
      return
      end

      subroutine put_free_stack(data, blkndx)
c----------------------------------------------------------------------------
c   Stores a word in the free_stack table at row "blkndx".
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer data, blkndx

      integer*8 get_index_from_base, ixx
      integer free_stack(1)
#ifdef ALTIX
      pointer (iptr, free_stack)
#else
      common free_stack
#endif
 
#ifdef ALTIX
      iptr = ishptr
#endif
      if (data .eq. 0) print *,'Task ',me,' put_free_stack: put data ',
     *     data,' into blkndx ',blkndx
      ixx = get_index_from_base(addr_free_stack, free_stack, 1)
      free_stack(ixx + blkndx - 1) = data 

      return
      end

      subroutine get_free_stack(data, blkndx)
c----------------------------------------------------------------------------
c   Retrieves a word from the free_stack table from row "blkndx".
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer data, blkndx

      integer*8 get_index_from_base, ixx
      integer free_stack(1)
#ifdef ALTIX
      pointer (iptr, free_stack)
#else
      common free_stack
#endif

#ifdef ALTIX
      iptr = ishptr
#endif  
      ixx = get_index_from_base(addr_free_stack, free_stack, 1)
      data = free_stack(ixx + blkndx - 1)
      return
      end

      subroutine put_blk_addr(data, blkndx)
c----------------------------------------------------------------------------
c   Stores a word in the blk_addr table at row "blkndx".
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer*8 data
      integer blkndx
 
      integer*8 get_index_from_base, ixx
      integer*8 blk_addr(1)
#ifdef ALTIX
      pointer (i8ptr, blk_addr) 
#else
      common blk_addr
#endif
 
#ifdef ALTIX
      i8ptr = i8shptr
#endif
      ixx = get_index_from_base(addr_blk_addr, blk_addr, 3)
      blk_addr(ixx + blkndx - 1) = data 
      return
      end

      subroutine get_blk_addr(data, blkndx)
c----------------------------------------------------------------------------
c   Retrieves a word from the blk_addr table from row "blkndx".
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      include 'trace.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer*8 data
      integer blkndx

      integer i, j
      integer*8 get_index_from_base, ixx
      integer*8 blk_addr(1)
#ifdef ALTIX
      pointer (i8ptr, blk_addr)
#else
      common blk_addr
#endif

#ifdef ALTIX
      i8ptr = i8shptr
#endif
      ixx = get_index_from_base(addr_blk_addr, blk_addr, 3)
      data = blk_addr(ixx + blkndx - 1)

c      if (data .lt. first_blk_addr .or.
c     *    data .gt. last_blk_addr) then
c          print *,'Task ',me,' GET_BLK_ADDR: invalid address'
c          print *,'addr for block ',blkndx,' is ',data
c          print *,'first_blk_addr ',first_blk_addr,
c     *            ' last_blk_addr ',last_blk_addr
c          do i = 1, blkmgr_blocks
c             print *,'address for block ',i,' is ',blk_addr(ixx+i-1)
c          enddo
c          call abort_job()
c      endif   

      return
      end

      block data blkmgr_data
c--------------------------------------------------------------------------
c   Block data subprogram to initialize the "blkmgr_init_flag" variable.
c--------------------------------------------------------------------------
      include 'maxdim.h'
      include 'blkmgr.h'
      data blkmgr_init_flag/.false./
      end

      subroutine set_block_id(array, array_block, stack, block)
c-------------------------------------------------------------------------
c   Encodes the array descriptor and array_block identifier into the
c   block_id_data portion of a block in stack "stack".
c-------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'machine_types.h'
       include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer array, array_block, stack, block
 
      integer ix(2)
#ifdef ALTIX
      pointer (iptr, ix)
#else
      common ix
#endif

      integer*8 ixblk, c_loc64, ixx
      integer*8 get_index_from_base
      integer*8 addr

#ifdef ALTIX
      iptr = ishptr
#endif

c--------------------------------------------------------------------------
c   Find the integer index of ix relative to x(indblk).
c   This enables us to reference the data as integer.
c--------------------------------------------------------------------------

      call get_blk_addr(addr, block)
      ixblk = get_index_from_base(addr, ix, 1)

c---------------------------------------------------------------------------
c   Encode the integer data into the block.
c---------------------------------------------------------------------------

      ix(ixblk) = array
      ix(ixblk+1) = array_block

      return
      end

      integer function get_block_array_handle(blkndx)
c-------------------------------------------------------------------------
c   Returns the array handle corresponding to blk number "blkndx".
c-------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
 
      integer blkndx
      integer data

      call get_blk_header(data, blkndx,c_array_handle)
      get_block_array_handle = data
      return
      end

      subroutine get_block_id(stack, block, array, array_block)
c-------------------------------------------------------------------------
c   Returns the array descriptor and array_block identifier from the
c   block_id_data portion of a block in stack "stack".
c-------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer array, array_block, stack, block
 
      integer ix(2)
#ifdef ALTIX
      pointer (iptr, ix)
#else
      common ix
#endif
      integer i   
      integer*8 ixblk
      integer*8 c_loc64, ixx
      integer*8 get_index_from_base
      integer*8 addr
      integer size

#ifdef ALTIX
      iptr = ishptr
#endif

       call get_blk_addr(addr, block)
       ixblk = get_index_from_base(addr, ix, 1)

c---------------------------------------------------------------------------
c   Return the integer data of the block.
c---------------------------------------------------------------------------

      array = ix(ixblk)
      array_block = ix(ixblk+1)
      return
      end

      subroutine validate_block(data, array, block, flag)
c--------------------------------------------------------------------------
c   Checks the block_id_data portion of the block against the expected 
c   array and block info.  If it matches, the subroutine returns, else  
c   the job is terminated (if the abort flag is set).
c--------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'machine_types.h'

      integer data(*)
      integer array, block
      integer*8 ixx, c_loc64
      logical flag
 
      if (data(1) .ne. array .or.
     *    data(2) .ne. block) then
         print *,'Task ',me,' ERROR: Invalid data in block, line ',
     *                   current_line
         print *,'       Expected array ',array,' block ',block
         print *,'       Data had array id ',data(1),' block id ',
     *                   data(2)
         ixx = 1
         print *,'       Data address of block ',
     *        c_loc64(data, ixx, intsize)
         if (flag) then
            call dump_block_ids()
            call abort_job()
         endif
      endif
      return
      end

      subroutine free_block(array_handle, array_block, blkndx,
     *                 array_table, narray_table, 
     *                 index_table, nindex_table, block_map_table)
c---------------------------------------------------------------------------
c   Marks the block associated with array "array_handle", block "array_block"
c   as free.  
c----------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'

      integer i, flag
      integer array_handle, array_block, stack, blkndx
      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer block_map_table(lblock_map_entry, *)
      integer dummy
      integer iseg(mx_array_index)
      integer nind, lookup, block_map_lookup

      call blkmgr_remove_block_from_list(
     *             array_table(c_block_list,array_handle), 
     *              dummy, blkndx, c_block_list_ptr)

      call get_blk_header(flag, blkndx,c_block_flags)
      if (and(flag, block_persistence_flag) .ne. 0) then
         call blkmgr_remove_block_from_list(persistent_list_head,
     *      persistent_list_tail, blkndx, c_persistent_ptr)
      endif

c---------------------------------------------------------------------------
c   If this is a distributed array, its entry must be removed from the
c   block_map_table.
c---------------------------------------------------------------------------

      if (array_table(c_array_type,array_handle) .eq. 
     *                                 distributed_array) then
                                          
         nind = array_table(c_nindex,array_handle)
         do i = 1, nind
            call get_blk_header(iseg(i), blkndx,c_block_segments+i-1)
         enddo
                                                                                
         lookup = block_map_lookup(iseg, nind, array_handle,
     *                          array_table(1,array_handle),
     *                          index_table, nindex_table)
                                                                                
         if (block_map_table(c_bmap_blkndx, lookup) .ne. blkndx) then
            print *,'Error: blkndx lookup value is ',
     *         block_map_table(c_bmap_blkndx, lookup),
     *         ' which does not ',
     *         'match blkndx ',blkndx,' for array ',array_handle,
     *         ' block ',array_block,'. Lookup index is ',
     *          lookup
            call array_block_summary(array_table, narray_table)
            call dump_block_ids()
            call abort_job()
         endif

         block_map_table(c_bmap_blkndx, lookup) = 0
      endif

c---------------------------------------------------------------------------
c   Put this block back on its proper free stack.
c--------------------------------------------------------------------------

       call get_blk_header(stack, blkndx,c_block_stack)
       call push_free_stack(stack, blkndx)
       call remove_block_from_use(blkndx)

c---------------------------------------------------------------------------
c   Destroy the block data info.
c---------------------------------------------------------------------------

       do i = 1, lblk_header
          call put_blk_header(0, blkndx,i)
       enddo

       call put_blk_header(MPI_REQUEST_NULL, blkndx,c_block_request)
       call put_blk_header(stack, blkndx,c_block_stack)

      return
      end

      integer*8 function get_block_index(array, block, stack, blkndx,
     *                                   x, data_flag)
c--------------------------------------------------------------------------
c   Returns the index (in double words) of the block id info for block 
c   "block" of
c   array "array".  The index is relative to the argument "x".
c   If the "data_flag" argument is .true., the index of the actual data in 
c   the block is returned, otherwise, the index of the block_id area is 
c   returned.
c
c   If the block is not in the blkmgr table, a negative value is returned.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'machine_types.h'

      integer array, block, stack, blkndx
      logical data_flag

      double precision x(*)
      integer i, j, k, size
      integer*8 block_index, c_loc64, ixx
      integer*8 get_index_from_base
      integer*8 addr
      integer blk_in_use(1)
      integer handle, myblock

      if (blkndx .le. 0) then
         ixx = get_index_from_base(addr_blk_in_use, blk_in_use, 1)

c--------------------------------------------------------------------------
c   Scan all blocks that are registered as "in use".
c--------------------------------------------------------------------------

         do k = 1, nblk_in_use
            j = blk_in_use(ixx+k-1)

            call get_blk_header(handle, j,c_array_handle)
            call get_blk_header(myblock, j,c_array_block)
            if (handle .eq. array .and.
     *          myblock .eq. block) then
               i = j
               go to 100
            endif
         enddo
 
         get_block_index = -1
         return
      else
         i = blkndx
      endif

100   continue
      call get_blk_addr(addr, i)
      block_index = get_index_from_base(addr, x, 2)

      if (data_flag) then
         get_block_index = block_index +lblock_id_data
      else
         get_block_index = block_index
      endif
      
      return
      end 

      integer*8 function get_block_data_index(array, block, stack,
     *                                        blkndx, x)
c--------------------------------------------------------------------------
c   Returns the index (in double words) of the data for block "block" of
c   array "array".  The index is relative to the argument "x".
c
c   If the block is not in the blkmgr table, a negative value is returned.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'

      integer array, block, stack, blkndx
      integer*8 block_index
      double precision x(*)

      integer*8 get_block_index

      block_index = get_block_index(array, block, stack, blkndx, 
     *                              x, .true.)
      if (block_index .ne. -1) then
         get_block_data_index = block_index 
      else
         get_block_data_index = -1
      endif
      return
      end 

      subroutine clear_block(array_handle, array_block, xstack, 
     *                       blkndx, size)
c--------------------------------------------------------------------------
c   Sets the content of a block to 0.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      integer array_handle, array_block, xstack, blkndx, size
      integer j
      integer*8 indblk, get_block_index
      integer stack

#ifdef ALTIX
      dptr = dshptr
#endif
      indblk = get_block_index(array_handle, array_block, stack, blkndx,
     *                         x, .false.)
      if (indblk .eq. -1) then
         print *,'Error: Block not found: array, block = ',
     *        array_handle, array_block
         call abort_job()
      endif

      call validate_block(x(indblk), array_handle, array_block, 
     *                       .true.)
      indblk = indblk + lblock_id_data

      call get_blk_header(stack, blkndx,c_block_stack)
      if (size .gt. stack_blocksize(stack)) then
         print *,'Error: Block size > stack blocksize'
         print *,'array_handle ',array_handle,' blkndx ',blkndx
         print *,'size, stack, stack_blocksize(stack) ',
     *       size, stack, stack_blocksize(stack)
         call dump_block_ids()
         call abort_job()
      endif

      do j = 1, size
         x(indblk+j-1) = 0.
      enddo

      return
      end

      subroutine fill_block(array_handle, array_block, stack,
     *                      blkndx, value)
c--------------------------------------------------------------------------
c   Sets the content of a block to a given value.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      double precision value

      integer array_handle, array_block, stack, blkndx
      integer j, n
      integer*8 indblk, get_block_index
#ifdef ALTIX
      dptr = dshptr
#endif

      indblk = get_block_index(array_handle, array_block, stack, blkndx,
     *                         x, .true.)
      if (indblk .eq. -1) then
         print *,'Error: Block not found: array, block = ',
     *        array_handle, array_block
         call abort_job()
      endif

      call get_blk_header(n, blkndx,c_block_size)
      do j = 1, n
         x(indblk+j-1) = value
      enddo

      return
      end

      subroutine multiply_block(array_handle, array_block, stack,
     *                          blkndx, value,
     *                          result_array, result_block,
     *                          result_stack, result_blkndx, 
     *                          blocksize)
c--------------------------------------------------------------------------
c   Sets the content of a block to 0.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision x(1)
#ifdef ALTIX
      pointer (dptr, x)
#else
      common x
#endif

      double precision value

      integer array_handle, array_block, stack, blkndx
      integer result_array, result_block, result_stack, result_blkndx
      integer j
      integer*8 indblk, iresult, get_block_index
      integer blocksize

#ifdef ALTIX
      dptr = dshptr
#endif

      indblk = get_block_index(array_handle, array_block, stack, blkndx, 
     *                         x, .true.)
      if (indblk .eq. -1) then
         print *,'Error: Block not found: array, block = ',
     *        array_handle, array_block
         call abort_job()
      endif

      iresult = get_block_index(result_array, result_block, 
     *                          result_stack, result_blkndx, 
     *                          x, .true.)
      if (indblk .eq. -1) then
         print *,'Error: Block not found: array, block = ',
     *        result_array, result_block
         call abort_job()
      endif

      do j = 1, blocksize
         x(iresult+j-1) = x(indblk+j-1) * value
      enddo

      return
      end

      subroutine dump_block_ids()
c--------------------------------------------------------------------------
c   Dump the current state of the blocks.                       
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer i, j, k, n, ierr
      integer seg(mx_array_index)
      integer*8 indblk, get_block_index
      integer*8 data_addr, block_addr
      integer size
      integer flag, flagval, persistent
      integer data_block, data_array
      integer blkndx, stack
      integer nused(nblkmgr_stacks)
      integer array_data, blk_data, nind, request_data
      double precision x(2)
#ifdef ALTIX
      pointer (dptr, x)

      dptr = dshptr
#endif
 
      print *,'Dump of busy blocks on processor ',me,
     *     ' max_blkmgr_blocks = ',max_blkmgr_blocks

      size = bytes_per_double
      i = 1
      do stack = 1, nblkmgr_stacks
      n = stack_size(stack)
      nused(stack)  = 0 
      do j = 1, n 
         call get_blk_header(flag, i,c_block_flags)
         flagval = and(flag, block_busy_flag)
         persistent = and(flag, block_persistence_flag)
         if (flagval .ne. 0 .or.
     *       persistent .ne. 0) then
            call get_blk_addr(block_addr,i)
            data_addr  = block_addr + lblock_id_data * size
            nused(stack) = nused(stack) + 1
 
            call get_block_id(stack, i,data_array, data_block)

            call get_blk_header(array_data, i,c_array_handle)
            call get_blk_header(blk_data, i,c_array_block)
            call get_block_segments(i, seg)
            call get_blk_header(nind, i, c_block_nindex)
            print *,'Block ',i,' ARRAY ID: ',array_data,
     *             ' BLOCK ID ',blk_data,
     *             ' FLAGS ',flag,
     *             ' SEGS ',(seg(k),k=1,nind),
     *               ' ADDR ',
     *               block_addr,' DADDR ',data_addr,
     *               'DATA ARRAY ',data_array,' DATA BLOCK ',
     *               data_block


            indblk = get_block_index(array_data, blk_data,
     *                               stack, i, x, .false.)
            call validate_block(x(indblk), array_data, blk_data,
     *                          .false.)
         endif

         i = i + 1
      enddo
      enddo

      print *,'Task ',me,' Stack block usage: ',
     *   (nused(i),i=1,nblkmgr_stacks)

      i = 1
      do stack = 1, nblkmgr_stacks
      n = stack_size(stack)
      do j = 1, n
        call get_blk_header(array_data, i,c_array_handle)
        call get_blk_header(blk_data, i,c_array_block) 
        call get_blk_header(flag, i,c_block_flags) 
        call get_blk_header(request_data, i, c_block_request)
        call get_blk_header(size, i,c_block_size)
        call get_blk_header(stack, i,c_block_stack)
        print *,' blkndx ',i,' array,block,flags ',
     *     array_data, blk_data, flag, 
     *     ' request ',request_data,
     *     ' size ', size, 
     *     ' stack ',stack

         i = i + 1
      enddo
      enddo
      return
      end

      subroutine validate_all_blocks()
c--------------------------------------------------------------------------
c   Validates the current state of the blocks.                       
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer i, j, k, n, ierr
      integer*8 indblk, get_block_index
      integer*8 data_addr, block_addr
      integer flag, array_data, blk_data, size
      integer flagval, persistent
      integer data_block, data_array
      integer blkndx, stack
      double precision x(2)
#ifdef ALTIX
      pointer (dptr, x)

      dptr = dshptr
#endif
 
      size = bytes_per_double
      i = 1
      do stack = 1, nblkmgr_stacks
      n = stack_size(stack)
      do j = 1, n 
         call get_blk_header(flag, i,c_block_flags)
         flagval = and(flag,  block_busy_flag)
         persistent = and(flag,
     *                    block_persistence_flag)
         if (flagval .ne. 0 .or.
     *       persistent .ne. 0) then
 
            blkndx = 0
            call get_blk_header(array_data, i,c_array_handle)
            call get_blk_header(blk_data, i,c_array_handle)
            indblk = get_block_index(array_data, blk_data,
     *                               stack, blkndx, x, .false.)
            call validate_block(x(indblk), array_data, blk_data,
     *                          .true.)
         endif

         i = i + 1
      enddo
      enddo

      return
      end

      subroutine dump_block_data(array, block)
c--------------------------------------------------------------------------
c   Dump the current state of the blocks.                       
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer i, j, ierr
      integer*8 indblk, get_block_index
      integer blkndx 
      integer*8 data_addr, block_addr
      integer size
      integer flagval
      integer data_block, data_array
      integer array, block
      integer array_data, blk_data, flag, nind
      integer seg(mx_array_index)
      integer get_block_number
      integer stack
      double precision x(2)
#ifdef ALTIX
      pointer (dptr, x)

      dptr = dshptr
#endif
 
      i = get_block_number(array, block)
      print *,'Task ',me,' DUMP_BLOCK_DATA: array, block, i = ',
     *   array, block, i
      if (i .lt. 0) return

      do j = 1, nblkmgr_stacks
         if (i .ge. stack_start(j) .and.
     *       i .le. stack_start(j)+stack_size(j) - 1) then
            stack = j
            go to 100
         endif
      enddo

      size = bytes_per_double
  100 continue
      call get_blk_addr(block_addr, i)
      data_addr  = block_addr + lblock_id_data * size
 
      call get_block_id(stack, i,data_array, data_block)

      call get_blk_header(array_data, i,c_array_handle)
      call get_blk_header(blk_data, i,c_array_block)
      call get_blk_header(flag, i,c_block_flags)
      call get_block_segments(i, seg)
      call get_blk_header(nind, i,c_block_nindex) 

      print *,'Task ',me,' Block ',i,' ARRAY ID: ',array_data,
     *             ' BLOCK ID ',blk_data,
     *             ' FLAGS ',flag,
     *             ' SEGS ',(seg(j),j=1,nind),
     *             ' ADDR ',  block_addr,' DADDR ',data_addr,
     *               'DATA ARRAY ',data_array,' DATA BLOCK ',
     *               data_block

      blkndx = 0
      indblk = get_block_index(array_data, blk_data,
     *                         stack, blkndx, x, .false.)

      call validate_block(x(indblk), array_data, blk_data,
     *                    .false.)
      return
      end

      integer function find_current_block(x, array_table,
     *                             index_table, nindex_table,
     *                             segment_table, nsegment_table,
     *                             block_map_table, blkndx)
c--------------------------------------------------------------------------
c   Determines whether the block with the current ranges of indices is 
c   present in the set of blocks being managed.  If the block is found,
c   the block index is returned as a negative index.  If the block is 
c   in the array mapping AND also in the table of managed blocks, the 
c   block index is returned.
c
c   x = array handle 
c   array_table = table of arrays entry for array x.
c---------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer x, array_table(larray_table_entry)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer nsegment_table
      integer segment_table(lsegment_table_entry, nsegment_table)
      integer i, j, nindex
      integer index(mx_array_index), seg(mx_array_index)
      integer block_map_table(lblock_map_entry,*)
      integer array, array_block, blkndx
      integer array_data, data
      integer flag, flagval, flagtest
      integer type
      integer next
      integer lookup, block_map_lookup
      integer stack
      integer seg_blk
      integer ierr
      integer orig_index
      logical match
      integer*8 ixx, get_index_from_base
      integer bhdr(1) 

c---------------------------------------------------------------------------
c   Get the current segment of each index.
c---------------------------------------------------------------------------

      find_current_block = -1
      blkndx = 0

      type = array_table(c_array_type)
      if (type .eq. scalar_value .or.
     *    type .eq. static_array) then
         return
      endif

      nindex = array_table(c_nindex)
      do i = 1, nindex
         index(i) = array_table(c_index_array1+i-1)
         if (index(i) .gt. nindex_table) then
            return
         endif

         seg(i)   = index_table(c_current_seg,index(i))

c---------------------------------------------------------------------------
c   If this is a subindex, do we look up the block defined by the subindex
c   itself, or the block defined by its superindex?
c---------------------------------------------------------------------------

         if (index_table(c_index_type, index(i)) .eq. subindex) then
            orig_index = array_table(c_index_original+i-1)
            if (orig_index .gt. 0) then
               if (index_table(c_index_type,orig_index) .ne. 
     *                                             subindex) then
         
c---------------------------------------------------------------------------
c   Find the segment of the superindex corresponding to the current block
c   of the subindex.
c---------------------------------------------------------------------------

                  call subindex_seg_match(index(i), index_table,
     *                                    nindex_table, segment_table,
     *                                    nsegment_table, seg(i))
                  index(i) = orig_index 
               endif
            endif
         endif

         if (seg(i) .eq. undefined_segment) then
            return
         endif
      enddo

      stack = array_table(c_array_stack)

      if (type .eq. distributed_array) then

c---------------------------------------------------------------------------
c   Perform a direct lookup in the block_map_table.
c---------------------------------------------------------------------------

         lookup = block_map_lookup(seg, nindex, x, array_table,
     *                         index_table, nindex_table)
         blkndx = block_map_table(c_bmap_blkndx, lookup)
         if (blkndx .eq. 0) then
            find_current_block = -1
         else
            call get_block_id(stack, blkndx, array, array_block)

            find_current_block = array_block
         endif

         return
      endif

c-------------------------------------------------------------------------
c   Search for a matching block in the most recent blkndx.
c-------------------------------------------------------------------------

       blkndx = array_table(c_current_blkndx)
       if (blkndx .ne. 0) then

c---------------------------------------------------------------------------
c   If the array handle stored in the block itself doesn't match the 
c   expected array number, the blkndx has become invalid and we must
c   search the entire block cache.
c----------------------------------------------------------------------------

         call get_block_id(stack, blkndx, array, array_block)
         if (array .ne. x) go to 1000

c---------------------------------------------------------------------------
c   Either the block must be "busy" or the "persistent" bit must be set. 
c   Otherwise, we must search the entire block cache.
c---------------------------------------------------------------------------

         call get_blk_header(flag, blkndx,c_block_flags)
         if (and(flag, block_busy_flag) .eq. 0) then
            if (and(flag, block_persistence_flag) .eq. 0) then
               array_table(c_current_blkndx) = 0
               go to 1000
            endif
         endif

         do j = 1, nindex
            call get_blk_header(seg_blk, blkndx,c_block_segments+j-1)
            if (seg(j) .ne. seg_blk) then
               array_table(c_current_blkndx) = 0
               blkndx = 0
               go to 1000
            endif
         enddo

         find_current_block = array_block
         return
       endif

c-------------------------------------------------------------------------
c   Search for a matching block in the blk_header table.
c-------------------------------------------------------------------------

 1000 continue

c-------------------------------------------------------------------------
c   Search the linked list for the blocks of this array.
c-------------------------------------------------------------------------

      next = array_table(c_block_list)
      flagtest = or(block_busy_flag,block_persistence_flag)
      ixx = get_index_from_base(addr_blk_header, bhdr, 1)

      if (next .ne. 0) then
 2000 continue
            i = next
            array_data = 
     *            bhdr(ixx+(c_array_handle-1)*blkmgr_blocks+i-1)
c            call get_blk_header(array_data, i, c_array_handle)
            if (array_data .eq. x) then 

c-------------------------------------------------------------------------
c   Make certain that this block has been busied.
c-------------------------------------------------------------------------

               flag = 
     *            bhdr(ixx+(c_block_flags-1)*blkmgr_blocks+i-1) 
c               call get_blk_header(flag, i,c_block_flags)
               if (and(flag, flagtest) .eq. 0) then
                   if (and(flag, block_scrub_flag) .ne. 0) then

c--------------------------------------------------------------------------
c   This block is neither busy nor persistent.  It is marked for scrubbing.
c   Attempt to scrub the block now.
c--------------------------------------------------------------------------
 
                      call scrub_this_block(i,
     *                        index_table, nindex_table,
     *                        block_map_table, ierr)
                    endif
                   go to 1500
               endif

               do j = 1, nindex
c                  call get_blk_header(data, i,c_block_segments+j-1)
                  data = 
     *               bhdr(ixx+(c_block_segments+j-2)*blkmgr_blocks+i-1)
                  if (seg(j) .ne. data)  go to 1500
               enddo

               call get_block_id(stack, i, array, array_block)

               find_current_block = array_block
               blkndx = i
               array_table(c_current_blkndx) = i   ! Update array table entry
               return
            endif

 1500       continue
            call get_blk_header(next, i,c_block_list_ptr)
            if (next .eq. i) then
               print *,'Error: Recursion in block list'
               print *,'x = ',x,' head = ',array_table(c_block_list)
               call abort_job()
            endif

            if (next .ne. 0) go to 2000
      endif     ! next .ne. 0

      find_current_block = -1
      blkndx = 0
      array_table(c_current_blkndx) = 0
      return
      end

      integer function get_blkmgr_blocksize(stack)
c--------------------------------------------------------------------------
c   Returns the size of the blocksize used for the data area of each blkmgr
c   block, in double precision words.
c   If stack = 0, the function returns the maximum blocksize of all stacks.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer stack
      integer i, mxblk
 
      if (stack .eq. 0) then
         mxblk = 0
         do i = 1, nblkmgr_stacks
            mxblk = max(mxblk, stack_datasize(i))
         enddo
         get_blkmgr_blocksize = mxblk
      else
         get_blkmgr_blocksize = stack_datasize(stack)
      endif
      return
      end

      integer function get_blkmgr_numblocks(stack)
c--------------------------------------------------------------------------
c   Returns the number of blocks managed by the blkmgr routines.
c   If stack = 0, we return the total number of blocks over all the stacks.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer stack
      integer i, total_blocks
 
      if (stack .eq. 0) then
         total_blocks = 0
         do i = 1, nblkmgr_stacks
            total_blocks = total_blocks + stack_size(i)
         enddo 
         get_blkmgr_numblocks = total_blocks
      else
         get_blkmgr_numblocks = stack_size(stack)
      endif
      return
      end

      subroutine set_block_indices(array, block, blkndx, 
     *                             array_table_entry)
c-----------------------------------------------------------------------------
c   Sets the block indices in the block header for a given array and block.
c-----------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'

      integer array, block, blkndx
      integer array_table_entry(larray_table_entry)
      integer i, j, n

      i = blkndx
      n = array_table_entry(c_nindex)
      call put_blk_header(n, i,c_block_nindex)
      do j = 1, n
         call put_blk_header(array_table_entry(c_index_array1+j-1),
     *             i,c_block_indices+j-1)
      enddo

      return
      end

      subroutine set_block_segments(array, block, blkndx, 
     *                              index_table, nindex_table)
c--------------------------------------------------------------------------
c   Sets the block segments to the current values of the segment for each 
c   index in the index_table.
c--------------------------------------------------------------------------
      
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'

      integer array, block, blkndx
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer i, j, n, indx

      i = blkndx
      call get_blk_header(n, i,c_block_nindex)

      do j = 1, n
         call get_blk_header(indx, i,c_block_indices+j-1)
         if (index_table(c_current_seg,indx) .le. 0) then
            print *,'Error: current segment of index ',indx,
     *                    ' is ',index_table(c_current_seg,indx)
            print *,'set_block_segments was setting array ',
     *                  array,' block ',block,' at time of error.'
            call abort_job()
         endif

         call put_blk_header(index_table(c_current_seg,indx), 
     *               i,c_block_segments+j-1)

      enddo

      return

      print *,'Error: block ',block,' of array ',array,
     *        ' was not found in set_block_segments.'
      call dump_block_ids()
      call abort_job()
      return
      end

      subroutine set_block_busy_flag(array, array_block, blkndx,
     *                               value)
c-------------------------------------------------------------------------
c   Sets or clears the block_busy_flag for a given array, block 
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx

      if (array_block .lt. 1) return

      if (value .eq. 1) then
         call set_block_flag(array, array_block, blkndx, 
     *                       block_busy_flag)
      else if (value .eq. 0) then
         call clear_block_flag(array, array_block, blkndx,
     *                         block_busy_flag) 
      endif
      return
      end

      subroutine set_block_persistence_flag(array, array_block,
     *                               blkndx, value)
c-------------------------------------------------------------------------
c   Sets or clears the block_persistence_flag for a given array, block 
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx

      if (array_block .lt. 1) return

      if (value .eq. 1) then
         call set_block_flag(array, array_block, blkndx, 
     *                       block_persistence_flag)
      else if (value .eq. 0) then
         call clear_block_persistence_flag(array, array_block, 
     *                         blkndx) 
      endif
      return
      end

      subroutine set_block_scrub_flag(array, array_block, blkndx,
     *                                value)
c-------------------------------------------------------------------------
c   Sets or clears the block_scrub_flag for a given array, block 
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx

      if (array_block .lt. 1) return

      if (value .eq. 1) then
         call set_block_flag(array, array_block, blkndx,
     *                       block_scrub_flag)
      else if (value .eq. 0) then
         call clear_block_flag(array, array_block, blkndx,
     *                          block_scrub_flag) 
      endif
      return
      end

      subroutine set_block_computed_flag(array, array_block, blkndx,
     *                                   value)
c-------------------------------------------------------------------------
c   Sets or clears the block_computed_flag for a given array, block 
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx

      if (array_block .lt. 1) return

      if (value .eq. 1) then
         call set_block_flag(array, array_block, blkndx,
     *                       block_computed_flag)
      else if (value .eq. 0) then
         call clear_block_flag(array, array_block, blkndx,
     *                         block_computed_flag) 
      endif
      return
      end

      subroutine set_block_created_flag(array, array_block, blkndx,
     *                                  value)
c-------------------------------------------------------------------------
c   Sets or clears the block_created_flag for a given array, block 
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx

      if (array_block .lt. 1) return

      if (value .eq. 1) then
         call set_block_flag(array, array_block, blkndx,
     *                       block_created_flag)
      else if (value .eq. 0) then
         call clear_block_flag(array, array_block, blkndx,
     *                         block_created_flag) 
      endif
      return
      end

      subroutine set_block_request_outstanding_flag(array, 
     *                                    array_block, blkndx)
c-------------------------------------------------------------------------
c   Sets the block_request_outstanding_flag for a given array, block 
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, blkndx

      if (array_block .lt. 1) return

      call set_block_flag(array, array_block, blkndx,
     *                    block_request_outstanding_flag)
      return
      end

      subroutine set_block_request_flag(array, array_block,
     *                               blkndx, value)
c-------------------------------------------------------------------------
c   Sets or clears the block_request_flag for a given array, block
c   combination.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      integer array, array_block, value, blkndx

      if (array_block .lt. 1) return

      if (value .eq. 1) then
         call set_block_flag(array, array_block, blkndx,
     *                       block_request_flag)
      else if (value .eq. 0) then
         call clear_block_request_flag(array, array_block,
     *                         blkndx)
      endif
      return
      end

      subroutine get_block_computed_flag(array, array_block, blkndx,
     *                                   value)
c--------------------------------------------------------------------------
c   The value of the block_computed_flag for a given array, block setting
c   is returned in the "value" argument.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx
      integer iblk, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id #1'
         call dump_block_ids()
         call abort_job()
      endif

      call get_blk_header(flagval, iblk,c_block_flags)
      value = and(flagval, block_computed_flag)

      return
      end

      subroutine get_block_scrub_flag(array, array_block, blkndx,
     *                                value)
c--------------------------------------------------------------------------
c   The value of the block_scrub_flag for a given array, block setting
c   is returned in the "value" argument.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx
      integer iblk, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id #2'
         call dump_block_ids()
         call abort_job()
      endif

      call get_blk_header(flagval, iblk,c_block_flags)
      value = and(flagval, block_scrub_flag)

      return
      end

      subroutine get_block_created_flag(array, array_block, blkndx,
     *                                  value)
c--------------------------------------------------------------------------
c   The value of the block_created_flag for a given array, block setting
c   is returned in the "value" argument.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx
      integer iblk, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 3'
         call dump_block_ids()
         call abort_job()
      endif

      call get_blk_header(flagval, iblk,c_block_flags)
      value = and(flagval, block_created_flag)

      return
      end

      subroutine get_block_persistence_flag(array, array_block,
     *                                blkndx, value)
c--------------------------------------------------------------------------
c   The value of the block_persistence_flag for a given array, block setting
c   is returned in the "value" argument. 
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      include 'trace.h'
      integer array, array_block, value, blkndx
      integer iblk, flagval
      integer array_data, blk_data

      iblk = blkndx
      call get_blk_header(array_data, iblk,c_array_handle)
      call get_blk_header(blk_data, iblk,c_array_block)
      if (array_data .ne. array .or.
     *    blk_data .ne. array_block) then
         print *,'Error: Inconsistent data: array, block, blkndx ',
     *             array, array_block, blkndx 
         print *,'Task ',me,' line ',current_line,' actual block was ',
     *      array_data, blk_data
         call dump_block_ids()
         call abort_job()
      endif

      call get_blk_header(flagval, iblk,c_block_flags)
      value = and(flagval, block_persistence_flag)

      return
      end

      subroutine check_block_consistency(array, array_block,
     *                                blkndx)
c--------------------------------------------------------------------------
c   Check to make sure the argument values match the values stored in the
c   data block.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      include 'trace.h'
      integer array, array_block, value, blkndx
      integer iblk, flagval
      integer array_data, blk_data

      iblk = blkndx
      call get_blk_header(array_data, iblk,c_array_handle)
      call get_blk_header(blk_data, iblk,c_array_block)
      if (array_data .ne. array .or.
     *    blk_data .ne. array_block) then
         print *,'Error: Inconsistent data: array, block, blkndx ',
     *             array, array_block, blkndx 
         print *,'Task ',me,' line ',current_line,' actual block was ',
     *      array_data, blk_data
         call dump_block_ids()
         call abort_job()
      endif

      return
      end

      subroutine get_block_request_outstanding_flag(array, 
     *                               array_block, blkndx, value)
c--------------------------------------------------------------------------
c   The value of the block_request_outstanding_flag for a given array, 
c   block setting is returned in the "value" argument.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, value, blkndx
      integer iblk, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 5'
         call dump_block_ids() 
         call abort_job()
      endif

      call get_blk_header(flagval, iblk,c_block_flags)
      value = and(flagval, 
     *             block_request_outstanding_flag)

      return
      end

      subroutine clear_block_computed_flag(array, array_block, 
     *                                     blkndx)
c--------------------------------------------------------------------------
c   Clears all the block_computed_flag values in the table of blocks.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'

      integer i,j,flag, flagval
      integer array, array_block, blkndx

      i = blkndx
      if (i .lt. 1) then
         print *,'Error: clear_block_computed_flag cannot find',
     *     ' array, block ',array, array_block
         call dump_block_ids()
         call abort_job()
      endif

      call get_blk_header(flagval, i,c_block_flags)
      flag = and(flagval, block_computed_flag)
      if (flag .ne. 0) 
     *   call put_blk_header(xor(block_computed_flag,flagval),
     *                       i,c_block_flags)    

      return
      end

      subroutine clear_block_created_flag(array, array_block, blkndx)
c--------------------------------------------------------------------------
c   Clears all the block_created_flag values in the table of blocks.
c--------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'

      integer i,flag, flagval
      integer array, array_block, blkndx

      i = blkndx
      if (i .lt. 1) then
         print *,'Error: clear_block_created_flag cannot find',
     *     ' array, block ',array, array_block
         call abort_job()
      endif

      call get_blk_header(flagval, i,c_block_flags)
      flag = and(flagval, block_created_flag)
      if (flag .ne. 0) 
     *   call put_blk_header(xor(block_created_flag,flagval),
     *                       i,c_block_flags)

      return
      end

      subroutine set_block_flag(array, array_block, blkndx, 
     *                           block_flag)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, block_flag, blkndx
      integer iblk, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 6'
         call abort_job()
      endif

      call get_blk_header(flagval, iblk,c_block_flags)
      flagval = or(flagval, block_flag)
      call put_blk_header(flagval, iblk,c_block_flags)
      return
      end

      subroutine clear_block_flag(array, array_block, blkndx,
     *                            block_flag)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, block_flag, blkndx
      integer iblk, flag, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 7'
         call abort_job()
      endif

      call get_blk_header(flag, iblk,c_block_flags)  
      flagval = and(flag, block_flag)
      call put_blk_header(xor(flagval, flag), iblk,c_block_flags) 
      return
      end

      subroutine clear_block_persistence_flag(array, array_block, 
     *                                        blkndx)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, blkndx
      integer iblk, flag, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 8'
         call abort_job()
      endif

      call get_blk_header(flag, iblk,c_block_flags)
      flagval = and(flag, 
     *                  block_persistence_flag)
      call put_blk_header(xor(flagval, flag), iblk,c_block_flags)

      return
      end

      subroutine clear_block_request_flag(array, array_block, 
     *                                        blkndx)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, blkndx
      integer iblk, flag, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 8'
         call abort_job()
      endif

      call get_blk_header(flag, iblk,c_block_flags)
      flagval = and(flag, 
     *                  block_request_flag)
      call put_blk_header(xor(flagval, flag), iblk,c_block_flags)

      return
      end

      subroutine clear_block_request_outstanding_flag(array, 
     *                                     array_block, blkndx)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      integer array, array_block, blkndx 
      integer iblk, flag, flagval

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: Block not found by blkmgr: array ',array,
     *           ' array_block = ',array_block,' id # 9'
         call abort_job()
      endif

      call get_blk_header(flag, iblk,c_block_flags)
      flagval = and(flag, 
     *              block_request_outstanding_flag)
      if (flagval .ne. 0)
     *   call put_blk_header(xor(flagval, flag), iblk,c_block_flags)
      return
      end

      subroutine set_block_size(array, array_block, blkndx, size)
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer array, array_block, size, blkndx
      integer iblk

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: set_block_size is unable to find array ',
     *           array,' block ',array_block
         call abort_job()
      endif

      call put_blk_header(size, iblk,c_block_size)
      return
      end 

      integer function get_block_size(array, array_block, blkndx)
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer array, array_block, blkndx
      integer iblk, n

      iblk = blkndx
      if (iblk .le. 0) then
         print *,'Error: get_block_size is unable to find array ',
     *           array,' block ',array_block
         call abort_job()
      endif

      call get_blk_header(n, iblk,c_block_size)
      get_block_size = n
      return
      end

      integer function get_block_number(array, array_block)
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'machine_types.h'

      integer i, j, array, array_block
      integer flag, array_data, blk_data
      integer blk_in_use(1)
      integer*8 addr, get_index_from_base
      integer*8 ixx

      ixx = get_index_from_base(addr_blk_in_use, blk_in_use, 1)
      
c--------------------------------------------------------------------------
c   Scan all blocks that are registered as "in use".
c--------------------------------------------------------------------------

      do  j = 1, nblk_in_use
         i = blk_in_use(ixx+j-1)
 
         call get_blk_header(flag, i,c_block_flags)
         if (and(flag, block_busy_flag) 
     *              .ne. 0 .or.
     *       and(flag,
     *                 block_persistence_flag) .ne. 0) then
            call get_blk_header(array_data, i,c_array_handle)
            call get_blk_header(blk_data, i,c_array_block)
            if (array_data .eq. array .and.
     *             blk_data .eq. array_block) then
               get_block_number = i
               return
            endif
         endif
      enddo

      get_block_number = -1 
      return
      end

      integer function get_block_number_distributed(array, array_block,
     *                     array_table, block_map_table)
c----------------------------------------------------------------------------
c   Fast lookup of a block number of a distributed array.
c   No checking is done to insure correctness, and the array is assumed
c   to have been created previously.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'machine_types.h'
      include 'mpif.h'
      integer me, ierr

      integer i, j, array, array_block
      integer array_table(larray_table_entry,*)
      integer block_map_table(lblock_map_entry,*)
 
      integer flag, array_data, blk_data
      
      j = array_table(c_block_map,array)   ! block_map entry of 1st block.
      get_block_number_distributed = 
     *    block_map_table(c_bmap_blkndx, j+array_block-1)
      
      return
      end

      subroutine set_block_request(array, block, blkndx, request)
c-------------------------------------------------------------------------
c   Stores the mpi request handle in the blkmgr data for the array, block.
c-------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'mpif.h'
      integer array, block, blkndx, request
      integer j
      integer request_data

      if (array .le. 0 .or. block .le. 0) then
         print *,'Error: set_block_request was called with ',
     *      ' array ',array,' block ',block
         call abort_job()
      endif

      j = blkndx
      if (j .gt. 0) then
         call get_blk_header(request_data, j,c_block_request)
         if (request_data .ne. mpi_request_null .and.
     *      request .ne. mpi_request_null) then
            print *,'Error in set_block_request: Previous request',
     *          ' has not been cleared. array, block = ',
     *          array,block
            print *,'Previous request = ',request_data,
     *               ' mpi_request_null = ',mpi_request_null
            call abort_job()
         endif

         call put_blk_header(request, j,c_block_request)
      else
         print *,'Error: set_block_request cannot find array ',
     *     array,' block ',block
         call dump_block_ids()
         call abort_job()
      endif 

      return
      end

      integer function get_block_request(array, block, blkndx)
c-------------------------------------------------------------------------
c   Returns the value of the request stored in the blkmgr's data for the
c   array, block.
c-------------------------------------------------------------------------

      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer array, block, blkndx
      integer j, request

      j = blkndx
      if (j .gt. 0) then
         call get_blk_header(request, j,c_block_request)
         get_block_request = request
      else
         print *,'Error: get_block_request cannot find array ',
     *     array,' block ',block
         call dump_block_ids()
         call abort_job()
      endif 

      return
      end

      subroutine scrub_blocks(size, must_wait, array_table,
     *                        narray_table, index_table, nindex_table, 
     *                        block_map_table, ierr)
c----------------------------------------------------------------------------
c   Attempt to scrub a block larger than a given size.
c   Searches all blocks marked as "scrubbable".
c   If a block has been marked to be scrubbed, and if the block is on a stack
c   with a large enough stack size, the communication request
c   is checked.  If the request has been completed, the block_busy_flag
c   is cleared, indicating that the block may be reclaimed for later use.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_monitor.h'

      logical must_wait
      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer ierr, mpi_ierr
      integer size

      integer array, block, stack
      integer flag, comp_flag
      integer flagval
      integer request
      integer next, ptr
      integer status(MPI_STATUS_SIZE)
      logical complete

      ierr = 1
      ptr = comm_list_head
      if (ptr .eq. 0) return

  100 continue
      call get_blk_header(flagval, ptr,c_block_flags)
      flag = and(flagval, block_scrub_flag)
      call get_blk_header(next, ptr,c_comm_list_ptr)

      if (flag .ne. 0) then

c---------------------------------------------------------------------------
c   Check the size requirement.
c---------------------------------------------------------------------------

         call get_blk_header(stack, ptr,c_block_stack)
         if (size .le. stack_datasize(stack)) then
            call get_blk_header(request, ptr,c_block_request)
            call mpi_test(request, complete,
     *                    status, mpi_ierr)
            if (complete) then
               call put_blk_header(MPI_REQUEST_NULL,ptr,c_block_request)
               call get_blk_header(array, ptr, c_array_handle)
               call get_blk_header(block, ptr, c_array_block)
               call blkmgr_remove_block_from_list(comm_list_head,
     *             comm_list_tail, ptr, c_comm_list_ptr)

               if (server_monitor_on) 
     *            call server_monitor_write_log(ptr, 'r', 
     *                            status(mpi_source))

               call free_block(array, block, ptr, array_table,
     *                narray_table, index_table, nindex_table,
     *                block_map_table)
               ierr = 0
               return
            endif      ! complete
	 endif      ! size
      endif      ! flag

      if (next .ne. 0) then
         ptr = next
         go to 100
      endif 

c      print *,'Task ',me,'   SCRUB: Returning ierr ',ierr
      return
      end

      subroutine scrub_this_block(blkndx, 
     *                        index_table, nindex_table, 
     *                        block_map_table, ierr)
c----------------------------------------------------------------------------
c   Attempt to scrub a particular block.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_monitor.h'
      include 'sip_tables.h'

      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer ierr, mpi_ierr
      integer blkndx

      integer*8 ixarray_table, get_index_from_base
      integer array_table(1)
      integer narray_table

      integer array, block, stack
      integer request
      integer status(MPI_STATUS_SIZE)
      logical complete

      ierr = 1

      call get_blk_header(request, blkndx,c_block_request)
      call mpi_test(request, complete,
     *                    status, mpi_ierr)
      if (complete) then
         call put_blk_header(MPI_REQUEST_NULL,blkndx,c_block_request)
         call get_blk_header(array, blkndx, c_array_handle)
         call get_blk_header(block, blkndx, c_array_block)
         call blkmgr_remove_block_from_list(comm_list_head,
     *             comm_list_tail, blkndx, c_comm_list_ptr)

         if (server_monitor_on) 
     *      call server_monitor_write_log(blkndx, 'r', 
     *                      status(mpi_source))

         ixarray_table = get_index_from_base(array_table_base,
     *                       array_table, 1)
         call free_block(array, block, blkndx, 
     *          array_table(ixarray_table),
     *          narray_table_sip, index_table, nindex_table,
     *          block_map_table)
         ierr = 0
      endif

      return
      end

      subroutine scrub_from_stack(scrub_stack, array_table,
     *                        narray_table, index_table, nindex_table, 
     *                        block_map_table, ierr)
c----------------------------------------------------------------------------
c   Attempt to scrub a block from a given stack.
c   Searches all blocks marked as "scrubbable".
c   If a block has been marked to be scrubbed, and if the block is on the 
c   specified stack, the communication request
c   is checked.  If the request has been completed, the block_busy_flag
c   is cleared, indicating that the block may be reclaimed for later use.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_monitor.h'

      integer scrub_stack
      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer ierr, mpi_ierr

      integer array, block, stack
      integer flag, comp_flag
      integer flagval
      integer request
      integer next, ptr
      integer status(MPI_STATUS_SIZE)
      logical complete

      ierr = 1
      ptr = comm_list_head
      if (ptr .eq. 0) return

  100 continue
      call get_blk_header(flagval, ptr,c_block_flags)
      flag = and(flagval,  block_scrub_flag)
      call get_blk_header(next, ptr,c_comm_list_ptr)

      if (flag .ne. 0) then

c---------------------------------------------------------------------------
c   Check the size requirement.
c---------------------------------------------------------------------------

         call get_blk_header(stack, ptr,c_block_stack)
         if (stack .eq. scrub_stack) then
            call get_blk_header(request, ptr,c_block_request) 
            call mpi_test(request, complete,
     *                    status, mpi_ierr)
            if (complete) then
               call put_blk_header(MPI_REQUEST_NULL, 
     *                            ptr,c_block_request)
               call get_blk_header(array, ptr, c_array_handle)
               call get_blk_header(block, ptr, c_array_block) 
               call blkmgr_remove_block_from_list(comm_list_head,
     *             comm_list_tail, ptr, c_comm_list_ptr)

               if (server_monitor_on)
     *            call server_monitor_write_log(ptr, 'r', 
     *                               status(mpi_source))

               call free_block(array, block, ptr, array_table,
     *                narray_table, index_table, nindex_table,
     *                block_map_table)
               ierr = 0
               return
            endif
         endif
      endif

      if (next .ne. 0) then
         ptr = next
         go to 100
      endif 

c      print *,'Task ',me,'   SCRUB: Returning ierr ',ierr
      return
      end

      subroutine reclaim_persistent_block(size, array_table,
     *                     narray_table, index_table, nindex_table,
     *                     block_map_table, ierr)
c---------------------------------------------------------------------------
c   Find a block marked as "persistent" exceeding a given "size", and 
c   frees the block.  This makes the block available for further use.
c
c   If no block is available to be reclaimed, the ierr argument is set to
c   -1.  Otherwise, ierr = 0 is returned.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_monitor.h'

      integer narray_table, ierr, mpi_ierr
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer size, request
      integer i, old_array, block, stack, ptr, next
      integer flagval
      integer status(MPI_STATUS_SIZE)
      logical flag
 
c-------------------------------------------------------------------------
c   Find a persistent block that is not busy.
c-------------------------------------------------------------------------

      ierr = -1
      if (persistent_list_head .eq. 0) return
      
      next = persistent_list_head
  100 continue
      ptr = next
      i   = next 
      call get_blk_header(next, ptr, c_persistent_ptr)
      call get_blk_header(stack, i,c_block_stack)
      if (size .le. stack_datasize(stack)) then

c--------------------------------------------------------------------------
c   We must check the block_busy_flag and block_computed_flag.  If either 
c   one is non-zero, we cannot reclaim the block.
c--------------------------------------------------------------------------

      call get_blk_header(flagval, i,c_block_flags)
      if (and(flagval, block_busy_flag) .eq. 0 .and.
     *    and(flagval, block_computed_flag) .eq. 0) then 

c--------------------------------------------------------------------------
c   This is the oldest such block in the list, since the list is a FIFO.
c   Make sure the block is not still engaged in communication.
c--------------------------------------------------------------------------

         call get_blk_header(request, i,c_block_request)
         if (request .ne. 
     *                                mpi_request_null) then

c--------------------------------------------------------------------------
c   Test the block to see if it has finished its communication.
c   If it has finished, remove it from the comm_list.
c--------------------------------------------------------------------------

            call mpi_test(request,
     *                         flag, status, mpi_ierr)

            if (flag) then 
               call put_blk_header(MPI_REQUEST_NULL, i,c_block_request)
               call blkmgr_remove_block_from_list(comm_list_head,
     *                 comm_list_tail, i, c_comm_list_ptr)
               if (server_monitor_on)
     *            call server_monitor_write_log(ptr, 'r',
     *                                  status(mpi_source))
            endif
         else

c--------------------------------------------------------------------------
c   The block is not found to be engaged in communication.
c   Thus it may be safely reclaimed.
c--------------------------------------------------------------------------

            flag = .true.
         endif

c-------------------------------------------------------------------------
c   Free the block. The free_block subroutine will also remove the 
c   block from the persistent_block list.
c-------------------------------------------------------------------------

         if (flag) then
            call get_blk_header(old_array, i,c_array_handle)
            call get_blk_header(block,i,c_array_block)
            call free_block(old_array, block,
     *                i, array_table, narray_table,
     *                index_table, nindex_table, block_map_table)
            ierr = 0
            return
         endif
      endif    ! flag
      endif    ! size <= stack_size(stack)

c---------------------------------------------------------------------------
c   Point to the next block in the list and check for a block to reclaim.
c---------------------------------------------------------------------------

      if (next .ne. 0) go to 100

      return
      end

      subroutine reclaim_persistent_block_from_stack(mystack, 
     *           array_table,narray_table, index_table, nindex_table,
     *           block_map_table, ierr)
c---------------------------------------------------------------------------
c   Finds a block marked as "persistent" on a given stack, and 
c   frees the block.  This makes the block available for further use.
c
c   If no block is available to be reclaimed, the ierr argument is set to
c   -1.  Otherwise, ierr = 0 is returned.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'parallel_info.h'
      include 'server_monitor.h'

      integer narray_table, ierr, mpi_ierr
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer mystack, flagval, request
      integer i, old_array, block, stack, ptr, next
      integer status(MPI_STATUS_SIZE)
      logical flag
      
c-------------------------------------------------------------------------
c   Find a persistent block that is not busy.
c-------------------------------------------------------------------------

      ierr = -1
      if (persistent_list_head .eq. 0) return
      
      next = persistent_list_head
  100 continue
      ptr = next
      i   = next 
      call get_blk_header(next, ptr, c_persistent_ptr)
      call get_blk_header(stack, i,c_block_stack)
      if (stack .ne. mystack) then
         if (next .gt. 0) go to 100   ! check next entry
         return   ! end of list
      endif

c--------------------------------------------------------------------------
c   We must check the block_busy_flag and block_computed_flag.  If either 
c   one is non-zero, we cannot reclaim the block.
c--------------------------------------------------------------------------

      call get_blk_header(flagval, i,c_block_flags) 
      if (and(flagval, block_busy_flag) 
     *                                     .eq. 0 .and.
     *    and(flagval, block_computed_flag) 
     *                                     .eq. 0) then 

c--------------------------------------------------------------------------
c   This is the oldest such block in the list, since the list is a FIFO.
c   Make sure the block is not still engaged in communication.
c--------------------------------------------------------------------------

         call get_blk_header(request, i,c_block_request)
         if (request .ne. mpi_request_null) then

c--------------------------------------------------------------------------
c   Test the block to see if it has finished its communication.
c   If it has finished, remove it from the comm_list.
c--------------------------------------------------------------------------

            call mpi_test(request,
     *                         flag, status, mpi_ierr)
            if (flag) then
               call put_blk_header(MPI_REQUEST_NULL, i,c_block_request)
               call blkmgr_remove_block_from_list(comm_list_head,
     *                 comm_list_tail, i, c_comm_list_ptr)
               if (server_monitor_on)
     *            call server_monitor_write_log(ptr, 'r',
     *                                status(mpi_source))
            endif
         else

c--------------------------------------------------------------------------
c   The block is not found to be engaged in communication.
c   Thus it may be safely reclaimed.
c--------------------------------------------------------------------------

            flag = .true.
         endif

c-------------------------------------------------------------------------
c   Free the block. The free_block subroutine will also remove the 
c   block from the persistent_block list.
c-------------------------------------------------------------------------

         if (flag) then
            call get_blk_header(old_array, i,c_array_handle)
            call get_blk_header(block, i,c_array_block)  
            call free_block(old_array, block,
     *                i, array_table, narray_table,
     *                index_table, nindex_table, block_map_table)
            ierr = 0
            return
         endif
      endif

c---------------------------------------------------------------------------
c   Point to the next block in the list and check for a block to reclaim.
c---------------------------------------------------------------------------

      if (next .ne. 0) go to 100

      return
      end

      subroutine free_persistent_blocks(array, array_table, 
     *                    narray_table, index_table, nindex_table,
     *                    block_map_table)
c---------------------------------------------------------------------------
c   Releases all persistent blocks of an array.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'parallel_info.h'
      include 'trace.h'

      integer narray_table
      integer array_table(larray_table_entry, narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer array, flagval, request, blk_data
      integer ptr, i, ierr
      integer status(MPI_STATUS_SIZE)

      logical flag 

      ptr = array_table(c_block_list,array)
      if (ptr .eq. 0) return    ! no blocks

  100 continue
      i = ptr 

c------------------------------------------------------------------------
c   Check the persistence flag.
c------------------------------------------------------------------------

      call get_blk_header(flagval, i,c_block_flags)
      if (and(flagval, block_persistence_flag)
     *                       .ne. 0) then

c-------------------------------------------------------------------------
c   The block may still be communicating.  In that case, we must wait.
c--------------------------------------------------------------------------

         call get_blk_header(request, i,c_block_request)
         if (request .ne. MPI_REQUEST_NULL) then
  200        continue
             call mpi_test(request, flag,
     *                       status, ierr)
             if (flag) then
                call put_blk_header(MPI_REQUEST_NULL, 
     *                      i,c_block_request)
             else 
                call exec_thread_server(0)
                go to 200
             endif

             call blkmgr_remove_block_from_list(comm_list_head,
     *                    comm_list_tail, i, c_comm_list_ptr)
         endif

c-------------------------------------------------------------------------
c   Free the block and remove it from the list.
c   The free_block subroutine also removes this block from the 
c   persistent block linked list.
c-----------------------------------------------------------------------

            call get_blk_header(ptr, i,c_block_list_ptr)
            call get_blk_header(blk_data, i,c_array_block)
            call free_block(array, blk_data, i,
     *                   array_table, narray_table, 
     *                   index_table, nindex_table,
     *                   block_map_table)
      else

c--------------------------------------------------------------------------
c   Find the next block in the list.
c--------------------------------------------------------------------------
      
         call get_blk_header(ptr, i,c_block_list_ptr)
      endif
      if (ptr .ne. 0) go to 100

      return
      end

      subroutine blkmgr_insert_block_in_list(head, tail, blkndx, 
     *                                       list, fifo)
c-------------------------------------------------------------------------
c   Enters a newly-allocated block in a linked list.
c   The "fifo" argument adds the block to the end of the list, otherwise 
c   it is inserted at the head.
c-------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer head, tail, list
      integer blkndx
      integer existing_head
      logical fifo

      if (fifo) then

         if (tail .eq. 0) then

c-------------------------------------------------------------------------
c   List is currently empty.  Create a new head.
c-------------------------------------------------------------------------

            head = blkndx
            tail = blkndx
            call put_blk_header(0, blkndx,list)
         else

c--------------------------------------------------------------------------
c   Add to end of existing list.
c--------------------------------------------------------------------------

            call put_blk_header(blkndx, tail, list)
            tail = blkndx                      ! update tail 
            call put_blk_header(0, blkndx,list)
         endif
      else
         existing_head = head
         head          = blkndx
         call put_blk_header(existing_head, blkndx,list)
      endif

      return
      end

      subroutine blkmgr_remove_block_from_list(head, tail, 
     *                                         blkndx, list)
c--------------------------------------------------------------------------
c   Removes a block from the linked list of blocks.
c--------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer head, tail, blkndx, list
      integer ptr, data
      integer array

      call get_blk_header(array, blkndx,c_array_handle)
      if (head .eq. blkndx) then

c---------------------------------------------------------------------------
c   Removing the head.
c---------------------------------------------------------------------------

         call get_blk_header(ptr, blkndx,list)
         head = ptr
         call put_blk_header(0, blkndx,list) 
         if (blkndx .eq. tail) tail = 0  
      else

c--------------------------------------------------------------------------
c   Find the block whose pointer points to blkndx.
c--------------------------------------------------------------------------

         ptr = head
  100    continue
         if (ptr .eq. 0) then
            print *,'Error: Cannot find block ',blkndx,
     *                   ' to remove from block list.'
            call print_list(head, list)
            call dump_block_ids()
            call abort_job()
         else
            call get_blk_header(data, ptr,list)
            if (data .eq. blkndx) then

c--------------------------------------------------------------------------
c   This block points to the block we wish to remove.  Change its pointer
c   and return.
c--------------------------------------------------------------------------

               call get_blk_header(data, blkndx,list)
               call put_blk_header(data, ptr,list)
               call put_blk_header(0, blkndx,list)

c--------------------------------------------------------------------------
c  If we just removed the tail, this block becomes the new tail.
c--------------------------------------------------------------------------

               if (blkndx .eq. tail) tail = ptr
               return
            endif
         endif

c--------------------------------------------------------------------------
c   Point to the next block in the list, try again.
c--------------------------------------------------------------------------

         call get_blk_header(data, ptr,list)
         ptr = data
         go to 100
      endif
      
      return
      end

      subroutine print_list(head, list)
      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'

      integer head, list
      integer ptr, next
      integer array_data, blk_data, flagval
   
      if (head .eq. 0) then
         print *,'Task ',me,' line ',current_line,
     *     ' List ',list,' has no blocks in list.'
      else
         print *,'Task ',me,' line ',current_line, 
     *     ' head ',head,' list ',list
         ptr = head

  100    continue
         call get_blk_header(next, ptr,list)
         call get_blk_header(array_data, ptr,c_array_handle)
         call get_blk_header(blk_data, ptr,c_array_block)
         call get_blk_header(flagval, ptr,c_block_flags)
         
         print *,'Task ',me,' line ',current_line,'   ptr ',ptr,
     *      '-->',next,' array ',array_data,
     *    ' block ',blk_data,
     *     ' block flags ',flagval
         if (ptr .eq. next) then
            print *,'Error: recursion in list detected.'
            call abort_job()
         endif

         ptr = next
         if (ptr .eq. 0) then
            print *,'Task ',me,' line ',current_line,' End of list'
            return
         endif

         go to 100
      endif

      return
      end 

      subroutine pop_free_stack(stack, blkndx)
c---------------------------------------------------------------------------
c   Pops a block from the "stack" free stack.
c   The block index number (blkndx) is returned.
c   If there is a stack underflow, a -1 is returned.
c---------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      integer stack
      integer blkndx

      if (free_stack_ptr(stack) .lt. stack_start(stack)) then
         blkndx = -1
      else

c---------------------------------------------------------------------------
c   Capture the value of the stack.
c---------------------------------------------------------------------------

         call get_free_stack(blkndx, free_stack_ptr(stack))
         free_stack_ptr(stack) = free_stack_ptr(stack) - 1
         if (blkndx .lt. stack_start(stack) .or. 
     *       blkndx .gt. stack_start(stack)+stack_size(stack)-1) then
            print *,'Error: Invalid data on free stack: ',
     *           ' stack = ',stack,' blkndx ',blkndx,
     *           ' stack data should be between ',
     *          stack_start(stack),' and ',
     *          stack_start(stack)+stack_size(stack)-1
            print *,'stack pointer ',free_stack_ptr(stack)
            call abort_job()
         endif
      endif

      return
      end

      subroutine push_free_stack(stack, blkndx)
c---------------------------------------------------------------------------
c   Pushes a block onto the free stack.
c   The block index number (blkndx) is stored on the stack.
c   If there is a stack overflow, the job is aborted.
c---------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'parallel_info.h'
      integer stack, blkndx

      if (free_stack_ptr(stack) .ge. stack_start(stack) +
     *                                  stack_size(stack)-1) then
         print *,'Task ',me,' line ',current_line,
     *      ' free stack overflow: stack ',stack
         call abort_job()
      else

      if (blkndx .lt. stack_start(stack) .or. 
     *    blkndx .gt. stack_start(stack)+stack_size(stack)-1) then
          print *,'Error: Invalid data submitted to free stack:',
     *     ' stack ',stack,' blkndx ',blkndx
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Increment the stack pointer, store the value.
c---------------------------------------------------------------------------

         free_stack_ptr(stack) = free_stack_ptr(stack) + 1
         call put_free_stack(blkndx, free_stack_ptr(stack))
      endif
      return
      end

      integer function allocate_block(array, array_block, size,
     *                       array_table, narray_table,
     *                       index_table, nindex_table, 
     *                       block_map_table)
c-------------------------------------------------------------------------
c     Finds a free block, allocates it to block "array_block" of "array".
c     The block index is returned in the return code if 
c     successful, otherwise a negative value is returned.
c-------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'

      integer array, array_block, size, ierr
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer i
      integer k, max_iterations, kscrub, scrub_iterations
      integer blkndx
      integer stack, istack

      allocate_block = -1
      max_iterations = 600000   
      scrub_iterations = 100

      if (array_block .lt. 1) then
         print *,'Error in allocate_block: Invalid block number ',
     *           array_block,
     *           ' for array ',array
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Find the smallest stack on which this block will fit.
c--------------------------------------------------------------------------

      do i = 1, nblkmgr_stacks
         if (size .le. stack_datasize(i)) then
            istack = i
            go to 100
         endif
      enddo

      print *,'Error: array, block ',array,array_block,
     *   ' will not fit on any stack.  Size = ',size
      call array_block_summary(array_table, narray_table)
      call dump_block_ids()
      call abort_job()
  100 continue

      do k = 1, max_iterations

c--------------------------------------------------------------------------
c   Pop a block off one of the free stacks.  We use the smallest stack
c   on which we can find a suitable block.
c--------------------------------------------------------------------------

         do i = istack, nblkmgr_stacks
            stack = i
            call pop_free_stack(stack, blkndx)
            if (blkndx .gt. 0) go to 200 

c---------------------------------------------------------------------------
c   Try to scrub a block from the current stack.
c---------------------------------------------------------------------------

            do kscrub = 1, scrub_iterations
               call exec_thread_server(0)  
               call scrub_from_stack(stack, array_table,
     *                     narray_table, index_table, nindex_table,
     *                     block_map_table, ierr)

               if (ierr .ne. 0) then

c-------------------------------------------------------------------------
c   We failed to scrub a block.  Try to reclaim a persistent block that
c   lives on this stack.
c-------------------------------------------------------------------------

                  call reclaim_persistent_block_from_stack(stack,
     *              array_table,narray_table, index_table, nindex_table,
     *              block_map_table, ierr)
               endif
            
               if (ierr .eq. 0) then
                  call pop_free_stack(stack,blkndx)
                  if (blkndx .gt. 0) go to 200
               endif
            enddo      ! kscrub_iterations

         enddo         ! istack
  200    continue

         if (blkndx .gt. 0) then
            i = blkndx

c--------------------------------------------------------------------------
c   Set up the newly allocated block.
c--------------------------------------------------------------------------

            call put_blk_header(array, i,c_array_handle)
            call put_blk_header(array_block,i,c_array_block) 
            call put_blk_header(mpi_request_null,i,c_block_request)
            call put_blk_header(block_busy_flag,i,c_block_flags)
            call put_blk_header(size,i,c_block_size)   
            call put_blk_header(0,i,c_block_list_ptr)
            allocate_block = i 

            call set_block_id(array, array_block, stack, i)
            call put_blk_in_use(i)
            return 
         else
            if (comm_list_head .ne. 0) then
               call exec_thread_server(0)
               call scrub_blocks(size, .false., array_table,
     *                   narray_table, index_table, nindex_table,
     *                   block_map_table, ierr)  
c               print *,'Task ',me,' ALLOC SCRUB array,block ',
c     *              array,array_block,' k, ierr ',k, ierr
            else
               ierr = 1
            endif

            if (ierr .ne. 0) then     ! could not scrub any.
                call reclaim_persistent_block(size, array_table, 
     *               narray_table, index_table, nindex_table,
     *               block_map_table, ierr)
c                 print *,'Task ',me,' ALLOC RECLAIM array,block ',
c     *             array,array_block,' k, ierr ',k, ierr
            endif

            if (ierr .ne. 0) then
               call exec_thread_server(0)   ! do one pass of thread server.
            endif
         endif
      enddo
 
      print *,'Task ',me,' line ',current_line,
     *  ' Cannot allocate block for array ',array,' stack ',
     *   stack
      call array_block_summary(array_table, narray_table)
      call dump_block_ids()
      call abort_job()

      return
      end

      integer function allocate_block_light(array, array_block, size,
     *                       array_table, narray_table,
     *                       index_table, nindex_table, 
     *                       block_map_table)
c-------------------------------------------------------------------------
c     Finds a free block, allocates it to block "array_block" of "array".
c     The block index is returned in the return code if 
c     successful, otherwise a negative value is returned.
c
c     This routine does not go through the entire allocation algorithm,
c     it may be used when an instruction does not requre the block to 
c     be allocated to run.
c-------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'trace.h'
      include 'mpif.h'
      include 'parallel_info.h'

      integer array, array_block, size, ierr
      integer narray_table
      integer array_table(larray_table_entry,narray_table)
      integer nindex_table
      integer index_table(lindex_table_entry,nindex_table)
      integer block_map_table(lblock_map_entry,*)
      integer i
      integer k, max_iterations, kscrub, scrub_iterations
      integer blkndx
      integer stack, istack

      allocate_block_light = -1

      if (array_block .lt. 1) then
         print *,'Error in allocate_block: Invalid block number ',
     *           array_block,
     *           ' for array ',array
         call abort_job()
      endif

c--------------------------------------------------------------------------
c   Find the smallest stack on which this block will fit.
c--------------------------------------------------------------------------

      do i = 1, nblkmgr_stacks
         if (size .le. stack_datasize(i)) then
            istack = i
            go to 100
         endif
      enddo

      print *,'Error: array, block ',array,array_block,
     *   ' will not fit on any stack.  Size = ',size
      call array_block_summary(array_table, narray_table)
      call dump_block_ids()
      call abort_job()
  100 continue

c--------------------------------------------------------------------------
c   Pop a block off one of the free stacks.  We use the smallest stack
c   on which we can find a suitable block.
c--------------------------------------------------------------------------

       do stack = istack, nblkmgr_stacks     
          call pop_free_stack(stack, blkndx)
c          print *,'TASK ',me,' ALLOC: stack , blkndx ',stack,blkndx
          if (blkndx .gt. 0) go to 200 

c---------------------------------------------------------------------------
c   Try to scrub a block from the current stack.
c---------------------------------------------------------------------------

          call exec_thread_server(0)  
          call scrub_from_stack(stack, array_table,
     *                     narray_table, index_table, nindex_table,
     *                     block_map_table, ierr)
c          print *,'TASK ',me,' ALLOC AFTER SCRUB ierr ',ierr
          if (ierr .eq. 0) then
             blkndx = ierr 
             call pop_free_stack(stack,blkndx)
c             print *,'TASK ',me,' ALLOC AFTER POP: stack, blkndx ',
c     *          stack,blkndx 
             go to 200
          endif
      enddo    ! stack

c---------------------------------------------------------------------------
c   Couldn't get a block off any of the free stacks.  Return without 
c   allocating the block.
c---------------------------------------------------------------------------

      allocate_block_light = -1
c      print *,'TASK ',me,' ALLOC_LIGHT RETURN -1'
      return
            
  200 continue 

c---------------------------------------------------------------------------
c   We got a block.  Set it up for use.
c---------------------------------------------------------------------------

      i = blkndx

c--------------------------------------------------------------------------
c   Set up the newly allocated block.
c--------------------------------------------------------------------------

      call put_blk_header(array, i,c_array_handle)
      call put_blk_header(array_block,i,c_array_block) 
      call put_blk_header(mpi_request_null,i,c_block_request)
      call put_blk_header(block_busy_flag,i,c_block_flags)
      call put_blk_header(size,i,c_block_size)   
      call put_blk_header(0,i,c_block_list_ptr)
      allocate_block_light = i 

      call set_block_id(array, array_block, stack, i)
      call put_blk_in_use(i)
c      print *,'TASK ',me,' ALLOC LIGHT SUCCEEDED: RETURN ',i  
      return 
      end

      integer function allocate_scratch_block(base, index, handle,
     *                         array_table, narray_table, 
     *                         index_table, nindex_table, 
     *                         block_map_table)
c-------------------------------------------------------------------------
c     Finds a free block, allocates it as a scratch block.
c     The index of the scratch area relative to "base" is returned along 
c     with a zero return code if 
c     successful, otherwise a non-zero value is returned.
c-------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      include 'blkmgr.h'
      include 'mpif.h'
      include 'trace.h'
      include 'machine_types.h'

      integer narray_table
      integer nindex_table
      integer array_table(larray_table_entry, narray_table)
      integer index_table(lindex_table_entry,nindex_table)
      integer block_map_table(*)
      integer handle
      integer i, j, k, flagval
      integer size
      integer*8 index, block_addr
      integer*8 get_index_from_base
      integer*8 c_loc64, ixx
      integer maxblocksize 
      integer stack, data
      integer kscrub, ierr
      integer scrub_iterations
      parameter (scrub_iterations = 100)

      double precision base(1)

      allocate_scratch_block = 1

      if (.not. blkmgr_init_flag) then
         print *,'Error: blkmgr has not been initialized.'
         call abort_job()
      endif

c---------------------------------------------------------------------------
c   Find the next available block.
c   Scratch blocks are always allocated off the largest (blocksize) stack.
c---------------------------------------------------------------------------

      maxblocksize = stack_blocksize(1)
      stack = 1
      do i = 2, nblkmgr_stacks
         if (stack_blocksize(i) .gt. maxblocksize) then
            maxblocksize = stack_blocksize(i)
            stack = i
         endif
      enddo

      call pop_free_stack(stack, handle)  
      if (handle .le. 0) then
         do kscrub = 1, scrub_iterations
            call scrub_from_stack(stack, array_table,
     *                  narray_table, index_table, nindex_table,
     *                  block_map_table, ierr)
            if (ierr .ne. 0) then

c-------------------------------------------------------------------------
c   We failed to scrub a block.  Try to reclaim a persistent block that
c   lives on this stack.
c-------------------------------------------------------------------------

               call reclaim_persistent_block_from_stack(stack,
     *           array_table,narray_table, index_table, nindex_table,
     *           block_map_table, ierr)
            endif

            if (ierr .eq. 0) then
               call pop_free_stack(stack,handle)
               if (handle .gt. 0) go to 200
            endif
         enddo      ! kscrub_iterations

         print *,'Error: Unable to allocate a scratch block, line ',
     *       current_line,' stack, handle ',stack, handle
         call abort_job()
      endif

  200 continue

c---------------------------------------------------------------------------
c   Turn on the block_scratch_flag
c---------------------------------------------------------------------------

      i = handle
      call put_blk_header(mpi_request_null,i,c_block_request)
      call get_blk_header(data, i,c_block_flags)
      flagval = or(data, block_scratch_flag)
      call put_blk_header(flagval,i,c_block_flags)

c---------------------------------------------------------------------------
c   Determine the index to return.
c---------------------------------------------------------------------------
      
      call get_blk_addr(block_addr, i)
      index = get_index_from_base(block_addr, base, 2) 
 
      allocate_scratch_block = 0
      return
      end

      subroutine free_scratch_block(handle)
c-------------------------------------------------------------------------
c   Releases the scratch block whose handle is "handle".
c-------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer i, stack, handle, flagval
      integer flag

c--------------------------------------------------------------------------
c   Determine the stack for this block.
c--------------------------------------------------------------------------

      stack = nblkmgr_stacks
      do i = 1, nblkmgr_stacks-1
         if (handle .ge. stack_start(i) .and. 
     *       handle .lt. stack_start(i+1)) stack = i
      enddo

c---------------------------------------------------------------------------
c   Make sure this is actually a scratch block.
c---------------------------------------------------------------------------

      call get_blk_header(flag, handle,c_block_flags)
      flagval =  and(flag, block_scratch_flag)
      if (flagval .eq. 0) then
         print *,'Error: Attempt to release scratch block ',handle
         print *,'The block is not a scratch block. flags = ',
     *        flag
         call abort_job()
      endif

      call put_blk_header(xor(flag,block_scratch_flag),
     *                    handle,c_block_flags)
      
      call push_free_stack(stack, handle)
      return
      end

      integer function get_stack_index(blkndx)
c---------------------------------------------------------------------------
c   Determines the stack index of a given block.
c---------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer blkndx

      integer i, stack

      stack = nblkmgr_stacks
      do i = 1, nblkmgr_stacks-1
         if (blkndx .ge. stack_start(i)) stack = i
      enddo

      get_stack_index = stack
      return
      end

      subroutine server_monitor_write_log(blkndx, optype, server)
c--------------------------------------------------------------------------
c   Write a record to the server monitor logfile.
c   These files are saved for later analysis.
c--------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'server_monitor.h' 
      include 'parallel_info.h'
      include 'trace.h'

      integer blkndx, server
      character*(*) optype
      integer nw, ierr
      integer flagval
      integer array, block
      
      character*24 timestamp, get_timestamp

      if (.not. server_monitor_on) return

c--------------------------------------------------------------------------
c   For "request" types, make sure the block was actually used in a request.
c--------------------------------------------------------------------------

      if (optype .eq. 'r') then
         call get_blk_header(flagval, blkndx, c_block_flags)
         if (and(flagval, block_request_flag) .eq. 0)
     *       return

c---------------------------------------------------------------------------
c   Turn off the request flag saved by blkmgr.
c---------------------------------------------------------------------------

         call get_blk_header(array,blkndx,c_array_handle)
         call get_blk_header(block, blkndx,c_array_block)
         call set_block_request_flag(array, block,
     *                               blkndx, 0)
      endif

      timestamp = get_timestamp()
      call get_blk_header(nw, blkndx, c_block_size)
      write (server_logunit,*) timestamp,me,server,nw,optype,
     *                         current_line
      return
      end

      subroutine get_block_indices(blkndx, ind)
c--------------------------------------------------------------------------
c   Returns the indices for block "blkndx" in the "ind" array.
c--------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer blkndx
      integer ind(*)

      integer i, nind

      call get_blk_header(nind, blkndx, c_block_nindex)
      do i = 1, nind
         call get_blk_header(ind(i), blkndx, c_block_indices+i-1)
      enddo

      return
      end

      subroutine get_block_segments(blkndx, seg)
c--------------------------------------------------------------------------
c   Returns the segments for block "blkndx" in the "seg" array.
c--------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      integer blkndx
      integer seg(*)
      integer i, nind

      call get_blk_header(nind, blkndx, c_block_nindex)
      do i = 1, nind
         call get_blk_header(seg(i),blkndx, c_block_segments+i-1)
      enddo

      return
      end

      subroutine free_all_blocks()
c----------------------------------------------------------------------------
c  Clears all remaining MPI requests from the blkmgr data blocks.
c----------------------------------------------------------------------------
      implicit none
      include 'maxdim.h'
      include 'blkmgr.h'
      include 'parallel_info.h'
      include 'mpif.h'
      integer i, ierr, request, array
      integer status(MPI_STATUS_SIZE)
      logical flag

#ifndef ALTIX
      return
#endif

      do i = 1, blkmgr_blocks
         call get_blk_header(request, i, c_block_request)
         if (request .ne. MPI_REQUEST_NULL) then
            call mpi_request_free(request, ierr)
         endif
      enddo

      return
      end

      subroutine subindex_seg_match(index, index_table,
     *                  nindex_table, segment_table,
     *                  nsegment_table, seg)
c--------------------------------------------------------------------------
c   Determines the superindex segment that contains the current segment 
c   of the subindex.
c--------------------------------------------------------------------------

      implicit none
      include 'interpreter.h'
      integer index
      integer nindex_table, nsegment_table, seg
      integer index_table(lindex_table_entry, nindex_table)
      integer segment_table(lsegment_table_entry, nsegment_table)

      integer superindex, r1, r2, r1sub, r2sub
      integer i, current_seg
 
c--------------------------------------------------------------------------
c   Look up the superindex and beginning and ending ranges of the subindex 
c   segment.
c--------------------------------------------------------------------------

      superindex = index_table(c_subindex_ptr, index)
      current_seg = index_table(c_current_seg, index)
      call get_index_segment(index, current_seg, segment_table,
     *                             nsegment_table, index_table,
     *                             nindex_table, r1sub, r2sub) 

c--------------------------------------------------------------------------
c   Find the superindex segment containing r1sub, r2sub.
c--------------------------------------------------------------------------

      do i = 1, nsegment_table
         if (segment_table(c_index,i) .eq. superindex) then
            r1 = segment_table(c_range1,i)
            r2 = segment_table(c_range2,i)
            if (r1 .le. r1sub .and. r2 .ge. r2sub) then
               seg = segment_table(c_segment,i)
               return
            endif 
         endif
      enddo

      print *,'Error: subindex_seg_match cannot find matching segment'
      print *,'index, superindex, r1sub, r2sub ',
     *      index, superindex, r1sub, r2sub
      call abort_job() 
      return
      end 
