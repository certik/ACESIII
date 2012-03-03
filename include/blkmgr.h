
      integer block_busy_flag
      integer block_data_present_flag 
      integer block_persistence_flag
      integer block_computed_flag
      integer block_created_flag
      integer block_request_outstanding_flag
      integer block_scratch_flag
      integer block_scrub_flag
      integer block_request_flag
      integer block_prefetch_flag
      parameter (block_busy_flag = 1)
      parameter (block_data_present_flag = 2)
      parameter (block_persistence_flag = 4)
      parameter (block_computed_flag = 8)
      parameter (block_created_flag = 16)
      parameter (block_request_outstanding_flag = 32)
      parameter (block_scratch_flag = 64)
      parameter (block_scrub_flag = 128)
      parameter (block_request_flag = 256)
      parameter (block_prefetch_flag = 512)

      integer max_blkmgr_blocks
      parameter (max_blkmgr_blocks = 20000)
      integer c_array_handle
      parameter (c_array_handle = 1)
      integer c_array_block 
      parameter (c_array_block = 2) 
      integer c_block_flags
      parameter (c_block_flags = 3)
      integer c_block_nindex
      parameter (c_block_nindex = 4)
      integer c_block_size
      parameter (c_block_size = 5)
      integer c_block_request
      parameter (c_block_request = 6)
      integer c_block_list_ptr
      parameter (c_block_list_ptr = 7)
      integer c_comm_list_ptr
      parameter (c_comm_list_ptr = 8)
      integer c_persistent_ptr
      parameter (c_persistent_ptr = 9)
      integer c_block_stack
      parameter (c_block_stack = 10)
      integer c_block_indices
      parameter (c_block_indices = 11)
      integer c_block_segments
      parameter (c_block_segments = c_block_indices+mx_array_index)
      integer lblk_header
      parameter (lblk_header = c_block_segments + mx_array_index-1)

      integer lblock_id_data
      parameter (lblock_id_data = 3)
      
      integer max_stacks
      parameter (max_stacks = 10)

      integer blkmgr_blocks
      integer nblk_in_use
      integer stack_blocksize   ! size of block_id area + data area
      integer stack_datasize    ! size of data area
      integer*8 stack_base_addr
      integer blkmgr_comm
      integer blkmgr_comm_rank
      integer blkmgr_next_block
      logical blkmgr_init_flag
      integer free_stack_ptr
      integer stack_start
      integer nblkmgr_stacks
      integer comm_list_head, comm_list_tail
      integer persistent_list_head, persistent_list_tail
      integer stack_size
      integer*8 addr_blk_header
      integer*8 addr_free_stack
      integer*8 addr_blk_addr
      integer*8 addr_blk_in_use
      integer*8 first_blk_addr, last_blk_addr

      common /blkmgr/addr_blk_header, addr_free_stack, addr_blk_addr,
     *               addr_blk_in_use, first_blk_addr, last_blk_addr,
     *               stack_base_addr(max_stacks),blkmgr_blocks, 
     *               nblk_in_use,
     *               blkmgr_init_flag, blkmgr_comm, blkmgr_comm_rank,
     *               blkmgr_next_block, 
     *               free_stack_ptr(max_stacks), nblkmgr_stacks,
     *               stack_start(max_stacks),
     *               stack_size(max_stacks), 
     *               stack_blocksize(max_stacks),
     *               stack_datasize(max_stacks),
     *               comm_list_head, comm_list_tail,
     *               persistent_list_head, persistent_list_tail
