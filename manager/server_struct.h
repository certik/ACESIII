
c------------------------------------------------------------------------------
c   Description of server_table.
c-----------------------------------------------------------------------------

      integer c_server_memloc
      integer c_server_diskloc
      integer c_server_size
      integer c_server_array
      integer c_server_nind
      integer c_server_flags
      integer c_server_busy_node
      integer c_server_file
      integer c_server_iblk
      integer c_server_bsegs
      integer c_server_esegs
      integer lserver_table_entry
      parameter (c_server_memloc    = 1)
      parameter (c_server_diskloc   = 2)
      parameter (c_server_size      = 3)
      parameter (c_server_array     = 4)
      parameter (c_server_nind      = 5)
      parameter (c_server_flags     = 6)
      parameter (c_server_busy_node = 7)
      parameter (c_server_file      = 8)
      parameter (c_server_iblk      = 9)
      parameter (c_server_bsegs     = 10)
      parameter (c_server_esegs  = c_server_bsegs + mx_array_index)
      parameter (lserver_table_entry =
     *            c_server_esegs + mx_array_index-1)

c---------------------------------------------------------------------------
c   Flags used in server data.
c---------------------------------------------------------------------------

      integer server_busy_flag
      parameter (server_busy_flag = 1)
      integer server_dirty_flag
      parameter (server_dirty_flag = 2)
      integer server_deleted_flag
      parameter (server_deleted_flag = 4)
      integer preparesum_copy_flag
      parameter (preparesum_copy_flag = 8)
      integer integral_calculation_flag
      parameter (integral_calculation_flag = 16)
      integer restore_flag
      parameter (restore_flag = 64)
      integer readonly_flag
      parameter (readonly_flag = 10101)
      integer writeonly_flag
      parameter (writeonly_flag = 20202)

c---------------------------------------------------------------------------
c   Common block used in server routines.
c---------------------------------------------------------------------------

      integer mx_server_blocksizes
      parameter (mx_server_blocksizes = 10)
      integer mx_served_arrays
      parameter (mx_served_arrays = 500)

      integer*8 base_mem_addr
      integer*8 server_table_base_addr
      integer nserver_table_entries
      integer nmessage_buffers
      integer nserver_blocksizes
      integer server_blocksizes(mx_server_blocksizes)
      integer server_unit(mx_server_blocksizes)
      integer next_server_diskloc(mx_server_blocksizes)
      integer server_mem_blocksize
      integer served_array_table, nserved_arrays
      integer served_array_status
      integer served_numblocks
      integer served_array_entry
      integer barrier_seqno
      integer barrier_count
      integer barrier_msg_count
      integer nbarrier_msgs
      integer request_buffer
      logical barrier_in_progress
      character*500 server_filename(mx_server_blocksizes)
      common /server/base_mem_addr, server_table_base_addr,
     *               nserver_table_entries, nserver_blocksizes,
     *               server_blocksizes, nmessage_buffers,
     *               server_unit, next_server_diskloc,
     *               server_mem_blocksize, server_filename,
     *               served_array_table(mx_served_arrays), 
     *               served_array_status(mx_served_arrays),
     *               served_numblocks(mx_served_arrays),
     *               served_array_entry(mx_served_arrays),
     *               nserved_arrays, barrier_seqno, barrier_count,
     *               barrier_msg_count, nbarrier_msgs,
     *               barrier_in_progress, request_buffer

c---------------------------------------------------------------------------
c   Message buffer stack
c---------------------------------------------------------------------------

      integer mx_server_msg
      parameter (mx_server_msg = 100)
      
      integer c_msg_type
      parameter (c_msg_type = 1)
      integer c_msg_array
      parameter (c_msg_array = 2)
      integer c_msg_source 
      parameter (c_msg_source = 3)
      integer c_msg_tag
      parameter (c_msg_tag = 4)
      integer c_msg_nind
      parameter (c_msg_nind = 5)
      integer c_msg_size
      parameter (c_msg_size = 6)
      integer c_msg_memptr
      parameter (c_msg_memptr = 7)
      integer c_msg_msgbuffer
      parameter (c_msg_msgbuffer = 8)
      integer c_msg_state
      parameter (c_msg_state = 9)
      integer c_msg_cause
      parameter (c_msg_cause = 10)
      integer c_msg_request
      parameter (c_msg_request = 11)
      integer c_msg_flag
      parameter (c_msg_flag = 12)
      integer c_server_list_ptr
      parameter (c_server_list_ptr = 13)
      integer c_msg_seqno
      parameter (c_msg_seqno = 14)
      integer c_msg_current_line 
      parameter (c_msg_current_line = 15)
      integer c_msg_stat_key
      parameter (c_msg_stat_key = 16)
      integer c_msg_stptr
      parameter (c_msg_stptr = 17)
      integer c_msg_iblk
      parameter (c_msg_iblk = 18)  
      integer c_msg_bsegs
      parameter (c_msg_bsegs = 19)
      integer c_msg_esegs
      parameter (c_msg_esegs = c_msg_bsegs+mx_array_index)
      integer c_msg_bsegs2
      parameter (c_msg_bsegs2 = c_msg_esegs + mx_array_index)
      integer c_msg_esegs2
      parameter (c_msg_esegs2 = c_msg_bsegs2 + mx_array_index)
      integer lserver_msg_entry
      parameter (lserver_msg_entry = c_msg_esegs2+mx_array_index-1) 

      integer server_node_ptr
      integer server_msg_node
      integer server_msg
      integer server_seqno
      integer*8 blk_to_list_offset
      integer server_work_list_head, server_work_list_tail
      common /server_messages/server_msg_node(mx_server_msg),
     *         server_node_ptr, server_seqno,
     *         server_work_list_head, server_work_list_tail,
     *         server_msg(lserver_msg_entry,mx_server_msg),
     *         blk_to_list_offset(2,mx_server_msg)

c---------------------------------------------------------------------------
c   Server memory block table
c---------------------------------------------------------------------------

      integer mx_server_memblocks
      parameter (mx_server_memblocks = 50000)
      
      integer nserver_memblocks
      integer nclean_blocks
      integer server_table_ptr
      integer ndirty, dirty_threshold, max_backup
      integer dirty_list_ptr,dirty_list_head, dirty_list_tail
      integer clean_blocks, clean_block_ptr

      common /server_memblock/nserver_memblocks, nclean_blocks,
     *      server_table_ptr(mx_server_memblocks), 
     *      dirty_list_ptr(mx_server_memblocks), dirty_list_head, 
     *      dirty_list_tail, ndirty, dirty_threshold, max_backup,
     *      clean_blocks(mx_server_memblocks), clean_block_ptr

c---------------------------------------------------------------------------
c   Server processing states
c---------------------------------------------------------------------------

      integer begin_state
      parameter (begin_state = 101)
      integer null_state
      parameter (null_state = 102)
      integer quit_state
      parameter (quit_state = 103)
      integer recv_block_state
      parameter (recv_block_state = 104)
      integer wait_for_backup_state
      parameter (wait_for_backup_state = 105)
      integer wait_for_restore_state
      parameter (wait_for_restore_state = 106)
      integer wait_for_block_state
      parameter (wait_for_block_state = 107)
      integer wait_for_send_state
      parameter (wait_for_send_state = 108)
      integer enter_backup_state
      parameter (enter_backup_state = 109)
      integer enter_restore_state
      parameter (enter_restore_state = 110)
      integer enter_backup_and_restore_state
      parameter (enter_backup_and_restore_state = 111)
      integer request_cleanup_state
      parameter (request_cleanup_state = 112)
      integer prequest_intermediate_state
      parameter (prequest_intermediate_state = 113)

      integer null_cause
      parameter (null_cause = 201)
      integer backup_agent_cause
      parameter (backup_agent_cause = 202)
      integer restore_cause
      parameter (restore_cause = 203)
      integer busy_cause
      parameter (busy_cause = 204)

      integer blocks_list_done
      parameter (blocks_list_done = 90909)
      integer list_to_blocks_done
      parameter (list_to_blocks_done = 80808)

c----------------------------------------------------------------------------
c  Tags and message types
c---------------------------------------------------------------------------

      integer server_readytag
      parameter (server_readytag = 3336)
      integer server_prepare_msgtype
      parameter (server_prepare_msgtype = 3333)
      integer server_request_msgtype
      parameter (server_request_msgtype = 3334)
      integer server_prepare_increment
      parameter (server_prepare_increment = 3340)
      integer server_quit_msgtype
      parameter (server_quit_msgtype = 3338)
      integer server_barrier_signal
      parameter (server_barrier_signal = 3341)
      integer server_stat_data_tag 
      parameter (server_stat_data_tag = 3342)
      integer server_copy_msg
      parameter (server_copy_msg = 3343)
      integer server_delete_msg
      parameter (server_delete_msg = 3344)
      integer server_prequest_msg
      parameter (server_prequest_msg = 3345)
      integer server_blocks_to_list_msg
      parameter (server_blocks_to_list_msg = 3346)
      integer server_checkpoint_msg 
      parameter (server_checkpoint_msg = 3347)
      integer server_restart_msg
      parameter (server_restart_msg = 3348)
      integer server_commit_msg
      parameter (server_commit_msg = 3349)
      integer server_list_to_blocks_msg 
      parameter (server_list_to_blocks_msg = 3350)

