      integer max_blocks
      parameter (max_blocks = 1000)
      include 'maxdim.h'

c--------------------------------------------------------------------------
c   Parameters describing the optable entries.
c--------------------------------------------------------------------------  
 
      integer c_opcode
      integer c_op1_array
      integer c_op2_array
      integer c_result_array
      integer c_scalar_op1
      integer c_scalar_op2
      integer c_scalar_result
      integer c_pardo_lock_index
      integer c_ind1
      integer c_ind2
      integer c_ind3
      integer c_ind4
      integer c_user_sub
      integer c_instr_timer
      integer c_opblock
      integer c_opblkndx
      integer c_oploop
      integer c_lineno
      integer c_pardo_batch
      integer c_pardo_max_batch
      integer c_pardo_chunk_size
      integer c_pardo_batch_end
      integer c_pardo_next_batch_start
      integer c_pardo_signal
      integer c_server_stat_key
      parameter (c_opcode = 1)
      parameter (c_op1_array = 2)
      parameter (c_op2_array = 3)
      parameter (c_result_array = 4)
      parameter (c_ind1   = 5)
      parameter (c_ind2   = 6)
      parameter (c_ind3   = 7)
      parameter (c_ind4   = 8)
      parameter (c_user_sub      = c_ind1+mx_array_index)
      parameter (c_instr_timer   = c_user_sub+1)
      parameter (c_opblkndx      = c_instr_timer+1)
      parameter (c_opblock       = c_opblkndx+1)
      parameter (c_oploop        = c_opblock+1)
      parameter (c_scalar_op1    = c_oploop+1)
      parameter (c_scalar_op2    = c_scalar_op1+1)
      parameter (c_scalar_result = c_scalar_op2+1)
      parameter (c_pardo_batch   = c_scalar_result+1)
      parameter (c_pardo_max_batch = c_pardo_batch+1)
      parameter (c_pardo_signal  = c_pardo_max_batch+1)
      parameter (c_server_stat_key = c_pardo_signal+1)
      parameter (c_lineno        = c_server_stat_key+1)

      integer loptable_entry
      parameter (loptable_entry = c_lineno)

c--------------------------------------------------------------------------
c   "Equivalenced" optable indices for pardo load-balancing data
c--------------------------------------------------------------------------

      integer c_pardo_info
      parameter (c_pardo_info       = c_scalar_op1)
      parameter (c_pardo_lock_index = c_user_sub)
      parameter (c_pardo_chunk_size = c_scalar_op1)
      parameter (c_pardo_batch_end  = c_scalar_op2)
      parameter (c_pardo_next_batch_start = c_scalar_result)

c---------------------------------------------------------------------------
c   Parameters describing the index_table entries.
c---------------------------------------------------------------------------

      integer lindex_table_entry
      parameter (lindex_table_entry = 8)
      integer c_index_size
      parameter (c_index_size = 1)
      integer c_nsegments
      parameter (c_nsegments = 2)
      integer c_current_seg
      parameter (c_current_seg = 3)
      integer c_bseg
      parameter (c_bseg = 4)
      integer c_eseg
      parameter (c_eseg = 5)
      integer c_index_type
      parameter (c_index_type = 6)
      integer c_next_seg
      parameter (c_next_seg = 7)
      integer c_subindex_ptr
      parameter (c_subindex_ptr = 8)

c---------------------------------------------------------------------------
c   Parameters describing the segment table entries.
c---------------------------------------------------------------------------

      integer lsegment_table_entry
      parameter (lsegment_table_entry = 6)
      integer c_index 
      parameter (c_index = 1)
      integer c_segment
      parameter (c_segment = 2)
      integer c_range1
      parameter (c_range1 = 3)
      integer c_range2
      parameter (c_range2 = 4)
      integer c_subseg1
      parameter (c_subseg1 = 5)
      integer c_subseg2
      parameter (c_subseg2 = 6)

c---------------------------------------------------------------------------
c   Parameters describing the array_table entries.
c---------------------------------------------------------------------------

      integer c_nindex
      integer c_array_type
      integer c_numblks
      integer c_index_array1
      integer c_index_array2
      integer c_index_array3
      integer c_index_array4
      integer c_index_original
      integer c_index_range1
      integer c_index_range2
      integer c_block_map
      integer c_scalar_index
      integer c_create_flag
      integer c_put_flag
      integer c_prepare_flag
      integer c_current_blkndx
      integer c_block_list
      integer c_array_stack
      integer c_array_status
      parameter (c_nindex = 1)
      parameter (c_array_type   = 2)
      parameter (c_numblks      = 3)
      parameter (c_index_array1 = 4) 
      parameter (c_index_array2 = 5) 
      parameter (c_index_array3 = 6) 
      parameter (c_index_array4 = 7) 
      parameter (c_index_range1 = c_index_array1+mx_array_index)
      parameter (c_index_original = c_index_range1)
      parameter (c_index_range2 = c_index_range1+mx_array_index)
      parameter (c_block_map    = c_index_range2+mx_array_index) 
      parameter (c_scalar_index = c_block_map+1)
      parameter (c_create_flag  = c_scalar_index+1)
      parameter (c_put_flag     = c_create_flag + 1)
      parameter (c_prepare_flag = c_put_flag + 1)
      parameter (c_current_blkndx = c_prepare_flag+1)
      parameter (c_block_list   = c_current_blkndx+1)
      parameter (c_array_stack  = c_block_list+1)
      parameter (c_array_status = c_array_stack+1)
      
      integer larray_table_entry
      parameter (larray_table_entry = c_array_status)
      
c--------------------------------------------------------------------------
c   Parameters describing the block map entry for each block.
c--------------------------------------------------------------------------

      integer lblock_map_entry
      parameter (lblock_map_entry = 2+mx_array_index)

      integer c_processor
      integer c_block_map_seg
      integer c_bmap_blkndx

      parameter (c_processor = 1)
      parameter (c_bmap_blkndx = 2)
      parameter (c_block_map_seg = 3)

c--------------------------------------------------------------------------
c   Parameters describing the eblock table entry for each block.
c--------------------------------------------------------------------------

      integer leblock_entry

      integer c_eblock_blksize
      integer c_eblock_memloc
      integer c_eblock_segs
      integer c_eblock_range1
      integer c_eblock_range2
      parameter (c_eblock_blksize = 1)
      parameter (c_eblock_segs = 2)
      parameter (c_eblock_range1 = c_eblock_segs + mx_array_index)
      parameter (c_eblock_range2 = c_eblock_range1 + mx_array_index)
      parameter (c_eblock_memloc = c_eblock_range2 + mx_array_index)
      parameter (leblock_entry = c_eblock_memloc)

c--------------------------------------------------------------------------
c   Miscellaneous parameters.
c--------------------------------------------------------------------------

      integer contraction_op
      integer sum_op
      integer reindex_op
      integer do_op
      integer enddo_op
      integer get_op
      integer user_sub_op
      integer put_op
      integer go_to_op
      integer create_op
      integer delete_op
      integer call_op
      integer return_op
      integer jz_op
      integer stop_op
      integer sp_add_op
      integer sp_sub_op   
      integer sp_mult_op   
      integer sp_div_op  
      integer sp_equal_op 
      integer sp_nequal_op
      integer sp_ge_op       
      integer sp_le_op  
      integer sp_gt_op 
      integer sp_lt_op   
      integer sp_ldi_op
      integer sp_ldindex_op
      integer sp_ldi_sym_op
      integer pardo_op
      integer endpardo_op
      integer exit_op
      integer assignment_op
      integer cycle_op
      integer self_multiply_op 
      integer subtract_op  
      integer collective_sum_op  
      integer divide_op  
      integer prepare_op
      integer request_op
      integer compute_integrals_op
      integer put_replace_op
      integer tensor_op
      integer fl_add_op
      integer fl_sub_op
      integer fl_mult_op
      integer fl_div_op
      integer fl_eq_op
      integer fl_ne_op
      integer fl_ge_op
      integer fl_le_op
      integer fl_gt_op
      integer fl_lt_op
      integer fl_load_value_op
      integer prepare_increment_op
      integer allocate_op
      integer deallocate_op
      integer destroy_op
      integer prequest_op
      integer where_op
      integer last_opcode

      parameter (contraction_op = 101)
      parameter (sum_op         = 102)
      parameter (reindex_op     = 103)
      parameter (do_op          = 104)
      parameter (enddo_op       = 105)
      parameter (get_op         = 106)
      parameter (user_sub_op    = 107)
      parameter (put_op         = 108)
      parameter (go_to_op       = 109)
      parameter (create_op      = 110)
      parameter (delete_op      = 111)
      parameter (call_op        = 112)
      parameter (return_op      = 113)
      parameter (jz_op          = 114)
      parameter (stop_op        = 115)
      parameter (sp_add_op      = 116)
      parameter (sp_sub_op      = 117)
      parameter (sp_mult_op     = 118)
      parameter (sp_div_op      = 119)
      parameter (sp_equal_op    = 120)
      parameter (sp_nequal_op   = 121)
      parameter (sp_ge_op       = 122)
      parameter (sp_le_op       = 123)
      parameter (sp_gt_op       = 124)
      parameter (sp_lt_op       = 125)
      parameter (sp_ldi_op      = 126)
      parameter (sp_ldindex_op  = 127)
      parameter (pardo_op       = 128)
      parameter (endpardo_op    = 129)
      parameter (exit_op        = 130)
      parameter (assignment_op  = 131)      ! =
      parameter (cycle_op       = 132)
      parameter (self_multiply_op = 134)    ! *= scalar
      parameter (subtract_op    = 135)      ! -
      parameter (collective_sum_op = 136)
      parameter (divide_op      = 137)
      parameter (prepare_op     = 138)
      parameter (request_op     = 139)
      parameter (compute_integrals_op = 140)
      parameter (put_replace_op = 141)
      parameter (tensor_op      = 142)
      parameter (fl_add_op      = 146)
      parameter (fl_sub_op      = 147)
      parameter (fl_mult_op     = 148)
      parameter (fl_div_op      = 149)
      parameter (fl_eq_op       = 150)
      parameter (fl_ne_op       = 151)
      parameter (fl_ge_op       = 152)
      parameter (fl_le_op       = 153)
      parameter (fl_gt_op       = 154)
      parameter (fl_lt_op       = 155)
      parameter (fl_load_value_op  = 157)
      parameter (prepare_increment_op = 158)
      parameter (allocate_op    = 159)
      parameter (deallocate_op  = 160)
      parameter (sp_ldi_sym_op  = 161)
      parameter (destroy_op     = 162)
      parameter (prequest_op    = 163)
      parameter (where_op       = 164)
      parameter (last_opcode    = where_op)

c----------------------------------------------------------------------------
c   Character descriptions of operations.
c----------------------------------------------------------------------------

      character*12 contraction_opdesc
      character*12 sum_opdesc
      character*12 reindex_opdesc
      character*12 do_opdesc
      character*12 enddo_opdesc
      character*12 get_opdesc
      character*12 user_sub_opdesc
      character*12 put_opdesc
      character*12 go_to_opdesc
      character*12 create_opdesc
      character*12 delete_opdesc
      character*12 call_opdesc
      character*12 return_opdesc
      character*12 jz_opdesc
      character*12 stop_opdesc
      character*12 sp_add_opdesc
      character*12 sp_sub_opdesc   
      character*12 sp_mult_opdesc   
      character*12 sp_div_opdesc  
      character*12 sp_equal_opdesc 
      character*12 sp_nequal_opdesc
      character*12 sp_ge_opdesc       
      character*12 sp_le_opdesc  
      character*12 sp_gt_opdesc 
      character*12 sp_lt_opdesc   
      character*12 sp_ldi_opdesc
      character*12 sp_ldindex_opdesc
      character*12 pardo_opdesc
      character*12 endpardo_opdesc
      character*12 exit_opdesc
      character*12 assignment_opdesc
      character*12 self_multiply_opdesc 
      character*12 subtract_opdesc  
      character*12 collective_sum_opdesc  
      character*12 divide_opdesc  
      character*12 prepare_opdesc
      character*12 request_opdesc
      character*12 prequest_opdesc
      character*12 compute_integrals_opdesc
      character*12 put_replace_opdesc
      character*12 tensor_opdesc
      character*12 prepare_increment_opdesc
      character*12 allocate_opdesc
      character*12 deallocate_opdesc
      character*12 sp_ldi_sym_op_opdesc
      character*12 create_window_op_opdesc
      character*12 where_op_opdesc

      parameter (contraction_opdesc = 'contraction')
      parameter ( sum_opdesc        = 'sum')
      parameter ( reindex_opdesc    = ' ')
      parameter ( do_opdesc         = ' ')
      parameter ( enddo_opdesc      = ' ')
      parameter ( get_opdesc        = 'get')
      parameter ( user_sub_opdesc   = 'user_sub')
      parameter ( put_opdesc        = 'put')
      parameter ( go_to_opdesc      = ' ')
      parameter ( create_opdesc     = 'create')
      parameter ( delete_opdesc     = 'delete')
      parameter ( call_opdesc       = ' ')
      parameter ( return_opdesc     = ' ')
      parameter ( jz_opdesc         = ' ')
      parameter ( stop_opdesc       = ' ')
      parameter ( sp_add_opdesc     = ' ')
      parameter (sp_sub_opdesc      = ' ')
      parameter ( sp_mult_opdesc    = ' ')
      parameter (sp_div_opdesc      = ' ')
      parameter ( sp_equal_opdesc   = ' ')
      parameter ( sp_nequal_opdesc  = ' ')
      parameter ( sp_ge_opdesc      = ' ') 
      parameter ( sp_le_opdesc      = ' ')
      parameter ( sp_gt_opdesc      = ' ')
      parameter ( sp_lt_opdesc      = ' ')
      parameter ( sp_ldi_opdesc     = ' ')
      parameter ( sp_ldindex_opdesc = ' ')
      parameter ( pardo_opdesc      = ' ')
      parameter ( endpardo_opdesc   = ' ')
      parameter ( exit_opdesc       = ' ')
      parameter ( assignment_opdesc = 'assignment')
      parameter ( self_multiply_opdesc = 'multiply')
      parameter ( subtract_opdesc   = 'subtract')
      parameter ( collective_sum_opdesc = 'collective') 
      parameter ( divide_opdesc     = 'divide')
      parameter ( prepare_opdesc    = 'prepare')
      parameter ( request_opdesc    = 'request')
      parameter ( prequest_opdesc   = 'prequest')
      parameter ( compute_integrals_opdesc = 'integrals')
      parameter ( put_replace_opdesc = 'put_replace')
      parameter ( tensor_opdesc      = 'tensor mult.')
      parameter ( prepare_increment_opdesc = 'prepare_inc')
      parameter ( allocate_opdesc    = 'allocate')
      parameter ( deallocate_opdesc  = 'deallocate')
      parameter ( sp_ldi_sym_op_opdesc = ' ')
      parameter ( create_window_op_opdesc = 'create_win')
      parameter ( where_op_opdesc    = ' ')

      integer served_array
      integer local_array
      integer static_array
      integer distributed_array
      integer temp_array
      integer scalar_value
      integer dummy_array_type
      parameter (served_array    = 201)
      parameter (static_array       = 202) 
      parameter (distributed_array = 203)
      parameter (temp_array        = 204)
      parameter (scalar_value      = 205)
      parameter (local_array       = 206)
      parameter (dummy_array_type  = 207)

      integer write_only_array_status
      integer read_only_array_status
      parameter (write_only_array_status = -1)
      parameter (read_only_array_status  = -2)

c--------------------------------------------------------------------------
c   Symbolic constants used to define index ranges.
c--------------------------------------------------------------------------

      integer norb
      integer nocc
      integer nvirt
      integer bocc
      integer eocc
      integer bvirt
      integer evirt
      integer naocc
      integer nbocc
      integer navirt
      integer nbvirt
      integer itrips_parm
      integer itripe_parm
      integer ihess1_parm
      integer ihess2_parm
      integer jhess1_parm
      integer jhess2_parm
      integer subb_parm
      integer sube_parm
      integer baocc
      integer bbocc
      integer eaocc
      integer ebocc
      integer bavirt
      integer bbvirt
      integer eavirt
      integer ebvirt
      integer noccorb, nvirtorb, boccorb,eoccorb,bvirtorb, evirtorb,
     *  naoccorb, nboccorb, navirtorb, nbvirtorb, baoccorb, bboccorb,
     *  eaoccorb, eboccorb, bavirtorb, bbvirtorb, eavirtorb, ebvirtorb
      integer cc_iter_cons
      integer cc_hist_cons
      integer cc_beg_cons
      integer scf_iter_cons
      integer scf_hist_cons
      integer scf_beg_cons
      integer natoms_cons
      integer aoindex
      integer moindex
      integer moaindex
      integer mobindex
      integer simple_index
      integer laindex
      integer subindex
      integer undefined_segment
      integer sip_sub_segsize_parm
      integer sip_sub_occ_segsize_parm
      integer sip_sub_virt_segsize_parm
      integer sip_sub_ao_segsize_parm
      parameter ( norb = 0)
      parameter ( nocc = -1)
      parameter ( nvirt = -2)
      parameter ( bocc = -3)
      parameter ( eocc = -4)
      parameter ( bvirt = -5)
      parameter ( evirt = -6)
      parameter ( naocc = -7)
      parameter ( nbocc = -8)
      parameter ( navirt = -9)
      parameter ( nbvirt = -10)
      parameter ( baocc = -11)
      parameter ( bbocc = -12)
      parameter ( eaocc = -13)
      parameter ( ebocc = -14)
      parameter ( bavirt = -15)
      parameter ( bbvirt = -16)
      parameter ( eavirt = -17)
      parameter ( ebvirt = -18)
      parameter ( noccorb = -19)
      parameter ( nvirtorb = -20)
      parameter ( boccorb  = -21)
      parameter ( eoccorb  = -22)
      parameter ( bvirtorb = -23) 
      parameter ( evirtorb = -24)
      parameter ( naoccorb = -25) 
      parameter ( nboccorb = -26) 
      parameter ( navirtorb = - 27)
      parameter ( nbvirtorb = -28)
      parameter ( baoccorb  = -29)
      parameter ( bboccorb  = -30)
      parameter ( eaoccorb  = -31)
      parameter ( eboccorb  = -32)
      parameter ( bavirtorb = -33)
      parameter ( bbvirtorb = -34)
      parameter ( eavirtorb = -35)
      parameter ( ebvirtorb = -36)
      parameter (cc_iter_cons  = -37)
      parameter (cc_hist_cons  = -38)
      parameter (cc_beg_cons   = -39)
      parameter (scf_iter_cons = -40)
      parameter (scf_hist_cons = -41)
      parameter (scf_beg_cons  = -42)
      parameter (natoms_cons   = -43)
      parameter (itrips_parm   = -44)
      parameter (itripe_parm   = -45)
      parameter (ihess1_parm   = -46)
      parameter (ihess2_parm   = -47)
      parameter (jhess1_parm   = -48)
      parameter (jhess2_parm   = -49)
      parameter (subb_parm     = -50)
      parameter (sube_parm     = -51)
      parameter (sip_sub_segsize_parm = -52)
      parameter (sip_sub_occ_segsize_parm = -53)
      parameter (sip_sub_virt_segsize_parm = -54)
      parameter (sip_sub_ao_segsize_parm = -55)

      parameter ( aoindex = 1001)
      parameter ( moindex = 1002)
      parameter ( moaindex = 1003)
      parameter ( mobindex = 1004) 
      parameter ( simple_index = 1005) 
      parameter ( laindex  = 1006)
      parameter ( subindex  = 1007)
      parameter ( undefined_segment = -90909)

      integer wildcard_indicator
      parameter (wildcard_indicator = 90909)

c-------------------------------------------------------------------------
c   Message tags for distributing SIP tables to processors.
c-------------------------------------------------------------------------

      integer sip_table_header_tag
      parameter (sip_table_header_tag = 7007)
      integer sip_index_table_send_tag
      parameter (sip_index_table_send_tag = 7008)
      integer sip_array_table_send_tag
      parameter (sip_array_table_send_tag = 7009)
      integer sip_optable_send_tag
      parameter (sip_optable_send_tag = 7010)
      integer sip_scalar_table_send_tag
      parameter (sip_scalar_table_send_tag = 7011)

c--------------------------------------------------------------------------
c   SIP Server parameters.
c--------------------------------------------------------------------------

      integer sip_server_prepare
      integer sip_server_request
      integer sip_server_data_message
      integer sip_server_message
      integer sip_server_response_message
      integer sip_mgr_ready_tag 
      integer sip_server_quit
      integer sip_server_prepare_increment
      integer len_sip_server_message
      integer len_sip_server_prequest
      integer sip_server_barrier_signal
      integer sip_server_stat_data_tag
      integer sip_server_copy_message
      integer sip_server_delete_message
      integer sip_server_prequest
      integer sip_server_blocks_to_list
      integer sip_server_checkpoint_signal
      integer sip_server_restart_signal
      integer sip_server_commit_signal
      integer sip_server_list_to_blocks
      parameter (sip_server_prepare           = 3333)     ! "prepare =" syntax
      parameter (sip_server_request           = 3334)
      parameter (sip_server_data_message      = 3335)
      parameter (sip_server_message           = 3336)
      parameter (sip_server_response_message  = 3337)
      parameter (sip_server_quit              = 3338)
      parameter (sip_mgr_ready_tag            = 3339)
      parameter (sip_server_prepare_increment = 3340)   ! "prepare +=" syntax
      parameter (sip_server_barrier_signal    = 3341)
      parameter (sip_server_stat_data_tag     = 3342)
      parameter (sip_server_copy_message      = 3343)
      parameter (sip_server_delete_message    = 3344)
      parameter (sip_server_prequest          = 3345)
      parameter (sip_server_blocks_to_list    = 3346)
      parameter (sip_server_checkpoint_signal = 3347)
      parameter (sip_server_restart_signal    = 3348)
      parameter (sip_server_commit_signal     = 3349)
      parameter (sip_server_list_to_blocks    = 3350)
      parameter (len_sip_server_message = 6+3*mx_array_index)
      parameter (len_sip_server_prequest = len_sip_server_message+
     *                                       3*mx_array_index)

c---------------------------------------------------------------------------
c   Message tags for load-balancing pardo.
c----------------------------------------------------------------------------

      integer pardo_job_tag
      integer pardo_job_request_tag
      integer thread_server_ready_tag
      parameter (pardo_job_tag         = 4449)
      parameter (pardo_job_request_tag = 4450)
      parameter (thread_server_ready_tag = 3456)

c--------------------------------------------------------------------------
c   Message tags for writing blocks to lists.
c--------------------------------------------------------------------------

      integer block_list_descriptor_tag
      parameter (block_list_descriptor_tag = 4446)
      integer block_list_array_tag
      parameter (block_list_array_tag = 4451)

