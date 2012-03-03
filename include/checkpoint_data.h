
      integer mx_allocates
      parameter (mx_allocates = 1000)
      integer mx_creates
      parameter (mx_creates = 1000)
      integer mx_ckpt_arrays 
      parameter (mx_ckpt_arrays = 1000)
      integer nactive_allocate_table, nactive_create_table,
     *        master_ckpt_unit     
      integer active_allocate_table, active_allocate_op
      integer active_create_table, active_create_op
      integer ckpt_arrays
      integer*8 ckpt_diskaddr, disk_sizes
      integer*8 free_space, free_space_sizes
      integer*8 free_space_candidates, candidate_sizes
      integer nckpt_arrays
      integer nfree_space, nfree_space_candidates
      logical restart_job, restart_status

      common /checkpoint_data/active_allocate_table(mx_allocates),
     *                        active_allocate_op(mx_allocates),
     *                        active_create_table(mx_creates),
     *                        active_create_op(mx_creates),
     *                        nactive_allocate_table,
     *                        nactive_create_table,
     *                        master_ckpt_unit,
     *                        ckpt_arrays(mx_ckpt_arrays),
     *                        nckpt_arrays, nfree_space, 
     *                        nfree_space_candidates, restart_job,
     *                        restart_status
      common /ckpt_diskloc/ckpt_diskaddr(mx_ckpt_arrays),
     *                     disk_sizes(mx_ckpt_arrays),
     *                     free_space(mx_ckpt_arrays),
     *                     free_space_sizes(mx_ckpt_arrays), 
     *                     free_space_candidates(mx_ckpt_arrays), 
     *                     candidate_sizes(mx_ckpt_arrays)
