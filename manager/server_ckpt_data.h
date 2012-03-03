      integer mx_ckpt_arrays
      parameter (mx_ckpt_arrays=1000)
      integer ckpt_dat_unit, ckpt_ndx_unit

      common /server_ckpt_data/ckpt_dat_unit, ckpt_ndx_unit

      integer index_file_header_size
      parameter (index_file_header_size = 4*mx_ckpt_arrays)

      integer*8 ndx_file_header(4,mx_ckpt_arrays)
      integer*8 data_file_sizes(mx_ckpt_arrays)
      integer*8 data_file_location(mx_ckpt_arrays)
      logical commit_flag
      common /ndx_file/ndx_file_header, data_file_sizes, 
     *                 data_file_location,commit_flag

      character*20 ckpt_ndx_filename
      character*20 ckpt_dat_filename
      common /server_ckpt_filenames/ckpt_ndx_filename,ckpt_dat_filename

      integer nfree_space_ndx, nfree_space_dat, nfree_ndx_can, 
     *        nfree_dat_can
      integer*8 free_space_ndx(2,mx_ckpt_arrays), 
     *          free_space_dat(2,mx_ckpt_arrays),
     *          free_space_ndx_candidate(2,mx_ckpt_arrays), 
     *          free_space_dat_candidate(2,mx_ckpt_arrays)

      common /server_ckpt_free_space/nfree_space_ndx, nfree_space_dat,
     *             nfree_ndx_can, nfree_dat_can,
     *             free_space_ndx, free_space_dat,
     *             free_space_ndx_candidate, free_space_dat_candidate 
