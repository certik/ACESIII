      integer me
      integer nprocs
      integer my_company_rank, my_company_size
      integer restart_file
      logical mpi_io_support

      common /parallel_info/me, nprocs, my_company_rank, 
     *                      my_company_size, restart_file,
     *                      mpi_io_support
