c-----------------------------------------------------------------------
c   Defines the columns of the PST (Process Status Table).
c-----------------------------------------------------------------------
      integer r_role
      integer r_company_id
      integer r_company_comm
      integer r_company_rank
      integer r_last

      parameter (r_role = 1)
      parameter (r_company_id = 2)
      parameter (r_company_comm = 3) 
      parameter (r_company_rank  = 4)
      parameter (r_last = r_company_rank)

      integer maxprocs
      integer pst
      parameter (maxprocs = 10000)
      common /pst/pst(maxprocs, r_last)

