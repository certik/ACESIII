c-------------------------------------------------------------------------
c
c   Table containing company configuration data.
c
c-------------------------------------------------------------------------

      integer max_company
      integer lcompany_entry
      integer inter_comm, inter_group
      parameter (max_company = 300)
      parameter (lcompany_entry = 7)

      integer c_company_id 
      parameter (c_company_id = 1)
      integer c_platoon_id
      parameter (c_platoon_id = 2)
      integer c_nmgr
      parameter (c_nmgr = 3)
      integer c_mgr_mem 
      parameter (c_mgr_mem = 4)
      integer c_nwrkr
      parameter (c_nwrkr = 5)
      integer c_wrkr_mem
      parameter (c_wrkr_mem = 6)
      integer c_iocompany
      parameter (c_iocompany = 7)

      common /company/c_table(max_company, lcompany_entry), inter_comm,
     *                inter_group, company_sial_prog(max_company)
      integer c_table
      character*80 company_sial_prog
