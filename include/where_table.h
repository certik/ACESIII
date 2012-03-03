      integer iwhere
      integer nwhere
      integer max_where_table
      parameter (max_where_table = 500)
      integer c_where_ind1
      parameter (c_where_ind1 = 1)
      integer c_where_cond
      parameter (c_where_cond = 2)
      integer c_where_ind2
      parameter (c_where_ind2 = 3)
      integer where_table

      common /where_table_block/iwhere, nwhere, 
     *                          where_table(max_where_table,3)
