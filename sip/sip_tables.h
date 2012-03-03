
      integer*8 optable_base
      integer*8 index_table_base
      integer*8 eindex_table_base
      integer*8 array_table_base
      integer*8 scalar_table_base
      integer*8 address_table_base
      integer*8 segment_table_base
      integer*8 bmap_table_base
      integer*8 t1_handle_base
      integer*8 proctab_base
      integer noptable_sip
      integer narray_table_sip
      integer nindex_table_sip
      integer nscalar_table_sip
      integer nsegment_table_sip 
      integer mx_noptable
      integer mx_nindex_table 
      integer mx_narray_table
      integer mx_scalar_table
      common /sip_tables/optable_base, index_table_base, 
     *                   eindex_table_base, 
     *                   array_table_base, scalar_table_base,
     *                   noptable_sip, 
     *                   nindex_table_sip,
     *                   narray_table_sip, nscalar_table_sip,
     *                   nsegment_table_sip,
     *                   mx_noptable, 
     *                   mx_nindex_table,
     *                   mx_narray_table , mx_scalar_table
      common /runtime_tables/bmap_table_base, segment_table_base,
     *                   proctab_base, address_table_base, 
     *                   t1_handle_base
