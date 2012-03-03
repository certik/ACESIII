      integer max_mem_handle
      integer max_config_handle
      integer memtype
      integer memflag
      integer memsize
      integer*8 max_memsize
      integer mem_context
      integer config_table
      integer*8 memtable
      integer anchor(1)
      integer*8 mem_index
      integer*8 perm_offset, temp_offset, temp_ptr
      integer perm_memory, temp_memory

      parameter (perm_memory = 1)
      parameter (temp_memory = 2)
      parameter (max_mem_handle = 200)
      parameter (max_config_handle = 50)
      parameter (memsize = 1)
      parameter (memtype = 2)
      parameter (memflag = 3)
      common /memtable/mem_index, temp_ptr, max_memsize,
     *                 perm_offset(max_mem_handle),
     *                 temp_offset(max_mem_handle), 
     *                 memtable(max_mem_handle,memflag), anchor,
     *                 config_table(max_config_handle,max_mem_handle),
     *                 mem_context

