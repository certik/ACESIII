      integer max_array_list
      parameter (max_array_list = 100)
      integer array_list
      integer narray_list
      common /blockdata_list/array_list(max_array_list), narray_list

c---------------------------------------------------------------------------
c   Signal definitions.
c---------------------------------------------------------------------------

      integer blocks_list_done
      parameter (blocks_list_done = 90909)
      integer list_to_blocks_done
      parameter (list_to_blocks_done = 80808)

