      common /sheap_block/ishptr, i8shptr, dshptr

      integer ishared_heap(*)
      pointer (ishptr, ishared_heap)
      integer i8shared_heap(*)
      pointer (i8shptr, i8shared_heap)
      double precision dshared_heap(*)
      pointer (dshptr, dshared_heap)
