#ifndef _RIC_HEAP_COM_
#define _RIC_HEAP_COM_
c ric_heap.com : begin

c This common block contains the heap address and array indices for
c processing RICs. New additions must be initialized to 1 in
c bd_ric_heap.F and set in init_ric_heap.F.

#ifndef NO_EXTERNAL
      external bd_ric_heap
#endif

      double precision dRICHeap(1)
      integer*8 z_RICHeap, z_DerBMat, z_BMat, z_GMat, z_BTGInv
      integer*8 heapptr

      common /ric_heap_com/ dRICHeap,
     &                      z_RICHeap, z_DerBMat, z_BMat, z_GMat,
     &                      z_BTGInv, heapptr
      save   /ric_heap_com/

c ric_heap.com : end
#endif /* _RIC_HEAP_COM_ */
