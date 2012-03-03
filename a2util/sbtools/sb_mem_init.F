
c This routine carves up the icore(i0) space for using setptr/relptr.
c It's clumsy, but you must communicate with this routine through
c sbcore.com by setting ineeded, dneeded, maxicore, etc.

c For dummy runs (when memknown=0), iSize sets maxmem.

c INPUT:
c int iHeap : the first integer in a heap (i.e. icore(i0) in ACES lingo)
c int iSize : the number of INTEGERS in the heap

      subroutine sb_mem_init(iHeap,iSize)
      implicit none

#include "f_types.h"
#include "aces.h"
#include "machsp.com"
#include "sbcore.com"
#include "sb_mem.com"

      integer iHeap(*), iSize
      F_ADR c_adr, zheap, ztmp
      external c_adr

      call callstack_push('SB_MEM_INIT')

#ifdef _ASSERT
      if (iSize.le.0) then
         print *, '@SB_MEM_INIT: Assertion failed.'
         print *, '              iSize = ',iSize
         call errex
      end if
#endif

      if (memknown.ne.0) then

         if (dynmem.eq.1) then
            maxmem=iSize
            if (memknown.eq.-1) then
               ineeded=maxicore
               dneeded=(maxmem-maxicore)/iintfp
            end if
            zheap = c_adr(iHeap)

c i0 is the icore index that addresses iHeap
#ifdef _PTRS_ARE_WORDS
            ztmp = (zheap-c_adr(icore(1)))
#else
c shift zheap up to a 64-bit boundary
            zheap = ifltln * ((zheap+ifltln-1)/ifltln)
            ztmp = (zheap-c_adr(icore(1)))/iintln
#endif /* _PTRS_ARE_WORDS */
            i0 = 1 + ztmp
            if (i0.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

c i1 is the first icore index that cannot be used (after the "icore" heap)
c (Yau: shift this up to 64-bit boundary)
            ztmp = ztmp + ineeded + iand(ineeded,iintfp-1)
            i1 = 1 + ztmp
            if (i1.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

c repeat the process for dcore
            zheap = c_adr(icore(i1))

c d0 is the dcore index that addresses zheap
#ifdef _PTRS_ARE_WORDS
            ztmp = (zheap-c_adr(dcore(1)))
#else
            ztmp = (zheap-c_adr(dcore(1)))/ifltln
#endif /* _PTRS_ARE_WORDS */
            d0 = 1 + ztmp
            if (d0.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

c d1 is the first dcore index that cannot be used (after the "dcore" heap)
            ztmp = ztmp + dneeded
            d1 = 1 + ztmp
            if (d1.ne.1+ztmp) then
               print *, '@SB_MEM_INIT: Integer overflow.'
               print *, '              The heap cannot be addressed.'
               call errex
            end if

         else
            maxmem=0
            maxicore=nondynimem
            ineeded=nondynimem
            dneeded=nondyndmem
            i0=1
            d0=1
            i1=nondynimem+1
            d1=nondyndmem+1
         end if
         iptr=i0
         dptr=d0
         call izero(icore(i0),ineeded)
         call dzero(dcore(d0),dneeded)

      else ! dummy run

         if (dynmem.eq.1) then
            maxmem=iSize
            i1=maxicore+1
            d1=(maxmem-maxicore)/iintfp+1
         else
            maxmem=0
            maxicore=nondynimem
            i1=nondynimem+1
            d1=nondyndmem+1
         end if
         i0=1
         d0=1
         iptr=1
         dptr=1
         ineeded=0
         dneeded=0

      end if

      call callstack_pop
      return
      end

