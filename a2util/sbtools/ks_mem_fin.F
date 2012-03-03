
c This routine summarizes the memory usage.

      subroutine sb_mem_fin
      implicit none

#include "aces.h"
#include "machsp.com"
#include "kscore.com"
#include "ks_mem.com"

      call callstack_push('SB_MEM_FIN')

      if (memknown.ne.1) then
#ifdef _DEBUG
         if (dynmem.eq.1) then
            write(*,9000) ineeded
            write(*,9010) dneeded*iintfp
            write(*,9020) ineeded+dneeded*iintfp
            write(*,9030) maxmem
         else
            write(*,9000) ineeded
            write(*,9010) dneeded
            write(*,9040) i1-1
            write(*,9050) d1-1
         end if
#endif
      end if
#ifdef _DEBUG
 9000 format(t3,'@SB_MEM, You need ',i9,' words of kscore memory.')
 9010 format(t3,'       , You need ',i9,' words of dcore memory.')
 9020 format(t3,'       , You need ',i9,' total words of memory.')
 9030 format(t3,'       , You have ',i9,' words of memory.')
 9040 format(t3,'       , You have ',i9,' words of kscore memory.')
 9050 format(t3,'       , You have ',i9,' words of dcore memory.')
#endif

      if (memknown.ne.0) then
c        Free memory if dynmem
      else
         if (dynmem.eq.1) then
            if (ineeded+dneeded*iintfp .gt. maxmem) then
               write(*,*) '@SB_MEM_FIN: exceeded available memory'
               call errex
               stop
            end if
         else
            if (ineeded.gt.i1-1 .or. dneeded.gt.d1-1) then
               write(*,*) '@SB_MEM_FIN: exceeded available memory'
               call errex
               stop
            end if
         end if
      end if

      call callstack_pop
      return
      end

