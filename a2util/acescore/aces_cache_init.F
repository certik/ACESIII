
c This routine reserves space at iCore(i0) for an I/O cache.

c The address of iCore(i0) is saved before returning; therefore, there is no
c significance placed on these two variables besides pointing to the usable
c memory heap.

c INPUT
c int ICORE(*) : the memory anchor used to place the cache

c INTPUT/OUTPUT
c int I0 : on input, this is the ICORE index that places the cache
c        : on output, this points to the first integer after the cache
c int IMEM : the amount of usable memory at ICORE(I0)
c int NPRCWD : the desired/actual integer-length of a physical record
c int NSLOTS : the desired/actual number of cache slots

#include "f_types.h"
c#define _DEBUG_ACES_CACHE_INIT

      subroutine aces_cache_init(iCore,i0,iMem,nPRcWd,nSlots)
      implicit none

c ARGUMENTS
      integer iCore(*), i0, iMem, nPRcWd, nSlots

c EXTERNAL FUNCTIONS
      F_ADR c_adr
      external c_adr

c INTERNAL VARIABLES
      double precision dTmp
      F_ADR zTmp
      integer i, ndx, iTmp

c COMMON BLOCKS
#include "cache.com"
#include "filspc.com"
#include "machsp.com" /* for iintln and ifltln */
#include "icdacc.com" /* for idaccm */

c ----------------------------------------------------------------------

c   o the cache subsystem never turns off
      if (bCacheUp) return

c   o turn on the cache subsystem flag
      bCacheUp = .true.

c   o make sure iMem is natural
      if (iMem.lt.1) then
         print *, '@ACES_CACHE_INIT: There is no place for the cache.'
         print *, '   iMem = ',iMem
         call aces_exit(1)
      end if

c   o initialize the filspc common block
c     (set the length of a physical record and cache slot)
      if (nPRcWd.lt.1) then
         print *, '@ACES_CACHE_INIT: The length of a cache slot must ',
     &            'be at least 1 integer long.'
         print *, '                  requested length = ',nPRcWd
         call aces_exit(1)
      end if
      iprcln = nPRcWd*iintln
      iprcwd = nPRcWd
      iprcrl = nPRcWd*idaccm
      if (iand(iprcln,ifltln-1).ne.0) then
         print '(/)'
         print *, '@ACES_CACHE_INIT: WARNING - The length of a cache ',
     &            'slot is not a whole multiple'
         print *, '                  of floating-point numbers. This ',
     &            'will negatively impact'
         print *, '                  performance.'
         print '(/)'
      end if

c   o set the number of cache slots
      if (nSlots.lt.1) then
         print *, '@ACES_CACHE_INIT: There must be at least 1 cache ',
     &            'slot.'
         print *, '                  requested number = ',nSlots
         call aces_exit(1)
      end if
      if (nSlots.gt._MAX_CACHE_SLOTS) then
         print *, '@ACES_CACHE_INIT: There can be at most ',
     &            _MAX_CACHE_SLOTS,' cache slots.'
         print *, '                  requested number = ',nSlots
         print *, '                  resetting to ',_MAX_CACHE_SLOTS
         nSlots = _MAX_CACHE_SLOTS
      end if
      cachnum = nSlots

c   o use ndx to point to the first double boundary in iCore(i0)
#ifdef _PTRS_ARE_WORDS
      zTmp = 1 + c_adr(iCore(i0)) - c_adr(icache(1))
#else
      zTmp = c_adr(iCore(i0))
      if (iand(zTmp,ifltln-1).ne.0) then
         i0 = i0 + 1
         zTmp = c_adr(iCore(i0))
         if (iand(zTmp,ifltln-1).ne.0) then
            print *, '@ACES_CACHE_INIT: The core memory is not aligned.'
            print *, '   &icore[0] = ',c_adr(iCore(1))
            call aces_exit(1)
         end if
         iMem = iMem - 1
      end if
      zTmp = zTmp - c_adr(icache(1))
      if (iand(zTmp,iintln-1).ne.0) then
         print *, '@ACES_CACHE_INIT: The cache anchor is not aligned.'
         print *, '   &icache[0] = ',c_adr(icache(1))
         call aces_exit(1)
      end if
      zTmp = 1+(zTmp/iintln)
#endif /* _PTRS_ARE_WORDS */
      ndx = zTmp
      if (ndx.ne.zTmp) then
         print *, '@ACES_CACHE_INIT: Assertion failed.'
         print *, '                  The I/O cache is not addressable.'
         call aces_exit(1)
      end if
      do i = 1, cachnum
         cachndx(i) = ndx
         ndx = ndx + iprcwd
      end do
#ifdef _DEBUG
      print *, '@ACES_CACHE_INIT: allocated >=',
     &         ishft(iprcln*cachnum,-20),
     &         ' MB of cache from core memory'
#endif

c   o move i0 to a double boundary and decrease iMem
      i = iprcwd * cachnum
      i0 = i0 + i
      zTmp = c_adr(iCore(i0))
#ifdef _PTRS_ARE_WORDS
      if (iand(zTmp,iintfp-1).ne.0) then
         i0 = i0 + 1
         iMem = iMem - 1
      end if
#else
      if (iand(zTmp,ifltln-1).ne.0) then
         i0 = i0 + 1
         iMem = iMem - 1
      end if
#endif /* _PTRS_ARE_WORDS */
      iMem = iMem - i
      if (iMem.lt.1) then
         print *, '@ACES_CACHE_INIT: The total cache size must be less',
     &            ' than the total core size.'
         print *, '                  cache slots = ',cachnum
         print *, '                  slot size   = ',iprcwd,' integers'
         print *, '                  core size   = ',iMem,' integers'
         call aces_exit(1)
      end if

c   o (re)set the cache statistics
      call aces_cache_reset

#ifdef _DEBUG_ACES_CACHE_INIT
      print *, '@ACES_CACHE_INIT: /filspc/ iprcln = ',iprcln,
     &                                  ', iprcwd = ',iprcwd,
     &                                  ', iprcrl = ',iprcrl
      print *, '                  /cache/ cachnum = ',cachnum,
     &                                 ', cachndx = ',cachndx
#endif

      return
c     end subroutine aces_cache_init
      end

