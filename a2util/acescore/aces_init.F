
c This routine initializes the ACES environment for use by any program.
c It or aces_init_rte should be called before any other ACES routine.

c INPUT
c logical BALLOCMEM : do (not) allocate memory according to IFLAGS(h_IFLAGS_mem)
c                     NOTE: If memory is allocated, then aces_init will
c                           initialize the I/O subsystem and the cache
c                           subsystem.

c INPUT/OUTPUT
c int ICORE(*) : a reference integer (address) used to anchor a malloc'ed heap,
c                the value of which is set to the byte address of the heap
c                (not valid on 32-bit builds with 64-bit pointers)

c OUTPUT
c int ICORENDX : the integer index that addresses the malloc'ed heap
c int ICOREDIM : the number of addressable integers at ICORE(ICORENDX)
c int IUHF : the spin state of the reference wavefunction
c            = 0; spin-  restricted (RHF)
c            = 1; spin-unrestricted (UHF)

#include "f_types.h"
c#define _DEBUG_ACES_INIT

      subroutine aces_init(iCore,iCoreNdx,iCoreDim,iUHF,bAllocMem)
      implicit none

c ARGUMENTS
      integer iCore(*), iCoreNdx, iCoreDim, iUHF
      logical bAllocMem

c EXTERNAL FUNCTIONS
      F_ADR c_adr
      external c_adr

c PARAMETERS
      integer iMemMin, iMemInc
      parameter (iMemMin=1*1024*1024,iMemInc=1*1024*1024)

c INTERNAL VARIABLES
      double precision dTmp
      integer i0, iMem, iTmp

c TIMING VARIABLES
      integer stime_sec, stime_usec
      integer utime_sec, utime_usec
      integer rtime_sec, rtime_usec

c COMMON BLOCKS
#include "flags.com" /* for iflags */
#include "flags2.com" /* for iflags2 */
#include "machsp.com" /* for iintfp */
#include "aces_time.com" /* for timing stats */
c BWCC needs aces_init to be reentrant, so an initialized common value must
c exist. The blockdata routine is appended to the end of this file.
      external aces_bd_aces_reflag
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/

c ----------------------------------------------------------------------

#ifdef _TRAPFPE
c   o enable floating-point exception handling if the compiler cannot
      call trapfpe
#endif

c   o initialize timing information
      call c_rutimes(utime_sec,utime_usec,stime_sec,stime_usec)
      call c_gtod   (rtime_sec,rtime_usec)
      ame_stime_in = (1.d-6 * stime_usec) + stime_sec
      ame_utime_in = (1.d-6 * utime_usec) + utime_sec
      ame_rtime_in = (1.d-6 * rtime_usec) + rtime_sec
      ame_timed = .true.

#ifdef _DEBUG_ACES_INIT
      print *, '@ACES_INIT: (start) initializing the ACES environment'
#endif

c   o setup stdout buffering to see output immediately
      call bufferlines

c   o initialize the runtime environment (held in common blocks)
      call aces_init_rte

c   o gather parallel statistics (needed by gfname in ja_init)
      call aces_com_parallel_aces

c   o initialize the job archive subsystem (getrec/putrec)
      call aces_ja_init
c???      call getrec(1,'JOBARC','IENDSTAT',1,iStat)

c   o load the ACES State Variables
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)

c   o allocate memory and initialize the I/O subsystem
      if (bAllocMem) then

c      o load the requested core size
         iMem = iflags(36)

c      o allocate core memory
         if (xFlag.eq.0) then
            iCore(1) = 0
            do while ((iCore(1).eq.0).and.(iMem.gt.iMemMin))
               call aces_malloc(iMem,iCore,i0)
               if (iCore(1).eq.0) iMem = iMem - iMemInc
            end do
#ifdef _DEBUG_ACES_INIT
            print *, '@ACES_INIT: heap size   = ',iMem,' integers'
            print *, '            iCore index = ',i0
#endif
            if (iMem.lt.iflags(36)) then
               print *, '@ACES_INIT: MEMORY WARNING!'
               print *, '            requested ',iflags(36),' integers'
               print *, '            allocated ',iMem,' integers'
            end if
            if (iCore(1).ne.0) then
#ifdef _DEBUG
c            o make sure iCore(i0) starts on a double word boundary
               print *, '@ACES_INIT: heap range',c_adr(iCore(i0)),
     &                                           c_adr(iCore(i0+iMem))-1
#ifdef _PTRS_ARE_WORDS
               iTmp = iand(c_adr(iCore(i0)),iintfp-1)
               if (iTmp.ne.0) then
                  print *, '@ACES_INIT: alignment defect in aces_malloc'
                  print *, '            absolute offset = ',iTmp,' W'
                  call aces_exit(1)
               end if
#else
               iTmp = iand(c_adr(iCore(i0)),ifltln-1)
               if (iTmp.ne.0) then
                  print *, '@ACES_INIT: alignment defect in aces_malloc'
                  print *, '            absolute offset = ',iTmp,' B'
                  call aces_exit(1)
               end if
#endif
#endif
               xFlag = 1
            else
c               print *, '@ACES_INIT: Request for ',iMem,
c     &                  ' integers of memory failed.'
               print *, '@ACES_INIT: unable to allocate at least ',
     &                  iMemMin,' integers of memory'
               call aces_exit(1)
            end if
         end if

c      o initialize the I/O subsystem ('T' creates and initializes the cache)
         call aces_io_init(iCore,i0,iMem,.true.)

c      o transfer the iCore statistics
         iCoreNdx = i0
         iCoreDim = iMem

c     else if (.not.bAllocMem)
      else

         iCore(1) = 0
         iCoreNdx = 1
         iCoreDim = 1

c     end if (bAllocMem)
      end if

c   o assign iUHF
      if (iflags(11).eq.0) then
         iUHF = 0
      else
         iUHF = 1
      end if

c   o initialize the chemical system (held in common blocks)
      call aces_init_chemsys

#ifdef _DEBUG_ACES_INIT
      print *, '@ACES_INIT: (done)  initializing the ACES environment'
#endif

      return
c     end subroutine aces_init
      end

c ----------------------------------------------------------------------
      blockdata aces_bd_aces_reflag
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/
      data                 xFlag /0/
      end
c ----------------------------------------------------------------------

