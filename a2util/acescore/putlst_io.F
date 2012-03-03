
c This routine stores physical records to list files.
c It is a primitive of putlst and should not be called directly.

      subroutine putlst_io(ZList,iFile,iRecNum,iWrdNdx,nWords)
      implicit none

c ARGUMENTS
      integer ZList(*), iFile, iRecNum, iWrdNdx, nWords

c EXTERNAL FUNCTIONS
      integer iieq, iiamin, iiamax

c INTERNAL VARIABLES
      integer iFileOut, iRecOut
      integer iTmp, iLoc, iNdx, iVictim
      logical bNotDone

c COMMON BLOCKS
#include "cache.com"
#include "filspc.com" /* for iprcwd */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert cache subsystem is up
      if (.not.bCacheUp) then
         print *, '@PUTLST_IO: Assertion failed.'
         print *, '   bCacheUp = ',bCacheUp
         iTmp = 1
      end if
c   o assert iFile, iRecNum, and iWrdNdx are all natural
      if ((iFile.lt.1).or.(iRecNum.lt.1).or.(iWrdNdx.lt.1)) then
         print *, '@PUTLST_IO: Assertion failed.'
         print *, '   iFile   = ',iFile
         print *, '   iRecNum = ',iRecNum
         print *, '   iWrdNdx = ',iWrdNdx
         iTmp = 1
      end if
c   o assert nWords >= 0
      if (nWords.lt.0) then
         print *, '@PUTLST_IO: Assertion failed.'
         print *, '   nWords = ',nWords
         iTmp = 1
      end if
c   o assert maximum element is in physical record
      if ((iWrdNdx-1+nWords).gt.iprcwd) then
         print *, '@PUTLST_IO: Assertion failed.'
         print *, '   iWrdNdx = ',iWrdNdx
         print *, '   nWords  = ',nWords
         print *, '   iprcwd  = ',iprcwd
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif /* _ASSERT */

      if ((iWrdNdx.lt.1).or.(nWords.lt.1)) return

c ----------------------------------------------------------------------

c   o increment the cachetime
      cachetime = cachetime + 1
      if (cachetime.lt.0) then
c      o scale-back the lrustats
         iLoc = iiamin(cachnum,lrustats,1)
         iTmp = lrustats(iLoc)
         do iLoc = 1, cachnum
            lrustats(iLoc) = lrustats(iLoc) - iTmp
         end do
c      o set cachetime to the next integer
         iLoc = iiamax(cachnum,lrustats,1)
         cachetime = lrustats(iLoc) + 1
      end if

c   o check to see if the file/record is in the cache
      bNotDone = .true.
      iLoc = cachnum
      do while (bNotDone.and.(iLoc.ne.0))
         if (cachrec(iLoc).ne.iRecNum) then
            iLoc = iLoc - 1
         else
            if (cachfil(iLoc).eq.iFile) then
               bNotDone = .false.
            else
               iLoc = iLoc - 1
            end if
         end if
      end do

c   o If the correct record is in the cache, then update the cache slot and
c     set the modification flag.
      if (iLoc.ne.0) then
#ifdef _CACHE_HIST
         cachewrite = cachewrite + 1
#ifdef _CACHE_HIST_VERBOSE
         print *, '@PUTLST_IO: cache write (',
     &            iFile,',',iRecNum,') -> ',iLoc
#endif /* _CACHE_HIST_VERBOSE */
#endif /* _CACHE_HIST */
         iNdx = cachndx(iLoc)-1+iWrdNdx
         call icopy(nWords,ZList,1,icache(iNdx),1)
         cachmod(iLoc) = cachmod(iLoc) + 1
         lrustats(iLoc) = cachetime
      else

#ifdef _CACHE_BYPASS
      if (nWords.eq.iprcwd) then
#ifdef _CACHE_HIST
         cacheskip = cacheskip + 1
#ifdef _CACHE_HIST_VERBOSE
         print *, '@PUTLST_IO: cache skip  (',iFile,',',iRecNum,')'
#endif /* _CACHE_HIST_VERBOSE */
#endif /* _CACHE_HIST */
         call aces_io_write(iFile,iRecNum,ZList,iprcwd)
      else
#endif /* _CACHE_BYPASS */

c   o the record must be added to the disk record
c     Dump the least recently used record to disk if it has been modified;
c     otherwise, just copy over it.
      iVictim = iiamin(cachnum,lrustats,1)
      iNdx = cachndx(iVictim)
      if (cachmod(iVictim).ne.0) then
         iFileOut = cachfil(iVictim)
         iRecOut  = cachrec(iVictim)
#ifdef _CACHE_HIST
         cachewriteback = cachewriteback + 1
#ifdef _CACHE_HIST_VERBOSE
         print *, '@PUTLST_IO: cache write-back ',iVictim,
     &            ' -> (',iFileOut,',',iRecOut,')'
#endif /* _CACHE_HIST_VERBOSE */
#endif /* _CACHE_HIST */
         call aces_io_write(iFileOut,iRecOut,icache(iNdx),iprcwd)
      end if
#ifdef _CACHE_HIST
      cachemiss = cachemiss + 1
#ifdef _CACHE_HIST_VERBOSE
      print *, '@PUTLST_IO: cache miss  (',
     &         iFile,',',iRecNum,') -> ',iVictim
#endif /* _CACHE_HIST_VERBOSE */
#endif /* _CACHE_HIST */
      if (nWords.lt.iprcwd) then
         call aces_io_read(iFile,iRecNum,icache(iNdx),iprcwd)
         iNdx = iNdx-1+iWrdNdx
      end if
      call icopy(nWords,ZList,1,icache(iNdx),1)
      cachfil(iVictim) = iFile
      cachrec(iVictim) = iRecNum
      cachmod(iVictim) = 1
      lrustats(iVictim) = cachetime

#ifdef _CACHE_BYPASS
c     end if (nWords.eq.iprcwd)
      end if
#endif

c     end if (iLoc.ne.0)
      end if

      return
c     end subroutine putlst_io
      end

