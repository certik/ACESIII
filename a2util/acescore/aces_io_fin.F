
c This routine finalizes the I/O subsystem.

      subroutine aces_io_fin
      implicit none

c EXTERNAL FUNCTIONS
      integer iiamax

c INTERNAL VARIABLES
      integer iUnit, iLenMOIO, iStat
      logical bOpened

c COMMON BLOCKS
#include "lists.com"
#include "sympop.com" /* for isytyp */
#include "cache.com" /* for cache histogram data */
#include "auxcache.com" /* for quikget */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iUnit = 0
c   o assert quikget is empty
      iLenMOIO = iiamax(_MAX_IO_GRPS*_MAX_IO_FAMS,quikget,1)
      if ((iLenMOIO.ne.1).or.(quikget(1,1).ne.0)) then
         print *, '@ACES_IO_FIN: Assertion failed.'
         print *, '   quikget col = ',1+((iLenMOIO-1)/_MAX_IO_GRPS)
         print *, '   quikget row = ',1+mod((iLenMOIO-1),_MAX_IO_GRPS)
         iUnit = 1
      end if
      if (iUnit.ne.0) call aces_exit(iUnit)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c   o make sure the I/O subsystem is up
      if (.not.bIOUp) return

c   o flush the automatic cache
      call aces_cache_flush

#ifdef _CACHE_HIST
c   o print the cache-effectiveness data
      print *, '@ACES_IO_FIN: caching effectiveness via total counts'
      print *, '              total misses      = ',cachemiss
      print *, '              total skips       = ',cacheskip
      print *, '              total reads  from = ',cacheread
      print *, '              total writes to   = ',cachewrite
      print *, '              total writes back = ',cachewriteback
#endif

c   o update storage records only if a list has been modified
      if (bIOMod) then

c      o get the total number of arrays
         iLenMOIO = _MAX_IO_GRPS * _MAX_IO_FAMS

c      o finalize the lists common block
         call putrec(1,'JOBARC','MOIOVEC',iLenMOIO,moio)
         call putrec(1,'JOBARC','MOIOWRD',iLenMOIO,moiowd)
         call putrec(1,'JOBARC','MOIOSIZ',iLenMOIO,moiosz)
         call putrec(1,'JOBARC','MOIODIS',iLenMOIO,moiods)
         call putrec(1,'JOBARC','MOIOFIL',iLenMOIO,moiofl)
         call putrec(1,'JOBARC','MOIOMXSZ',iLenMOIO,moiomxsz)

c      o write out the distribution types
         call putrec(1,'JOBARC','ISYMTYP',2*_MAX_IO_FAMS,isytyp)

c      o finalize the io_ptrs common block
         call putrec(1,'JOBARC','TOTRECMO',_MAX_IO_LUNS,pRec)
         call putrec(1,'JOBARC','TOTWRDMO',_MAX_IO_LUNS,iIntOff)

c      o reset the I/O modification flag
         bIOMod = .false.

c     end if (bIOMod)
      end if

#ifdef _DEBUG_ACES_IO
c   o reset the I/O stats so nothing gets "lost"
      call aces_cache_reset
      call aces_io_reset
#endif

#ifndef OPEN_LISTS_AS_NEEDED
c   o close the external file units
      do iUnit = _FIRST_IO_LUN, _FIRST_IO_LUN-1+_MAX_IO_LUNS
         inquire(unit=iUnit,opened=bOpened,err=666,iostat=iStat)
         if (bOpened) then
            close(unit=iUnit,status='KEEP',err=666,iostat=iStat)
         end if
      end do
#endif

c   o turn off the I/O subsystem flag
      bIOUp = .false.

      return

c   o I/O error
 666  print *, '@ACES_IO_FIN: I/O error'
      print '(/)'
      call aces_io_error('ACES_IO_FIN',iUnit,iStat)

c     end subroutine aces_io_fin
      end

