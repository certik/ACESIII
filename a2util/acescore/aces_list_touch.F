
c This routine creates storage statistics for a new double-precision array,
c and extends the proper storage unit by the requested amount of data.
c It will bomb if there is not enough disk space.

c INPUT
c int NCOLS : the number of columns in the array
c int NROWS : the number of rows    in the array
c int ILEFT  : the left  index in MOIO
c int IRIGHT : the right index in MOIO
c int INOFRAG : a behavior flag
c               >= 0; start the array at the next available location
c                < 0; start the array at the next available physical record
c               NOTE: A few routines (e.g., aces_family_redist) require all the
c                     arrays in one family to be packed together with no space
c                     between them.

c#define _FLUSH_WITH_NaN
c#define _TRAP_LIST_FAM 6 /* iRight */
c#define _TRAP_LIST_GRP 1 /* iLeft */
c#define _WARN_LIST_FAM 165
c#define _WARN_LIST_GRP 1

      subroutine aces_list_touch(nCols,nRows,iLeft,iRight,iNoFrag)
      implicit none

c ARGUMENTS
      integer nCols, nRows, iLeft, iRight, iNoFrag

c INTERNAL VARIABLES
#include "szfiles.dat"
      character*80 szFile
      integer      iLength
      integer ihFile, iUnit, pRecOld, nWords
      integer iStat, iJunk
      logical bExist, bOpened

c COMMON BLOCKS
#include "lists.com"
#include "machsp.com" /* for iintfp */
#include "icdacc.com" /* for idaccm */
#include "filspc.com" /* for iprcwd */
#include "auxcache.com" /* for quikget */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iStat = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_LIST_TOUCH: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         iStat = 1
      end if
c   o assert columns and rows are whole
      if ((nCols.lt.0).or.(nRows.lt.0)) then
         print *, '@ACES_LIST_TOUCH: Assertion failed.'
         print *, '   nRows = ',nRows
         print *, '   nCols = ',nCols
         iStat = 1
      end if
c   o assert left and right indices are properly bound
      if ((iLeft .lt.1).or.(_MAX_IO_GRPS.lt.iLeft ).or.
     &    (iRight.lt.1).or.(_MAX_IO_FAMS.lt.iRight)    ) then
         print *, '@ACES_LIST_TOUCH: Assertion failed.'
         print *, '   iLeft  = ',iLeft
         print *, '   iRight = ',iRight
         print *, '   MOIO rows = ',_MAX_IO_GRPS
         print *, '   MOIO cols = ',_MAX_IO_FAMS
         iStat = 1
      end if
c   o assert quikget is not set if the list does not exist yet
      if (moio(iLeft,iRight).eq.0) then
      if (quikget(iLeft,iRight).ne.0) then
         print *, '@ACES_LIST_TOUCH: Assertion failed.'
         print *, '   iLeft  = ',iLeft
         print *, '   iRight = ',iRight
         print *, '   quikget(,) = ',quikget(iLeft,iRight)
         iStat = 1
      end if
      end if
      if (iStat.ne.0) call aces_exit(iStat)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

#if defined(_TRAP_LIST_FAM) && defined(_TRAP_LIST_GRP)
      if ((iRight.eq._TRAP_LIST_FAM).and.
     &    (iLeft .eq._TRAP_LIST_GRP)     ) call c_abort
#endif
#if defined(_WARN_LIST_FAM) && defined(_WARN_LIST_GRP)
      if ((iRight.eq._WARN_LIST_FAM).and.
     &    (iLeft .eq._WARN_LIST_GRP)     ) then
         print *, '@ACES_LIST_TOUCH: touching ',iLeft,iRight
      end if
#endif

c   o if the array already exists, check to make sure it is still the same size
      if (moio(iLeft,iRight).ne.0) then
         if (moiosz(iLeft,iRight).ne.nRows.or.
     &       moiods(iLeft,iRight).ne.nCols    ) then
            print *, '@ACES_LIST_TOUCH: list (',iLeft,',',iRight,
     &               ') already exists and is of different dimensions'
            print *, '   currently: ',moiosz(iLeft,iRight),' x ',
     &                                moiods(iLeft,iRight)
            print *, '   requested: ',nRows,' x ',nCols
            call aces_exit(1)
         end if
         return
      end if

c   o ihFile is the integer handle of the file storing the array
      ihFile = 1 + (iRight-1)/_FAMS_PER_FILE

c   o open and/or create the unit
      call aces_io_open(ihFile)

c   o save the current record pointer
      pRecOld = pRec(ihFile)

c   o shift the free-space pointers if desired
      if ((iNoFrag.lt.0).and.(iIntOff(ihFile).ne.0)) then
         pRec(ihFile) = pRec(ihFile) + 1
         iIntOff(ihFile) = 0
      end if

c   o convert the handle to an external unit
      iUnit = ihFile-1+_FIRST_IO_LUN

c   o assign storage statistics
      moio  (iLeft,iRight) = pRec(ihFile)
      moiowd(iLeft,iRight) = 1 + iIntOff(ihFile)
      moiosz(iLeft,iRight) = nRows
      moiods(iLeft,iRight) = nCols
      moiofl(iLeft,iRight) = iUnit

c   o record the number of doubles as an upper bound for resizing
      moiomxsz(iLeft,iRight) = nRows*nCols

cYAU - There should be an isytyp modifier here.
c      isytyp(1,iRight) = ???
c      isytyp(2,iRight) = ???

c   o re-position the free-space pointers
      nWords = nRows*nCols*iintfp + iIntOff(ihFile)
      pRec(ihFile) = pRec(ihFile) + (nWords/iprcwd)
      iIntOff(ihFile) = mod(nWords,iprcwd)

#ifdef _2GB_FILES
      if (iprcln*pRec(ihFile).lt.0) then
         print *, '@ACES_LIST_TOUCH: Files over 2GB are not supported',
     &            ' by this binary.'
         call aces_exit(1)
      end if
#endif

c   o arrest disk space by writing to the highest record
      if (pRec(ihFile).gt.pRecOld) then
         write(unit=iUnit,rec=pRec(ihFile),err=666,iostat=iStat) 0
      end if

#ifdef _DEBUG_ACES_IO
c   o do nasty things to trap calling getlst before putlst
#ifdef _FLUSH_WITH_NaN
      call aces_list_memset(iLeft,iRight,255)
#else
      call aces_list_memset(iLeft,iRight,0)
#endif
#endif

c   o mark MOIO as modified
      bIOMod = .true.

      return

c   o write-error jump
 666  call gfname(szFiles(ihFile),szFile,iLength)
      print *, '@ACES_LIST_TOUCH: An error occurred while attempting ',
     &         'to extend the size of the file:'
      print *, '                  ',szFile(1:iLength)
      call aces_io_error('ACES_LIST_TOUCH',iUnit,iStat)

c     end subroutine aces_list_touch
      end

