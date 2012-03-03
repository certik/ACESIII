
c This routine prints a summary of all the lists in the storage files.

      subroutine aces_io_summary
      implicit none

c INTERNAL VARIABLES
      integer ihFile, iUnit
      integer i, max_grp, iGrp, iFam
      logical bTmp, bDoFam
#include "szfiles.dat"

c COMMON BLOCKS
#include "lists.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert I/O subsystem is up
      if (.not.bIOUp) then
         print *, '@ACES_IO_SUMMARY: Assertion failed.'
         print *, '   bIOUp = ',bIOUp
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      print '(/)'
      print '(5x,70a)', ('-',i=1,70)
      print '(5x,23x,a24)', 'SUMMARY OF STORAGE FILES'
      print '(5x,70a)', ('-',i=1,70)
      print '(5x,a5,8x,a4,14x,a4,11x,a8,7x,a9)',
     &          'GROUP', 'ROWS', 'COLS', 'PHYS REC', 'INT INDEX'

      do ihFile = 1, _MAX_IO_LUNS
         print '(5x,70a)', ('-',i=1,70)
         print '(5x,2a)', 'file name: ', szFiles(ihFile)
         iUnit = _FIRST_IO_LUN - 1 + ihFile
         do iFam = 1+(ihFile-1)*_FAMS_PER_FILE, ihFile*_FAMS_PER_FILE

c         o see if this family has any lists
            bDoFam = .false.
            max_grp = 0
            do iGrp = 1, _MAX_IO_GRPS
               bTmp = (moio  (iGrp,iFam).gt.0) .or.
     &                (moiowd(iGrp,iFam).gt.0) .or.
     &                (moiosz(iGrp,iFam).gt.0) .or.
     &                (moiods(iGrp,iFam).gt.0) .or.
     &                (moiofl(iGrp,iFam).gt.0)
               if (bTmp) then
                  bDoFam = .true.
                  max_grp = iGrp
               end if
            end do

c         o print this family
            if (bDoFam) then
               print '(/)'
               print '(5x,a1,i3,a1)', '(',iFam,')'
               do iGrp = 1, max_grp
                  if (moiofl(iGrp,iFam).gt.0) then
                     if (moiofl(iGrp,iFam).ne.iUnit) then
                        print '(5x,a)', 'WARNING: Array is misplaced!'
                     end if
                  end if
                  print '(5x,i2,7x,i8,10x,i8,10x,i9,10x,i6)',
     &                       iGrp, moiosz(iGrp,iFam), moiods(iGrp,iFam),
     &                       moio(iGrp,iFam), moiowd(iGrp,iFam)
               end do
            end if

c        end do iFam = 1+(ihFile-1)*_FAMS_PER_FILE, ihFile*_FAMS_PER_FILE
         end do
c     end do ihFile = 1, _MAX_IO_LUNS
      end do

      print '(5x,70a)', ('-',i=1,70)
      print '(/)'

      return
c     end subroutine aces_io_summary
      end

