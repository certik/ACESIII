
c This routine prints a summary of all the records in the JOBARC file up to
c the next open slot.

      subroutine aces_ja_summary
      implicit none

c INTERNAL VARIABLES
      integer i
      logical bNotDone

c COMMON BLOCKS
#include "jobarc.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert job archive subsystem is up
      if (.not.bJAUp) then
         print *, '@ACES_JA_SUMMARY: Assertion failed.'
         print *, '   bJAUp = ',bJAUp
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif /* _ASSERT */

c ----------------------------------------------------------------------

      print '(/)'
      print '(8x,64a)', ('-',i=1,64)
      print '(8x,21x,a22)', 'SUMMARY OF JOB ARCHIVE'
      print '(8x,64a)', ('-',i=1,64)
      print '(8x,a5,8x,3x,a4,3x,8x,a13,8x,a12)',
     &          'INDEX', 'NAME', 'ADDRESS (INT)', 'LENGTH (INT)'
      print '(8x,64a)', ('-',i=1,64)
      i = 1
      bNotDone = .true.
      do while (bNotDone)
      print '(8x,i5,8x,a1,a8,a1,8x,i13,8x,i12)',
     &           i, '"',marker(i),'"', rloc(i), rsize(i)
      if (iand(i,3).eq.0) print '(/)'
      bNotDone = ((marker(i).ne.'OPENSLOT').and.(i.le._MAX_JA_RECS))
      i = i+1
      end do
      print '(8x,64a)', ('-',i=1,64)
      print '(/)'

      return
c     end subroutine aces_ja_summary
      end

