
c This routine initializes the ACES runtime environment. It or aces_init
c should be called before any other ACES routine.

      subroutine aces_init_rte
      implicit none

c INTERNAL VARIABLES
      character*80 szJOBARC
      integer       iJOBARC
      logical bTmp

c COMMON BLOCKS
#include "machsp.com"
#include "icdacc.com"

c ----------------------------------------------------------------------

c   o initialize the machsp common block
      call aces_com_machsp
#ifdef _ASSERT
      if ((iintln.ne.4).and.(iintln.ne.8)) then
         print *, '@ACES_INIT_RTE: Assertion failed.'
         print *, '                iintln = ',iintln
         print *, '                ifltln = ',ifltln
         call aces_exit(1)
      end if
#endif /* _ASSERT */

c   o initialize the icdacc common block

cYAU - Okay, what is this "idacc?" all about...
c
c    First off, the name refers to (D)irect (ACC)ess. The ultimate problem
c is that the Fortran standard lets RECL take any unit when specifying
c the record length of direct access, unformatted files. It's not
c necessarily compiler, OS, or architecture dependent. In fact, some
c configurations let you specify the units of RECL with an environment
c variable!
c    To aid portability, idaccm was introduced by Nevin as a scaling
c factor for RECL. In essence, idaccm is the number of RECL units
c per integer; therefore, it could take any of the values: 1, 2, 4, or 8.
c    If there is an accompanying Fortran 9x compiler, then a short
c program will tell you whether to define _RECL_IS_WORDS_:
c       program test
c          integer :: i=0, len
c          inquire(iolength=len) i
c          print *, len
c       end program test
c Whatever number is printed, that should be the value of idaccm. For
c the time being, if RECL is not in bytes, then we will assume 4-byte
c words. It is possible RECL is 8-byte words, but that is very rare.

#ifdef _RECL_IS_WORDS_
#ifdef F_64BIT
      idaccm = 2
#else
      idaccm = 1
#endif /* F_64BIT */
#else
      idaccm = iintln
#endif /* _RECL_IS_WORDS_ */

#ifdef _DEBUG
c   o try to print this info just once
c     (should we bother trapping this with err/iostat?)
      call gfname('JOBARC',szJOBARC,iJOBARC)
      inquire(file=szJOBARC(1:iJOBARC),exist=bTmp)
      if (.not.bTmp) then
         print *, '@ACES_INIT_RTE: Build-specific parameters:'
         print *, '                bytes per integer = ',iintln
         print *, '                bytes per double  = ',ifltln
         print *, '                recls per integer = ',idaccm
      end if
#endif

      return
c     end subroutine aces_init_rte
      end

