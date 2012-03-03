
c This routine should be the common exit handler for integer overflows.
c Since this exception is not usually trapped, the tests must be hand-coded.

c INPUT
c char*(*) SZCALLER : name of calling routine
c integer  ITAG     : a caller-specific integer tag

c USAGE
c   o send the current line number to the exit routine (requires cpp processing)
c     if (i.lt.0) call trap_intovf('MY_ROUTINE',__LINE__)
c
c   o send the overflowed variable and its value to the exit routine
c     if (i.lt.0) call trap_intovf('MY_ROUTINE:I',i)

      subroutine trap_intovf(szCaller,iTag)
      implicit none

      character*(*) szCaller
      integer iTag

      print *, '@',szCaller,',',iTag,': integer overflow'
      call c_exit(1)

      end

