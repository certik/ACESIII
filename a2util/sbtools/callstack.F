
c The call stack is a list of routine names. Every time you call a routine,
c the name of that routine should be pushed onto the stack. When you exit,
c you should remove the name from the stack.

c This is useful for tracking down where errors occur, but it requires
c discipline since stack management is not automatic.

c The "stack" has two components: the array of strings ("invisible" to
c the rest of the code) and callstack_curr (from include/callstack.com).
c The reason is that in many routines, the overhead from calling push/pop
c would have a major impact on performance. Other times, it is completely
c negligible. You should try to use push/pop whenever possible and revert
c to setting callstack_curr manually when performance is critical. Using
c callstack_curr too deeply defeats the purpose of a stack, but it can
c be used as a tracer.

c Prototypes:
c    subroutine callstack_init(szRoutine)
c    subroutine callstack_push(szRoutine)
c    subroutine callstack_print
c    subroutine callstack_pop
c    subroutine callstack_term

#define _CALLSTACK_LEN 32 /* max length of routine name */
#define _CALLSTACK_NUM 100 /* max depth of stack */
c#define _CALLSTACK_TRACE /* for verbose program tracing */

#ifdef _FAIL_TO_COMPILE_
c Copy this block to each routine:
#undef _CALLSTACK_COM_
#include "callstack.com" /* for callstack_curr */
      external callstack_bd
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
#endif /* _FAIL_TO_COMPILE_ */

c ----------------------------------------------------------------------

      blockdata callstack_bd
      implicit none
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
      data callstack /_CALLSTACK_NUM*' '/
      data csn /0/
      end

c ----------------------------------------------------------------------

      subroutine callstack_init(szRoutine)
      implicit none
#undef _CALLSTACK_COM_
#include "callstack.com" /* for callstack_curr */
      external callstack_bd
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
      character*(*) szRoutine
      integer i
      callstack_curr = szRoutine
      callstack(1)   = szRoutine
      do i = 2, _CALLSTACK_NUM
         callstack(i) = ' '
      end do
      csn = 1
#ifdef _CALLSTACK_TRACE
      write(*,*) '@CALLSTACK_INIT: initialized callstack for ',szRoutine
#endif
      return
      end

c ----------------------------------------------------------------------

      subroutine callstack_push(szRoutine)
      implicit none
#undef _CALLSTACK_COM_
#include "callstack.com" /* for callstack_curr */
      external callstack_bd
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
      character*(*) szRoutine
#ifdef _CALLSTACK_TRACE
      write(*,*) '@CALLSTACK_PUSH: callstack(',csn+1,') = ',szRoutine
#endif
      callstack_curr = szRoutine
      csn = csn+1
      if (csn.le._CALLSTACK_NUM) then
         callstack(csn) = szRoutine
      else
         write(*,*) '@CALLSTACK_PUSH: stack length exceeded by ',
     &              szRoutine
      end if
      return
      end

c ----------------------------------------------------------------------

      subroutine callstack_print
      implicit none
#undef _CALLSTACK_COM_
#include "callstack.com" /* for callstack_curr */
      external callstack_bd
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
      integer i, j
      write(*,*) '@CALLSTACK: ',callstack_curr
      j = csn
      if (callstack(csn).eq.callstack_curr) j = csn-1
      do i = j, 1, -1
         write(*,*) ' called by: ',callstack(i)
      end do
      return
      end

c ----------------------------------------------------------------------

      subroutine callstack_pop
      implicit none
#undef _CALLSTACK_COM_
#include "callstack.com" /* for callstack_curr */
      external callstack_bd
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
      if (csn.le.0) then
         write(*,*) '@CALLSTACK_POP: callstack is empty.'
         call c_exit(1)
      end if
#ifdef _CALLSTACK_TRACE
      if (csn.le._CALLSTACK_NUM) then
         write(*,*) '@CALLSTACK_POP: popped ',callstack(csn)
      else
         write(*,*) '@CALLSTACK_POP: popping overflow level ',csn
      end if
#endif
      csn = csn-1
      callstack_curr = 'UNKNOWN'
      if (csn.le._CALLSTACK_NUM) callstack_curr = callstack(csn)
      return
      end

c ----------------------------------------------------------------------

      subroutine callstack_term
      implicit none
#undef _CALLSTACK_COM_
#include "callstack.com" /* for callstack_curr */
      external callstack_bd
      character*(_CALLSTACK_LEN) callstack(_CALLSTACK_NUM)
      integer csn
      common /callstack_com/ callstack, csn
      save   /callstack_com/
      if (csn.gt.1) then
         write(*,*) '@CALLSTACK_TERM: callstack is not empty.'
         call callstack_print
         call c_exit(1)
      end if
#ifdef _CALLSTACK_TRACE
      write(*,*) '@CALLSTACK_TERM: done with ',callstack(1)
#endif
      callstack(1) = ' '
      csn = 0
      return
      end

