
c This routine checks to see if any invalid keywords are left in line.

      subroutine nl_term
      implicit none
#include "namelist.com"
#include "callstack.com"
      integer i
      callstack_prev=callstack_curr
      callstack_curr='NL_TERM'
      if (prt_nl) call nl_prtbot
      if (nlnumline.ne.0) then
         write(*,*) '@NL_TERM: invalid namelist entries'
         write(*,*) '  Errors somewhere around '
         write(*,*)
         do i = 1, nlnumline
            write(*,*) nltext(i)
         end do
         call errex
      end if
      callstack_curr=callstack_prev
      return
      end

