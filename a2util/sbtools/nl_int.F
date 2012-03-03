
c This routine reads an integer value.  If it is not in the namelist,
c it defaults to the value passed in.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_int(key,def,var)
      implicit none

      character*(*) key
      integer def,var

#include "namelist.com"
#include "callstack.com"

      character*(nllinelen) val
      logical nl_key,present
      integer err

      callstack_prev=callstack_curr
      callstack_curr='NL_INT'

      var=def

      val=' '
      present=nl_key(key,val)
      if (present) then
         call str2int(val,var,err)
         if (err.ne.0) then
            write(*,*) '@NL_INT: invalid integer for ',key
            call errex
         end if
      end if

      if (prt_nl) then
         if (var.eq.def) then
            write(*,900) key,'integer',var
  900       format(a20,2x,a10,2x,i20)
         else
            write(*,910) key,'integer',var,def
  910       format(a20,2x,a10,2x,i20,2x,i20)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_int
      end

