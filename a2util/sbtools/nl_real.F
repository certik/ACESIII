
c This routine reads an real value.  If it isn't in the namelist,
c it defaults to the value passed in.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_real(key,def,var)
      implicit none

      character*(*) key
      double precision def,var

#include "namelist.com"
#include "callstack.com"

      character*(nllinelen) val
      logical nl_key,present
      integer err

      callstack_prev=callstack_curr
      callstack_curr='NL_REAL'

      var=def

      val=' '
      present=nl_key(key,val)
      if (present) then
         call str2dbl(val,var,err)
         if (err.ne.0) then
            write(*,*) '@NL_REAL: invalid real for ',key
            call errex
         end if
      end if

      if (prt_nl) then
         if (var.eq.def) then
            write(*,900) key,'real',var
  900       format(a20,2x,a10,2x,f20.8)
         else
            write(*,910) key,'real',var,def
  910       format(a20,2x,a10,2x,f20.8,2x,f20.8)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_real
      end

