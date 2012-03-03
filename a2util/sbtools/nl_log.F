
c This routine reads a logical value.  If the keyword has a value,
c the value must be one of the (case-insensitive) following:
c   yes, true,   1,  on,  y,  t
c   no,  false,  0,  off  n,  f
c If it does not have a value, it is set to .true.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_log(key,def,var)
      implicit none

      character*(*) key
      logical var,def

#include "namelist.com"
#include "callstack.com"

      character*5 val,defval
      logical nl_key,present,streq
      integer strlen

      callstack_prev=callstack_curr
      callstack_curr='NL_LOG'

      var=def

      val=' '
      present=nl_key(key,val)
      if (present) then
         var=.true.
         if (strlen(val).gt.0) then
            if (streq(val,'yes',.true.) .or.
     &          streq(val,'true',.true.) .or.
     &          streq(val,'y',.true.) .or.
     &          streq(val,'t',.true.) .or.
     &          streq(val,'1',.true.) .or.
     &          streq(val,'on',.true.)) then
               var=.true.
            else if (streq(val,'no',.true.) .or.
     &               streq(val,'false',.true.) .or.
     &               streq(val,'n',.true.) .or.
     &               streq(val,'f',.true.) .or.
     &               streq(val,'0',.true.) .or.
     &               streq(val,'off',.true.)) then
               var=.false.
            else
               write(*,*) '@NL_LOG: invalid value for ',key
               call errex
            end if
         end if
      end if

      if (prt_nl .and. .not.streq(key,'print_nl',.true.)) then
         val='false'
         if (var) val=' true'
         defval='false'
         if (def) defval=' true'
         if (val.eq.defval) then
            write(*,900) key,'logical',val
  900       format(a20,2x,a10,2x,a20)
         else
            write(*,910) key,'logical',val,defval
  910       format(a20,2x,a10,2x,a20,2x,a20)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_log
      end

