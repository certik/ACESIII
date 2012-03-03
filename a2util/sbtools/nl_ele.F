
c This routine reads a string which MUST be one of the strings given in list.
c list contains n values numbered 1..n.
c
c Variables which have default values may safely be passed in both
c as the default value and as the value read in.

      subroutine nl_ele(key,list,n,def,var)
      implicit none

      character*(*) key,list(*)
      integer n,def,var

#include "namelist.com"
#include "callstack.com"

      character*(nllinelen) val
      logical nl_key,present,streq
      integer i,strlen

      callstack_prev=callstack_curr
      callstack_curr='NL_ELE'

      var=def

      val=' '
      present=nl_key(key,val)
      if (present) then
         var=0
         do i=1,n
            if (streq(val,list(i),.true.)) var=i
         end do
      end if
      if (var.lt.1 .or. var.gt.n) then
         write(*,*) '@NL_ELE: invalid element'
         write(*,*)
         write(*,*) ' valid input key: ', key
         write(*,*)
         write(*,*) ' but invalid value: ', val
         write(*,*)
         write(*,*) ' possible values: '
         write(*,*)
         do i = 1, n
            write(*,*) list(i)
         end do
         call errex
      end if

      if (prt_nl) then
         if (var.eq.def) then
            write(*,900) key,'element',list(var)(1:strlen(list(var)))
 900        format(a20,2x,a10,2x,a20)
         else
            write(*,910) key,'element',
     &                   list(var)(1:strlen(list(var))),
     &                   list(def)(1:strlen(list(def)))
 910        format(a20,2x,a10,2x,a20,2x,a20)
         end if
      end if

      callstack_curr=callstack_prev
      return
c     end subroutine nl_ele
      end

