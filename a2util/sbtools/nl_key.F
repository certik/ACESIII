
c This searches through the namelist for a key.  It removes the key and
c value from the namelist and returns the value as a string.  It returns
c .true. if the key was found.

      function nl_key(fkey,val)
      implicit none

      logical nl_key
      character*(*) fkey,val

#include "namelist.com"
#include "callstack.com"

      character*(nllinelen) fullkey,key,testkey,testval,before,after
      integer ast,ind,loc1,loc2,locsp,loceq,l,j,strlen,i
      integer  nindex,ncindex,findex,fcindex
      external nindex,ncindex,findex,fcindex
      character*1 c,sp,tab

      callstack_prev = callstack_curr
      callstack_curr = 'NL_KEY'

      fullkey=fkey
      sp=char(32)
      tab=char(9)

c First, check to see if fullkey has an '*' in it.  When done:
c  fullkey  : the name of the fullkey
c  key      : the smallest abbreviation allowed
      call upcase(fullkey)
      ast=fcindex(fullkey,'*')
      if (ast.eq.1) then
         write(*,*) '@NL_KEY-F:  a key cannot begin with *'
         call errex
      end if
      key=fullkey
      if (ast.gt.1) then
         fullkey=fullkey(1:ast-1)//fullkey(ast+1:strlen(fullkey))
         key=fullkey(1:ast-1)
      end if

c Go though each line and look for the key
      do i=1,nlnumline
         ind=0
 10      continue
         loc1=nindex(nltext(i),key,ind)
         if (loc1.gt.0) then
c           Key found... check to make sure that it starts at the beginning
c           of a line or is preceeded by spaces
            if (loc1.gt.1) then
               c=nltext(i)(loc1-1:loc1-1)
               if (c.ne.sp .and. c.ne.tab) then
                  ind=loc1
                  goto 10
               end if
            end if
c           Get the key, value, and start and end locations
            l=strlen(nltext(i))
            locsp=ncindex(nltext(i),sp//tab,loc1)
            loceq=ncindex(nltext(i),'=',loc1)
            if (locsp.gt.0.and.loceq.gt.0.and.loceq.gt.locsp) loceq=0
            if (locsp.eq.0 .and. loceq.eq.0) then
               testkey=nltext(i)(loc1:l)
               testval=' '
               loc2=l+1
            else if (locsp.eq.0) then
               testkey=nltext(i)(loc1:loceq-1)
               testval=nltext(i)(loceq+1:l)
               loc2=l+1
            else if (loceq.eq.0) then
               testkey=nltext(i)(loc1:locsp-1)
               testval=' '
               loc2=locsp+1
            else
               testkey=nltext(i)(loc1:loceq-1)
               testval=nltext(i)(loceq+1:locsp-1)
               loc2=locsp+1
            end if
c           Check to see that this key is correct
            if (findex(fullkey,testkey).ne.1) then
               ind=loc1
               goto 10
            end if
c           Remove this key from the namelist text
            if (loc1.eq.1) then
               before=' '
            else
               before=nltext(i)(1:loc1-1)
            end if
            if (loc2.gt.l) then
               after=' '
            else
               after=nltext(i)(loc2:l)
            endif
            nltext(i)=before(1:strlen(before))//' '//
     &                after(1:strlen(after))
c           Remove this line if it is now empty
            if (strlen(nltext(i)).eq.0) then
               do j=i,nlnumline-1
                  nltext(j)=nltext(j+1)
               end do
               nltext(nlnumline)=' '
               nlnumline=nlnumline-1
            end if
            val=testval
            nl_key=.true.
            callstack_curr = callstack_prev
            return
         end if
      end do

      nl_key=.false.

      callstack_curr = callstack_prev
      return
c     end function nl_key
      end

