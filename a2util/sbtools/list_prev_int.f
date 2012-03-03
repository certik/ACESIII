
      subroutine list_prev_int(str,ptr,i,err)
      implicit none
      character*(*) str
      integer ptr,err,i
      integer l,pos
      integer  strlen,pcindex
      external strlen,pcindex
      err=0
      l=strlen(str)
      ptr=max(ptr,0)
      ptr=min(ptr,l+1)
      if (ptr.eq.0) then
         err=-1
         i=0
         return
      end if
      pos=pcindex(str,',',ptr)
      if (pos.eq.0) then
         call str2int(str(1:ptr-1),i,err)
      else
         call str2int(str(pos+1:ptr-1),i,err)
      end if
      if (err.eq.0) ptr=pos
      return
      end

