
      subroutine list_next_int(str,ptr,i,err)
      implicit none
      character*(*) str
      integer ptr,err,i
      integer l,pos
      integer  strlen,ncindex
      external strlen,ncindex
      err=0
      l=strlen(str)
      ptr=max(ptr,0)
      ptr=min(ptr,l+1)
      if (ptr.eq.l+1) then
         err=-1
         i=0
         return
      end if
      pos=ncindex(str,',',ptr)
      if (pos.eq.0) then
         call str2int(str(ptr+1:l),i,err)
         pos=l+1
      else
         call str2int(str(ptr+1:pos-1),i,err)
      end if
      if (err.eq.0) ptr=pos
      return
      end

