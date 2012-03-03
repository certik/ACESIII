
      subroutine list_next_ele(str,ptr,ele,err)
      implicit none
      character*(*) str,ele
      integer ptr,err
      integer l,pos
      integer  strlen,ncindex
      external strlen,ncindex
      err=0
      l=strlen(str)
      ptr=max(ptr,0)
      ptr=min(ptr,l+1)
      if (ptr.eq.l+1) then
         err=-1
         ele=' '
         return
      end if
      pos=ncindex(str,',',ptr)
      if (pos.eq.0) then
         ele=str(ptr+1:l)
         ptr=l+1
      else
         ele=str(ptr+1:pos-1)
         ptr=pos
      end if
      return
      end

