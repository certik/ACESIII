
      subroutine list_prev_ele(str,ptr,ele,err)
      implicit none
      character*(*) str,ele
      integer ptr,err
      integer l,pos
      integer  strlen,pcindex
      external strlen,pcindex
      err=0
      l=strlen(str)
      ptr=max(ptr,0)
      ptr=min(ptr,l+1)
      if (ptr.eq.0) then
         err=-1
         ele=' '
         return
      end if
      pos=pcindex(str,',',ptr)
      if (pos.eq.0) then
         ele=str(1:ptr-1)
         ptr=0
      else
         ele=str(pos+1:ptr-1)
         ptr=pos
      end if
      return
      end

