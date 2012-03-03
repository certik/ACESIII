
      subroutine list_prev_num(str,ptr,num,err)
      implicit none
      character*(*) str
      integer ptr,err
      double precision num
      integer l,pos
      integer  strlen,pcindex
      external strlen,pcindex
      err=0
      l=strlen(str)
      ptr=max(ptr,0)
      ptr=min(ptr,l+1)
      if (ptr.eq.0) then
         err=-1
         num=0.d0
         return
      end if
      pos=pcindex(str,',',ptr)
      if (pos.eq.0) then
         call str2dbl(str(1:ptr-1),num,err)
      else
         call str2dbl(str(pos+1:ptr-1),num,err)
      end if
      if (err.eq.0) ptr=pos
      return
      end

