
      subroutine list_next_num(str,ptr,num,err)
      implicit none
      character*(*) str
      integer ptr,err
      double precision num
      integer l,pos
      integer  strlen,ncindex
      external strlen,ncindex
      err=0
      l=strlen(str)
      ptr=max(ptr,0)
      ptr=min(ptr,l+1)
      if (ptr.eq.l+1) then
         err=-1
         num=0.d0
         return
      end if
      pos=ncindex(str,',',ptr)
      if (pos.eq.0) then
         call str2dbl(str(ptr+1:l),num,err)
         pos=l+1
      else
         call str2dbl(str(ptr+1:pos-1),num,err)
      end if
      if (err.eq.0) ptr=pos
      return
      end

