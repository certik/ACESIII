      function ncindex(str,c,ind)
      implicit none
      integer ncindex
      character*(*) str
      character*1 c
      integer ind,strlen,i,slen
      ncindex=0
      if (ind.lt.0) return
      slen=strlen(str)
      i=ind+1
      do while (.true.)
         if (i.gt.slen) return
         if (str(i:i).eq.c) then
            ncindex=i
            return
         end if
         i=i+1
      end do
      end
