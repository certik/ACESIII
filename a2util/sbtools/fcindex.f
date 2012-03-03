      function fcindex(str,c)
      implicit none
      integer fcindex
      character*(*) str
      character*1 c
      integer i,strlen,slen
      i=1
      slen=strlen(str)
      do while (str(i:i).ne.c .and. i.le.slen)
         i=i+1
      end do
      fcindex=0
      if (i.le.slen) fcindex=i
      return
      end
