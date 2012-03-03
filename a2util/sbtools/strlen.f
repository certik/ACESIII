      function strlen(str)
      implicit none
      integer strlen
      character*(*) str
      character*1 spc,tab,zero
      integer i
      i=len(str)
      spc=achar(32)
      tab=achar(9)
      zero=achar(0)
      do while (i.gt.0 .and.
     &          (str(i:i).eq.spc .or.
     &           str(i:i).eq.tab .or.
     &           str(i:i).eq.zero)
     &         )
         i=i-1
      end do
      strlen=i
      return
      end
