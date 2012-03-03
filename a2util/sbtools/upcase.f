      subroutine upcase(str)
      implicit none
      character*(*) str
      integer l,i,n,d
      integer iachar, iBigA, iSmlA, iSmlZ
      character*(1) achar
      iBigA = iachar('A')
      iSmlA = iachar('a')
      iSmlZ = iachar('z')
      d = iBigA-iSmlA
      l=len(str)
      do i=1,l
         n = iachar(str(i:i))
         if (iSmlA.le.n .and. n.le.iSmlZ) str(i:i) = achar(n+d)
      end do
      return
      end
