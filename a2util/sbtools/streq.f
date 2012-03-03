      function streq(str1,str2,insensitive)
      implicit none
      logical streq
      character*(*) str1,str2
      logical insensitive
      integer l1,l2,i,strlen,ic1,ic2,d
      integer iachar, iBigA, iSmlA, iSmlZ
      character*1 achar
      iBigA=iachar('A')
      iSmlA=iachar('a')
      iSmlZ=iachar('z')
      d = iBigA-iSmlA
      streq=.false.
      l1=strlen(str1)
      l2=strlen(str2)
      if (l1.ne.l2) return
      do i=1,l1
         ic1=iachar(str1(i:i))
         ic2=iachar(str2(i:i))
         if (insensitive) then
            if (iSmlA.le.ic1 .and. ic1.le.iSmlZ) ic1=ic1+d
            if (iSmlA.le.ic2 .and. ic2.le.iSmlZ) ic2=ic2+d
         end if
         if (ic1.ne.ic2) return
      end do
      streq=.true.
      return
      end
