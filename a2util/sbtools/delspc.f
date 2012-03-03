      subroutine delspc(instr,outstr)
      implicit none
      character*(*) instr,outstr
      integer j,islen,oslen,strlen,i
      character*1 achar,spc,tab
      spc=achar(32)
      tab=achar(9)
      j=1
      islen=strlen(instr)
      oslen=strlen(outstr)
      do i=1,islen
         if (instr(i:i).ne.spc .and. instr(i:i).ne.tab) then
            outstr(j:j)=instr(i:i)
            j=j+1
         end if
      end do
      if (j.gt.oslen+1) then
         print *, '@DELSPC: wrote out of bounds'
         call c_exit(1)
      end if
cYAU: this loop used to run to oslen+1, which makes no sense
      do i=j,oslen
         outstr(i:i)=spc
      end do
      return
      end
