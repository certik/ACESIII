
c This routine initializes the syminf common block.
c The code for the direct product table was taken straight from VMOL.

      subroutine aces_com_syminf
      implicit none

c INTERNAL VARIABLES
      integer i, j, k
      integer zMin, zMax
      logical bTmp

c COMMON BLOCKS
#include "syminf.com"

c INLINE FUNCTIONS
      integer ka
      ka(i,k) = (i-1)/(2**(k-1)) - 2*((i-1)/(2**k))

c ----------------------------------------------------------------------

c   o nirrep
      call getrec(-1,'JOBARC','COMPNIRR',1,nirrep)

c   o irrep a/b
c???      norbs=nocco(1)+nvrto(1)
c???      call getrec(-1,'JOBARC','IRREPALP',norbs,irrepa)
c???      call getrec(-1,'JOBARC','IRREPBET',norbs,irrepb)

c   o dirprd
      do i = 1, 8
      do j = 1, 8
         zMin = 1
         zMax = 1
         do k = 1, 3
            zMin = zMin + ( min(ka(i,k),ka(j,k)) * 2**(k-1) )
            zMax = zMax + ( max(ka(i,k),ka(j,k)) * 2**(k-1) )
         end do
         dirprd(i,j) = 1 + zMax - zMin
      end do
      end do

#ifdef _DEBUG
c   o try to print this info just once
c     (should we bother trapping this with err/iostat?)
      inquire(file='JAINDX',exist=bTmp)
      if (.not.bTmp) then
         if (nirrep.gt.0) then
            print *, '@ACES_COM_SYMINF: direct product table ='
            do i = 1, nirrep
               print '(8i3)', (dirprd(i,j),j=1,nirrep)
            end do
         end if
      end if
#endif

      return
c     end subroutine aces_com_syminf
      end

