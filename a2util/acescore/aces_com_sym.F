
c This routine initializes the sym common block.

      subroutine aces_com_sym
      implicit none

c INTERNAL VARIABLES
      integer irrep

c COMMON BLOCKS
#include "sym.com"
#include "syminf.com" /* for nirrep */

c ----------------------------------------------------------------------

c   o get the number of orbitals per irrep
      call getrec(-1,'JOBARC','SYMPOPOA',nirrep,pop(1,1))
      call getrec(-1,'JOBARC','SYMPOPOB',nirrep,pop(1,2))
      call getrec(-1,'JOBARC','SYMPOPVA',nirrep,vrt(1,1))
      call getrec(-1,'JOBARC','SYMPOPVB',nirrep,vrt(1,2))

c   o calculate the number of orbital pairs
      nt(1)   = 0
      nt(2)   = 0
      nfmi(1) = 0
      nfmi(2) = 0
      nfea(1) = 0
      nfea(2) = 0
      do irrep = 1, nirrep
c      o T1(V,O)
         nt(1)   = nt(1)   + vrt(irrep,1)*pop(irrep,1)
         nt(2)   = nt(2)   + vrt(irrep,2)*pop(irrep,2)
c      o F(M,I)
         nfmi(1) = nfmi(1) + pop(irrep,1)*pop(irrep,1)
         nfmi(2) = nfmi(2) + pop(irrep,2)*pop(irrep,2)
c      o F(E,A)
         nfea(1) = nfea(1) + vrt(irrep,1)*vrt(irrep,1)
         nfea(2) = nfea(2) + vrt(irrep,2)*vrt(irrep,2)
      end do

      return
c     end subroutine aces_com_sym
      end

