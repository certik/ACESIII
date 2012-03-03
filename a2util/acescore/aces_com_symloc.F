
c This routine initializes the symloc common block.

      subroutine aces_com_symloc
      implicit none

c INTERNAL VARIABLES
      integer p(8), q(8), irp_p, irp_q, irp_pq
      integer iSpin, iType, irrep, iLoc
      integer iSpins(2,25)
      data iSpins / 1,1 , 2,2 , 1,1 , 2,2 ,
     &              1,1 , 2,2 , 1,1 , 2,2 ,
     &              1,1 , 2,2 , 1,2 , 2,1 ,
     &              1,2 , 1,2 , 1,2 , 1,1 ,
     &              2,2 , 1,2 , 1,1 , 2,2 ,
     &              1,1 , 2,2 , 2,1 , 2,1 ,
     &              2,1 /
      character*4 szPacks(25)
      data szPacks / 'PACK', 'PACK', 'PACK', 'PACK',
     &               'PCK2', 'PCK2', 'PCK2', 'PCK2',
     &               'FULL', 'FULL', 'FULL', 'FULL',
     &               'FULL', 'FULL', 'FULL', 'FULL',
     &               'FULL', 'FULL', 'FULL', 'FULL',
     &               'FULL', 'FULL', 'FULL', 'FULL',
     &               'FULL' /
      character*1 czTypes(2,25)
      data czTypes / 'V','V' , 'V','V' , 'O','O' , 'O','O' ,
     &               'V','V' , 'V','V' , 'O','O' , 'O','O' ,
     &               'V','O' , 'V','O' , 'V','O' , 'V','O' ,
     &               'V','V' , 'O','O' , 'V','V' , 'O','V' ,
     &               'O','V' , 'O','V' , 'V','V' , 'V','V' ,
     &               'O','O' , 'O','O' , 'V','V' , 'O','O' ,
     &               'O','V' /

c COMMON BLOCKS
#include "symloc.com"
#include "syminf.com" /* for nirrep */
#include "sym.com" /* for pop and vrt */

c INLINE FUNCTIONS
      integer nnm1o2, nnp1o2, i
      nnm1o2(i)=(i*(i-1))/2
      nnp1o2(i)=(i*(i+1))/2

c ----------------------------------------------------------------------

c   o make sure the orbital vectors have been initialized
      i = 0
      do iSpin = 1, 2
      do irrep = 1, nirrep
         i = i + pop(irrep,iSpin) + vrt(irrep,iSpin)
      end do
      end do
      if (i.eq.0) return

c   o loop over pair types
      do iType = 1, 25

c      o load the proper number of orbitals into left (p) and right (q) spaces
         if (czTypes(1,iType).eq.'O') then
            call icopy(nirrep,pop(1,iSpins(1,iType)),1,p,1)
         else
            call icopy(nirrep,vrt(1,iSpins(1,iType)),1,p,1)
         end if
         if (czTypes(2,iType).eq.'O') then
            call icopy(nirrep,pop(1,iSpins(2,iType)),1,q,1)
         else
            call icopy(nirrep,vrt(1,iSpins(2,iType)),1,q,1)
         end if

c      o loop over total symmetry
         do irp_pq = 1, nirrep
            iLoc = 1
            if (szPacks(iType).eq.'FULL') then
               do irp_q = 1, nirrep
                  irp_p = dirprd(irp_q,irp_pq)
                  isymoff(irp_q,irp_pq,iType) = iLoc
                  iLoc = iLoc + p(irp_p)*q(irp_q)
               end do
            else if (szPacks(iType).eq.'PACK') then
               do irp_q = 1, nirrep
                  irp_p = dirprd(irp_q,irp_pq)
                  if (irp_p.lt.irp_q) then
                     isymoff(irp_q,irp_pq,iType) = iLoc
                     iLoc = iLoc + p(irp_p)*q(irp_q)
                  else if (irp_p.eq.irp_q) then
                     isymoff(irp_q,irp_pq,iType) = iLoc
                     iLoc = iLoc + nnm1o2(p(irp_p))
                  end if
               end do
            else if (szPacks(iType).eq.'PCK2') then
               do irp_q = 1, nirrep
                  irp_p = dirprd(irp_q,irp_pq)
                  if (irp_p.lt.irp_q) then
                     isymoff(irp_q,irp_pq,iType) = iLoc
                     iLoc = iLoc + p(irp_p)*q(irp_q)
                  else if (irp_p.eq.irp_q) then
                     isymoff(irp_q,irp_pq,iType) = iLoc
                     iLoc = iLoc + nnp1o2(p(irp_p))
                  end if
               end do
            end if
c        end do irp_pq = 1, nirrep
         end do

c     end do iType = 1, 25
      end do

      return
c     end subroutine aces_com_symloc
      end

