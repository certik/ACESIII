C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.

c This routine fills in the /int_gen_com/ common block and prints a
c summary of the basis set.


      subroutine int_gen_init
      implicit none

      integer i, j, n, iCenter, iShell, iOffset, iTmp
      integer previous_end

      include 'int_gen_parms.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer          temp(nshells)
      double precision xtemp(nshells)

      double precision dbuf(1)
#ifdef ALTIX
      pointer (dptr, dbuf)
#else
      common dbuf
#endif
      integer*8        iii
      integer          icomp

      integer nCarts, nSphrs, l

#ifdef ALTIX
      dptr = dshptr
#endif
c   o /int_gen_com/ nFpS(*), ivShlOff(*), nContAOs
      if (iSpherical.ne.0) then
         iOffset = 0
         do iShell = 1, nShells
            ivShlOff(iShell) = iOffset
            nSphrs = 2*ivAngMom(iShell)+1 
            nFpS(iShell)     = nSphrs*nCFpS(iShell)
            iOffset          = iOffset + nFpS(iShell)
         end do
         nContAOs = iOffset
      else
         iOffset = 0
         do iShell = 1, nShells
            ivShlOff(iShell) = iOffset
            nCarts = (ivAngMom(iShell)+1)*(ivAngMom(iShell)+2)/2 
            nFpS(iShell)     = nCarts*nCFpS(iShell)
            iOffset          = iOffset + nFpS(iShell)
         end do
         nContAOs = iOffset
      end if

c---------------------------------------------------------------------------
c   Perform better fit of shells to segment size for better performance.
c---------------------------------------------------------------------------

         if (intpkg .eq. gamess_package .or.
     *       intpkg .eq. gamess_derivative_package) go to 1000
         if (aat_database(1:1) .eq. ' ') then
            call binpack(nshells, nfps, sip_mx_ao_segsize, ixshells)
         else

c---------------------------------------------------------------------------
c   No sorting of shells.  Shells are processed in their original order.
c---------------------------------------------------------------------------

            do i = 1, nshells
               ixshells(i) = i
            enddo
         endif

c---------------------------------------------------------------------------
c   Reorder data according to the index array returned by binpack.
c---------------------------------------------------------------------------

         call reorder_coefficients(dbuf(iscfa), nbasis, dbuf(iepsa),
     *                  nfps, ixshells, nshells)
         if (nscfb .gt. 0) then
            call reorder_coefficients(dbuf(iscfb), nbasis,dbuf(iepsb),
     *                  nfps, ixshells, nshells)
         endif

         if (iuhf .eq. 2) then    ! ROHF
            call reorder_rohf(dbuf(ifockrohfa), nbasis, nfps, ixshells,
     *                       nshells)
            call reorder_rohf(dbuf(ifockrohfb), nbasis, nfps, ixshells,
     *                       nshells)
         endif

         do i = 1, nshells
            temp(i) = nfps(i)
         enddo

         do i = 1, nshells
            nfps(i) = temp(ixshells(i))
         enddo

         do i = 1, nshells
            temp(i) = ixalpha(i)
         enddo

         do i = 1, nshells
            ixalpha(i) = temp(ixshells(i))
         enddo

         do i = 1, nshells
            temp(i) = ixpcoef(i)
         enddo

         do i = 1, nshells
            ixpcoef(i) = temp(ixshells(i))
         enddo

         do i = 1, nshells
            temp(i) = ncfps(i)
         enddo

         do i = 1, nshells
            ncfps(i) = temp(ixshells(i))
         enddo

         do i = 1, nshells
            temp(i) = npfps(i)
         enddo

         do i = 1, nshells
            npfps(i) = temp(ixshells(i))
         enddo

         do i = 1, nshells
            temp(i) = ivangmom(i)
         enddo

         do i = 1, nshells
            ivangmom(i) = temp(ixshells(i))
         enddo

         do i = 1, nshells
            temp(i) = atom(i)
         enddo

         do i = 1, nshells
            atom(i) = temp(ixshells(i))
         enddo

c--------------------------------------------------------------------------
c   Reorder the shell x,y,z coordinates.
c--------------------------------------------------------------------------

         do icomp = 1, 3
            iii = master_icoord + icomp - 1

c--------------------------------------------------------------------------
c   Reorder component "icomp" of the coordinate array.
c--------------------------------------------------------------------------

            do i = 1, nshells
               xtemp(i) = dbuf(iii+3*(i-1))
            enddo

            do i = 1, nshells
               dbuf(iii+3*(i-1)) = xtemp(ixshells(i))
            enddo
         enddo

1000  continue

      previous_end = 0
      do iShell = 1, nShells
         end_nfps(iShell) = previous_end + nfps(iShell) 
         previous_end     = end_nfps(iShell)
      enddo


      if (intpkg .eq. flocke_package) then
         
c---------------------------------------------------------------------------
c   Shell offsets are different for Flocke package.
c---------------------------------------------------------------------------

         ivShlOff(1) = 0
         do iShell = 2, nShells
            ivShlOff(iShell) = end_nfps(iShell-1)
         enddo 
      endif
 
      return
c     end subroutine int_gen_init
      end
