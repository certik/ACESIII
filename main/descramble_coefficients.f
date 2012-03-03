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
      subroutine descramble_coefficients(ca, na, epsa, nfps, 
     *                ixshells, nshells)
c--------------------------------------------------------------------------
c   Unscrambles the transformation coefficients and eigenvalues 
c   prior to writing them out.  After the descramble procedure, the
c   data is back into VMOL order, and is usable by other ACES member
c   executables.
c--------------------------------------------------------------------------

      implicit none
      integer na, nshells
      integer nfps(nshells)
      integer ixshells(nshells)
      double precision ca(na,na)
      double precision epsa(na)
      double precision xa(na)
      integer ixa(na)
      integer nfps_orig(nshells)

      integer i, j, istart, ishell, n

      if (na .gt. 0) then

c-------------------------------------------------------------------------
c   Descramble the number of functions per shell.
c-------------------------------------------------------------------------

         do i = 1, nshells
            nfps_orig(ixshells(i)) = nfps(i)
         enddo

c--------------------------------------------------------------------------
c   Build a basis function index array from the shell index array.
c--------------------------------------------------------------------------

         istart = 1
         do i = 1, nshells
            ishell = ixshells(i)

c---------------------------------------------------------------------------
c   Calculate the starting basis function.
c---------------------------------------------------------------------------

            n = 0
            do j = 1, ishell-1
               n = n + nfps_orig(j)
            enddo

c----------------------------------------------------------------------------
c   Store the next "nfps_orig(ishell)" indices in the index array.
c----------------------------------------------------------------------------

            do j = 1, nfps_orig(ishell)
               ixa(istart) = n + j
               istart = istart + 1
            enddo
         enddo

         do j = 1, na

c--------------------------------------------------------------------------
c   Save column "j".
c--------------------------------------------------------------------------

            do i = 1, na
               xa(i) = ca(i,j)
            enddo

            do i = 1, na
               ca(ixa(i),j) = xa(i)
            enddo

         enddo

c         do i = 1, na
c            xa(i) = epsa(i)
c         enddo

c         do i = 1, na
c            epsa(ixa(i)) = xa(i)
c         enddo
      endif

      return
      end
