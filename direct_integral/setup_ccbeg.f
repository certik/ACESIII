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
      subroutine setup_ccbeg(alpha, ixalpha, pcoeff, ixpcoeff, 
     *                       ncfps, npfps, nshells, ccbeg, ccend,
     *                       indx_cc)
c---------------------------------------------------------------------------
c   Computes the beginning and ending non-zero contraction coefficients,
c   used to gain efficiency in the ERD integral package.
c---------------------------------------------------------------------------
      implicit none
      
      integer nshells
      integer ncfps(*), npfps(*), ixalpha(*), ixpcoeff(*)
      integer nalpha, npcoeff
      double precision alpha(*), pcoeff(*)
      integer ccbeg(*), ccend(*), indx_cc(nshells)
      integer i, j, ishell
      integer k, l, icc

      icc = 1
      do ishell = 1, nshells
         indx_cc(ishell) = icc
         do k = 1, ncfps(ishell)
           l = (k-1)*npfps(ishell) + 1
           do i = 1, npfps(ishell)
              if (pcoeff(ixpcoeff(ishell)+l-1) .ne. 0.d0) then
                 ccbeg(icc) = i
                 go to 10
              endif
              l = l + 1
           enddo
   10      continue

           l = k*npfps(ishell)
           do i = npfps(ishell), 1, -1
              if (pcoeff(ixpcoeff(ishell)+l-1) .ne. 0.d0) then
                 ccend(icc) = i
                 go to 20
              endif
              l = l - 1
           enddo
   20      continue 

           icc = icc + 1
         enddo   ! k
      enddo      ! ishell

      return
      end  
