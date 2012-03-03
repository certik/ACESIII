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
      double precision function nn_repulsion_energy(coord, charge, 
     *                           natoms)
c-------------------------------------------------------------------------
c   Returns the nuclear-nuclear repulsion energy term.
c-------------------------------------------------------------------------
      implicit none
      integer natoms
      double precision coord(3,natoms), charge(natoms)

      integer i, j
      double precision sum, x,y,z,r

      sum = 0.
      do i = 1, natoms
      do j = i + 1, natoms
         x = coord(1,i)-coord(1,j)
         y = coord(2,i)-coord(2,j)
         z = coord(3,i)-coord(3,j)
         r = dsqrt(x*x + y*y + z*z)
         sum = sum + charge(i) * charge(j) / r
      enddo
      enddo

      nn_repulsion_energy = sum
      return
      end
