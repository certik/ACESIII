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
      double precision function nuclear_nuclear_repulsion_energy(natoms)
c---------------------------------------------------------------------------
c   Returns the nuclear-nuclear repulsion energy of a system of atoms.
c   The charge and geometry data is found in the NUCLEAR common block.
c---------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'

      integer natoms
      integer i, j
      double precision vnn, x, y, z, r

      VNN = 0.0D0
      DO I = 1, natoms
         DO J = I + 1, natoms
            X = acenter(I,1) - acenter(J,1)
            Y = acenter(I,2) - acenter(J,2)
            Z = acenter(I,3) - acenter(J,3)
            R = DSQRT(X**2 + Y**2 + Z**2)
            VNN = VNN + charge(i)*charge(j)/R
         ENDDO ! J
      ENDDO ! I

      nuclear_nuclear_repulsion_energy = vnn
      return
      end 
