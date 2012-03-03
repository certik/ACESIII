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
      subroutine set_geometry(iatom, coord, atom_charge)
c--------------------------------------------------------------------------
c   Stores the atomic charge and geometry data.
c--------------------------------------------------------------------------

      implicit none
      include 'int_gen_parms.h'

      integer iatom
      double precision coord(3)
      double precision atom_charge  
      integer i

      if (iatom .gt. max_centers .or. iatom .le. 0)  then
         print *,'Error in set_geometry: Attempt to set iatom ',
     *     iatom,' but max_centers is ',max_centers
         call abort_job()
      endif

      charge(iatom) = atom_charge
      do i = 1, 3
         acenter(iatom,i) = coord(i)
      enddo   

      return
      end

      subroutine get_geometry(iatom, atom_charge, coord)
c-------------------------------------------------------------------------
c   Retrieves atomic coordinate and charge info for atom number "iatom".
c-------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'

      integer iatom
      double precision atom_charge, coord(3)

      integer i

      if (iatom .gt. max_centers .or. iatom .le. 0)  then
         print *,'Error in get_geometry: Attempt to get iatom ',
     *     iatom,' but max_centers is ',max_centers
         call abort_job()
      endif

      atom_charge = charge(i)
      do i = 1, 3
         coord(i) = acenter(iatom,i)
      enddo 

      return
      end

      integer function get_ncenters()
c--------------------------------------------------------------------------
c   Retrieves the number of atomic centers.
c--------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'

      get_ncenters = ncenters
      return
      end
