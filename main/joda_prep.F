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
      subroutine joda_prep()
      implicit none
      include 'machine_types.h'

      integer i, natom_jobarc, iordgp

c--------------------------------------------------------------------------
c   Prepare JOBARC for joda processing at the end of a geometry opt. step. 
c--------------------------------------------------------------------------

      Call igetrec (20, 'JOBARC', 'NATOMS  ', 1, natom_jobarc)
      Call igetrec (20, 'JOBARC', 'COMPORDR', 1, iordgp) 

      call built_cmp2zmat(natom_jobarc, iordgp)
      return
      end
