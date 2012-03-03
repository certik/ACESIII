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
      subroutine compute_integrals_wrapper(bval, eval, nind, x)
c-----------------------------------------------------------------------------
c   Wrapper routine to provide an interface for C and C++ code to 
c   use the "compute_integrals" subroutine.
c-----------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
#ifdef ALTIX
      include 'sheap.h'
#endif
      integer nind
      integer bval(nind), eval(nind)
      double precision x(*)

      integer i, nsend, maxblk
      integer*8 integral_scr
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif

      if (nind .ne. 4) then
         print *,'Error: compute_integrals_wrapper called with ',
     *       'an incorrect number of indices: ',nind
         call abort_job()
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      if (intpkg .eq. flocke_package) then
         integral_scr = dpkgscr
         maxblk  = zmax
      endif

      call compute_integrals(bval(1),eval(1), bval(2),eval(2),
     *                       bval(3),eval(3),bval(4),eval(4),
     *                 dbuf(integral_scr),maxblk, ibuf(ipkgscr),
     *                 dbuf(icoord),dbuf(ipcoeff),dbuf(ialpha),
     *                 ibuf(iccbeg), ibuf(iccend),
     *                 x, nsend)

      return
      end
