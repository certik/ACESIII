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
      integer*8 function get_index_from_base(base, anchor, type)
      implicit none
      include 'machine_types.h'

      double precision x(2)
      integer ix(2)
      integer size

      integer anchor(*)
      integer*8 base
      integer*8 c_loc64
      integer*8 anchor_addr
      integer*8 ixx
      integer type
      integer*8 unmask_addr
      integer*8 ubase, uanchor

      if (type .eq. 1) then
         size = intsize
      else if (type .eq. 2) then
         size = bytes_per_double
      else if (type .eq. 3) then
         size = 8   ! 64 bit integer
      else
         print *,'Invalid size in get_index_from_base: size ',size
      endif

      ixx = 1
      anchor_addr = c_loc64(anchor, ixx, 1)
      ubase = unmask_addr(base)
      uanchor = unmask_addr(anchor_addr)
c      print *,'Task ',my_pe(),' GET_INDEX: base ',ubase,
c     *     ' anchor_addr ',uanchor,
c     *     ' size ', size,' return index ',
c     *          (ubase-uanchor)/size+1
     
      get_index_from_base = (ubase - uanchor)/size+1
      return
      end
      
