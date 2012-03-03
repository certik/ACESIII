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
      integer*8 function to_bytes(nwords, type)
c-------------------------------------------------------------------------
c   Converts words to bytes.  Type = 1 --> integer, type = 2 --> float,
c   type = 3 --> double precision.
c-------------------------------------------------------------------------
      implicit none
      include 'machine_types.h'

      integer*8 nwords, type
      integer*8 len

      if (type .eq. 1) len = intsize     
c      if (type .eq. 2) len = bytes_per_double 
      if (type .eq. 3) len = bytes_per_double

      len = len * nwords
      to_bytes = len
      return
      end

