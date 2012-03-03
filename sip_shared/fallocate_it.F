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
      subroutine fallocate_it(sheap_flag, nw, ierr)
c--------------------------------------------------------------------------
c   Allocates memory from the shared heap, and sets up pointers to access
c   it.  The "ishptr" (ishared_heap array) is used to refer to the memory
c   as integers, and "dshptr" (dshared_heap array) is used to refer to
c   the shared heap as a double precision array.
c-------------------------------------------------------------------------

      implicit none
      integer nw
      integer ierr
      integer sheap_flag
#ifdef ALTIX
      include 'parallel_info.h'
      include 'machine_types.h'
      include 'sheap.h'
      common /addrchk/last_sheap_addr
      integer*8 last_sheap_addr
      
      ishptr = malloc(nw*intsize)
      if (ishptr .eq. 0) then
         ierr = -1
      else
         ierr = 0
      endif

      dshptr = ishptr
      i8shptr = ishptr
#else
      print *,'fallocate_it call is illegal on this system!'
      call abort_job()
#endif
      return
      end

