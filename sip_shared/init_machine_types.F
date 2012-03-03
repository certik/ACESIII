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
      subroutine init_machine_types()
      implicit none
      include 'machine_types.h'
      include 'parallel_info.h'

      integer iii(2)
      double precision x(2)

      intsize = loc(iii(2)) - loc(iii(1))
      bytes_per_double = loc(x(2)) - loc(x(1))
#ifdef ACES64
      aces64 = .true.
      aces_ints_per_double = bytes_per_double/8
#else
      aces64 = .false.
      aces_ints_per_double = bytes_per_double/intsize
#endif

#ifdef __crayx1
      mpi_int_ratio = 1
#else
      mpi_int_ratio = intsize / 4
#endif

#ifdef NO_MPI_IO
      mpi_io_support = .false.
#else
      mpi_io_support = .true.
#endif
      return
      end
