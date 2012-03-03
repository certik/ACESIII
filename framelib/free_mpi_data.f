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
      subroutine free_mpi_data()
c---------------------------------------------------------------------------
c   Frees the communicators and groups set up for a previous SIAL execution.
c---------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'parallel_info.h'
      include 'company.h'
      integer pst_get_company_comm

      integer comm
      integer group
      integer ierr

      comm = pst_get_company_comm(me)
      call mpi_comm_free(comm, ierr)
     
      return
      end
