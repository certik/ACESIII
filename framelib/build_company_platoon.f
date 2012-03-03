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
      subroutine build_company_platoon(id, platoon_id, nwrkr, wrkr_mem,
     *                   iocompany_flag, sial_program)
c-------------------------------------------------------------------------
c   
c   Configures the data for a company, and saves it in a table.
c
c   Arguments:
c	id			Integer identifier for the company.
c       platoon_id              Integer identifier of the platoon within the
c                               company.
c	nwrkr			Number of workers in the company.
c	wrkr_mem		Memory (in Mbytes) for each worker.
c       iocompany_flag          Logical flag indicating whether this is 
c                               the i/o company.
c       sial_program            SIAL program to be run by all workers within 
c                               the company
c
c-------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'company.h'

      integer id, platoon_id, nwrkr, wrkr_mem
      integer i, ierr
      logical iocompany_flag
      character*(*) sial_program

      save

c--------------------------------------------------------------------------
c   Find the first unused entry in c_table.
c--------------------------------------------------------------------------

      do i = 1, max_company
         if (c_table(i,c_company_id) .eq. MPI_UNDEFINED) then

c--------------------------------------------------------------------------
c   Store the arguments in the table's entry.
c--------------------------------------------------------------------------

           c_table(i,c_company_id) = id
           c_table(i,c_platoon_id) = platoon_id
           c_table(i,c_nwrkr) = nwrkr
           c_table(i,c_wrkr_mem) = wrkr_mem
           if (iocompany_flag) then
              c_table(i,c_iocompany) = 1
           else
              c_table(i,c_iocompany) = 0
           endif
           company_sial_prog(i) = sial_program

           return
         endif
      enddo
 
      print *,'Error: Limit of ', max_company,
     *        ' companies has been exceeded.'
      call mpi_abort(mpi_comm_world, 1, ierr)
      return
      end
