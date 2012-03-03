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
      subroutine init_company
c-------------------------------------------------------------------------
c   Initializes the table of company requirements.
c-------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'company.h'

      integer i, j

      do i = 1, max_company
         do j = 1, lcompany_entry
            c_table(i,j) = 0
         enddo

         c_table(i,c_company_id) = MPI_UNDEFINED
         company_sial_prog(i) = ' '
      enddo

      return
      end
