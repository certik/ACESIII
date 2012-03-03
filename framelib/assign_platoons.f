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
      subroutine assign_platoons(companies, nprocs, platoons)
c--------------------------------------------------------------------------
c   Assigns platoons according to the c_table specifications and previous
c   company assignments.
c--------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
      include 'company.h'
      
      integer nprocs
      integer companies(nprocs), platoons(nprocs)
      integer i, j, company, platoon_id, nw, np

      do i = 1, nprocs
         platoons(i) = 0
      enddo

      do i = 1, max_company
         company = c_table(i,c_company_id)
         if (company .ne. MPI_UNDEFINED) then 
            platoon_id = c_table(i,c_platoon_id)
            if (platoon_id .eq. 0) go to 100
            nw         = c_table(i,c_nwrkr)

c----------------------------------------------------------------------
c   Find the next "nw" workers in the company.
c----------------------------------------------------------------------
            
            np = 0
            do j = 1, nprocs
               if (companies(j) .eq. company) then
                  if (platoons(j) .eq. 0) then
                     platoons(j) = platoon_id
                     np          = np + 1
                     if (np .eq. nw) go to 100
                  endif 
               endif
            enddo

            print *,'Error in platoon configuration.'
            print *,'Company ',company,' platoon ',platoon_id,
     *               ' mapped ',np,' workers, needed ',nw
            call abort_job()

  100       continue
         endif
      enddo

      return
      end
