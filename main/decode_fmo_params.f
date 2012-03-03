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
      subroutine decode_fmo_params()
      implicit none
      include 'fmo.h'
      integer max_table, ntable
      parameter (max_table = 1000)

      character*80 table(max_table)
      common /parameter_table/ntable, table

      integer i, j, k, ii, ifmo, n, nvals, str_trimlen

c----------------------------------------------------------------------------
c   Search parameter table for the FMO keyword.
c----------------------------------------------------------------------------
  
      ifmo = 1
      nfmo = 0

      do i = 1, ntable
         if (table(i)(1:3) .eq. 'FMO') then

c----------------------------------------------------------------------------
c   Beginning of FMO array.
c----------------------------------------------------------------------------

            k = i
            n = str_trimlen(table(i))
            
c---------------------------------------------------------------------------
c   Find the '='.
c---------------------------------------------------------------------------

            do j = 1, len(table(i))
               if (table(i)(j:j) .eq. '=') then

c---------------------------------------------------------------------------
c   Decode the 1st FMO line.
c---------------------------------------------------------------------------

                  call c_decode_csv_integer(table(i)(j+1:n) // char(0),
     *                       fmo(ifmo), nvals)
                  go to 100
               endif
            enddo  

            print *,'Missing "=" on FMO parameter.'
            call abort_job()

  100       continue 
            if (nvals .gt. 0) then
               ifmo = ifmo + nvals
               nfmo = nfmo + nvals

c----------------------------------------------------------------------------
c   Decode remaining lines until a blank line is encountered (or end 
c   of table).    
c----------------------------------------------------------------------------

               do ii = k+1,ntable
                  n = str_trimlen(table(ii))
                  if (n .gt. 0) then
                     call c_decode_csv_integer(table(ii)(1:n)//char(0),
     *                       fmo(ifmo), nvals)
                     if (nvals .gt. 0) then
                        ifmo = ifmo + nvals
                        nfmo = nfmo + nvals
                     else
                        go to 200
                     endif 
                  else
                     go to 200
                  endif    ! n .gt. 0
               enddo    ! do ii
            endif    ! nvals .gt. 0   

         endif
      enddo
 
  200 continue
      return
      end

