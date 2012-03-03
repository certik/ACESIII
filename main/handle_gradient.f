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
      subroutine handle_gradient()
      implicit none 
      include 'int_gen_parms.h'

      integer unique(3*Ncenters), equal(3*Ncenters)

      integer ns
      integer iatom, jatom, jcomp, icomp
      integer i, j, xx, yy, zz
      double precision diff, avg, z

        do i = 1, 3*Ncenters 
           unique(i) = 0 
        enddo 

        do i = 1, 3*Ncenters 

c----------------------------------------------------------------------------
c   Find the next component of the nuclear-nuclear gradient for which we 
c   have not got a match yet.
c----------------------------------------------------------------------------

           if (unique(i) .eq. 0) then 
              icomp = mod(i-1,3)+1
              iatom = (i-1)/3+1
              z = dabs(NNgrad(icomp,iatom))

              ns = 1
              avg = dabs(gradient_data(i))

c----------------------------------------------------------------------------
c   Find and average all matching components of the nuclear-nuclear gradient.
c----------------------------------------------------------------------------   

              do j = i+1, 3*Ncenters 
                 jcomp = mod(j-1,3)+1
                 jatom = (j-1)/3+1
                 diff = z - dabs(NNgrad(jcomp, jatom)) 
                 if (dabs(diff) .lt. 1.0d-10) then 
                    unique(j)  = 1    ! flag as matching
                    equal(j)   = 1    ! flag for averaging
                    ns = ns + 1
                    avg = avg + dabs(gradient_data(j))
                 else
                    equal(j) = 0
                 endif 
              enddo ! j                

c----------------------------------------------------------------------------
c   Replace the gradient value with the average, retaining the sign of the 
c   original gradient value.
c----------------------------------------------------------------------------

              avg = avg / ns

              if (gradient_data(i) .gt. 0)
     &           gradient_data(i) = avg
              if (gradient_data(i) .lt. 0)
     &           gradient_data(i) = -avg

c-----------------------------------------------------------------------------
c   Replace all the gradient values which were used in the averaging process.
c-----------------------------------------------------------------------------

              do j = i + 1, 3*NCenters
                 if (equal(j) .eq. 1) then
                    if (gradient_data(j) .gt. 0)
     &                  gradient_data(j) = avg
                    if (gradient_data(j) .lt. 0)
     &                  gradient_data(j) = -avg
                 endif
              enddo

           endif 
        enddo ! i                

c--------------------------------------------------------------------------
c    Check for zero contribution 
c--------------------------------------------------------------------------

        do i = 1, 3*Ncenters 
           icomp = mod(i-1,3)+1
           iatom = (i-1)/3+1  
           if (dabs(NNgrad(icomp,iatom)) .lt. 1.0d-10) 
     &        gradient_data(i) = 0.0d0 
        enddo ! i 

c--------------------------------------------------------------------------
c   Write out symmetrized gradient data 
c--------------------------------------------------------------------------

        write(6,*) ' ' 
        write(6,*) ' Symmetrized gradient information ' 
        do i = 1, Ncenters 
           xx = 3*(i-1) + 1  
           yy = 3*(i-1) + 2  
           zz = 3*(i-1) + 3  
           write(6,*) i, gradient_data(xx),  
     &                   gradient_data(yy),  
     &                   gradient_data(zz)  
        enddo ! i 

      return
      end
