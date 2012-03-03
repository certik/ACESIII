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
      subroutine handle_hessian(hess, nnhess, ncenters, write_hess)
      implicit none

      integer ncenters
      double precision hess(3*ncenters,3*ncenters)
      double precision NNhess(3*ncenters,3*ncenters)
      logical write_hess
 
      integer i, j, ix, jx, ihess, jhess, ihess0, jhess0
      integer nunique
      double precision num, sum, res, unique 

c--------------------------------------------------------------------------
c   Require the symmetry of the Hessian to be the same as that of the
c   NN hessian. --> Zero out 'small' residuals due to integral accuracy,
c   etc.
c--------------------------------------------------------------------------

      do i = 1, Ncenters
      do j = 1, Ncenters
         do ix = 1, 3
         do jx = 1, 3
           ihess = (i-1)*3 + ix
           jhess = (j-1)*3 + jx

c           if (dabs(NNhess(jhess,ihess)) .lt. 1.0d-10) then
c              hess(jhess,ihess) = 0.0d0
c           else
              hess(jhess,ihess) = hess(jhess,ihess)
     *                          + NNhess(jhess,ihess)
c           endif
         enddo
         enddo
      enddo
      enddo

c--------------------------------------------------------------------------
c   Print Hessian data to stdout.
c--------------------------------------------------------------------------

      if (write_hess) then
         write(6,*) ' '
         write(6,*) ' Final Hessian data '

         do i = 1, Ncenters
         do j = 1, Ncenters

            write(6,*) ' ATOM A', i, ' ATOM B', j
            write(6,*) ' '

            write(6,33) hess((j-1)*3+1,(i-1)*3+1),
     *                   hess((j-1)*3+2,(i-1)*3+1),
     *                   hess((j-1)*3+3,(i-1)*3+1)

            write(6,33) hess((j-1)*3+1,(i-1)*3+2),
     *                   hess((j-1)*3+2,(i-1)*3+2),
     *                   hess((j-1)*3+3,(i-1)*3+2)

            write(6,33) hess((j-1)*3+1,(i-1)*3+3),
     *                   hess((j-1)*3+2,(i-1)*3+3),
     *                   hess((j-1)*3+3,(i-1)*3+3)

         enddo ! i
         enddo ! j
      endif
33    format(3F16.8)

c--------------------------------------------------------------------------
c   Make sure equivalent Hessian elements are exactly the same
c--------------------------------------------------------------------------

      do ihess0 = 1, 3*ncenters
      do jhess0 = 1, 3*ncenters 
         unique = nnhess(jhess0, ihess0)

         num = 0.d0
         sum = 0.d0
         do ihess = 1, 3*ncenters
         do jhess = 1, 3*ncenters
            if (nnhess(jhess,ihess) .eq. unique) then
               if (dabs(hess(jhess,ihess)-hess(jhess0,ihess0)) .lt. 
     *                           1.d-05) then
                  num = num + 1.d0
                  sum = sum + hess(jhess,ihess)
               endif
            endif
         enddo
         enddo

         res = sum / num
         do ihess = 1, 3*ncenters
         do jhess = 1, 3*ncenters
            if (nnhess(jhess,ihess) .eq. unique) then
               if (dabs(hess(jhess,ihess)-hess(jhess0,ihess0)) .lt.
     *                           1.d-05) then
                  hess(jhess,ihess) = res
               endif
            endif
         enddo
         enddo
      enddo
      enddo

      return
      end
