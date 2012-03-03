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
c
      subroutine form_R(b_array,
     *                   c1,c2,c3,c4,c5, c6,c7,c8,c9,c10)  
c
      double precision R(11,11), rpack(121), v(11) 
      double precision b_array(10,10)
      double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      double precision temp
      integer i, ndim 
#ifdef __crayx1
      integer*8 arg64(5), iP(11)
#else
      integer iP(11)
#endif

      integer maxdim
      parameter (maxdim = 10)

c--------------------------------------------------------------------------
c   Symmetrize the b_array into the R array
c--------------------------------------------------------------------------

      do i = 1, maxdim
      do j = i, maxdim
         R(i,j) = b_array(i,j)
         if (i .ne. j) R(j,i) = b_array(i,j) 
      enddo
      enddo

      ndim = maxdim
      do i = 1, maxdim
         nzero = maxdim
         do j = 1, maxdim
            if (dabs(r(i,j)) .gt. 0.d0) then
               nzero = nzero - 1
            endif
         enddo
         if (nzero .eq. maxdim) ndim = ndim - 1
      enddo

      do i = 1, ndim
         R(i,ndim+1) = -1.0d0 
         R(ndim+1,i) = -1.0d0 
      enddo ! i 
      R(ndim+1,ndim+1) =  0.0d0 

      n = 0
      do j = 1, ndim+1
      do i = 1, ndim+1
         n = n + 1
         rpack(n) = r(i,j)
      enddo
      enddo

c
c Form the v-vector. 
c ------------------ 
c
      do i = 1, ndim
         v(i) = 0.
      enddo
      v(ndim+1) = -1.d0

c
c Solve the eigenvalue problem R * C = v for the vector C.
c -------------------------------------------------------- 
c
#ifdef __crayx1
      arg64(1) = ndim+1
      arg64(2) = 1
      arg64(3) = ndim+1
      arg64(4) = ndim+1
      arg64(5) = i
      call dgesv(arg64(1), arg64(2), rpack, arg64(3), iP, v,
     *           arg64(4), arg64(5))
      i = arg64(5)  
#else
      call dgesv(ndim+1,1,rpack,ndim+1,iP,v,ndim+1,i)  
#endif

      if (i .ne. 0) then
         print *,'Error in dgesv: info = ',i
         call abort_job()
      endif

c
c Put solutions into c for use in the SIAL program. 
c ------------------------------------------------------------------------ 
c
      do i = ndim+1, maxdim+1
         v(i) = 0.
      enddo

      c1 = v(1) 
      c2 = v(2) 
      c3 = v(3) 
      c4 = v(4) 
      c5 = v(5) 
      c6 = v(6)
      c7 = v(7)
      c8 = v(8)
      c9 = v(9)
      c10 = v(10)
 
c
      end 
c 

