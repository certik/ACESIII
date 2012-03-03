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
      SUBROUTINE IDNMAT(A,N,NTIME)
C
C SUBROUTINE RETURNS THE N X N IDENTITY MATRIX.  NTIME UTILITY IS
C  USED IN SYMUNQ.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision a(*)
      if (n.eq.3) then
         a(1) = 1.d0
         a(2) = 0.d0
         a(3) = 0.d0
         a(4) = 0.d0
         a(5) = 1.d0
         a(6) = 0.d0
         a(7) = 0.d0
         a(8) = 0.d0
         a(9) = 1.d0
      else
         do i = 1, n*n
            a(i) = 0.d0
         end do
         ndx = 1
         do i = 1, n
            a(ndx) = 1.d0
            ndx = ndx + n + 1
         end do
      end if
      NTIME=NTIME+1
      RETURN
      END
