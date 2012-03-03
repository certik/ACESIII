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

cproto      subroutine dcopy(n,src,lda,dst,ldb)

cproto      subroutine trnsps_ket(nrows,ncols,obj_size,
cproto     &                      src,lda1,lda2,dst,ldb1,ldb2)
      subroutine trnsps_ket(nrows,ncols,obj_size,
     &                      src,lda1,lda2,dst,ldb1,ldb2)
      implicit none
      integer nrows, ncols, obj_size, lda1, lda2, ldb1, ldb2
      double precision src(lda1,lda2,*), dst(ldb1,ldb2,*)
      integer i, j, k
c#ifdef ASSERT
      i = 0
      if (nrows.lt.1.or.ncols.lt.1.or.obj_size.lt.1) then
         print *, '@TRNSPS_KET: Assertion failed.'
         print *, '             obj_size = ',obj_size
         print *, '             nrows    = ',nrows
         print *, '             ncols    = ',ncols
         i = 1
      end if
      if (lda1.lt.obj_size.or.ldb1.lt.obj_size) then
         print *, '@TRNSPS_KET: Assertion failed.'
         print *, '             obj_size = ',obj_size
         print *, '             lda1     = ',lda1
         print *, '             ldb1     = ',ldb1
         i = 1
      end if
      if (lda2.lt.nrows.or.ldb2.lt.ncols) then
         print *, '@TRNSPS_KET: Assertion failed.'
         print *, '             nrows = ',nrows
         print *, '             lda2  = ',lda2
         print *, '             ncols = ',ncols
         print *, '             ldb2  = ',ldb2
         i = 1
      end if
      if (i.ne.0) call c_exit(1)
c#endif
      do k = 1, nrows
      do j = 1, ncols
         do i = 1, obj_size
            dst(i,j,k) = src(i,k,j)
         end do
      end do
      end do
      return
      end

c ----------------------------------------------------------------------

cproto      subroutine trnsps_bra(nrows,ncols,times,
cproto     &                      src,lda1,lda2,dst,ldb1,ldb2)
      subroutine trnsps_bra(nrows,ncols,times,
     &                      src,lda1,lda2,dst,ldb1,ldb2)
      implicit none
      integer nrows, ncols, times, lda1, lda2, ldb1, ldb2
      double precision src(lda1,lda2,*), dst(ldb1,ldb2,*)
      integer i, j, k
c#ifdef _ASSERT
      i = 0
      if (nrows.lt.1.or.ncols.lt.1.or.times.lt.1) then
         print *, '@TRNSPS_BRA: Assertion failed.'
         print *, '             nrows = ',nrows
         print *, '             ncols = ',ncols
         print *, '             times = ',times
         i = 1
      end if
      if (lda1.lt.nrows.or.ldb2.lt.nrows) then
         print *, '@TRNSPS_BRA: Assertion failed.'
         print *, '             nrows = ',nrows
         print *, '             lda1  = ',lda1
         print *, '             ldb2  = ',ldb2
         i = 1
      end if
      if (lda2.lt.ncols.or.ldb1.lt.ncols) then
         print *, '@TRNSPS_BRA: Assertion failed.'
         print *, '             ncols = ',ncols
         print *, '             lda2  = ',lda2
         print *, '             ldb1  = ',ldb1
         i = 1
      end if
      if (i.ne.0) call c_exit(1)
c#endif
      do k = 1, times
         do j = 1, nrows
         do i = 1, ncols
            dst(i,j,k) = src(j,i,k)
         end do
         end do
      end do
      return
      end

c ----------------------------------------------------------------------

cproto      subroutine trnsps(nrows,ncols,src,lda,dst,ldb)
      subroutine trnsps(nrows,ncols,src,lda,dst,ldb)
      implicit none
      integer nrows, ncols, lda, ldb
      double precision src(lda,*), dst(ldb,*)
      integer tile
      parameter (tile=64)
      integer i, j, bi, bj
c#ifdef _ASSERT
      i = 0
      if (nrows.lt.1.or.ncols.lt.1) then
         print *, '@TRNSPS: Assertion failed.'
         print *, '         nrows = ',nrows
         print *, '         ncols = ',ncols
         i = 1
      end if
      if (lda.lt.nrows.or.ldb.lt.ncols) then
         print *, '@TRNSPS: Assertion failed.'
         print *, '         nrows = ',nrows
         print *, '         ncols = ',ncols
         print *, '         lda   = ',lda
         print *, '         ldb   = ',ldb
         i = 1
      end if
      if (i.ne.0) call c_exit(1)
c#endif
c      do j = 1, nrows
c      do i = 1, ncols
c         dst(i,j) = src(j,i)
c      end do
c      end do
c 1x1   -> 3.61
c 4x4   -> 1.18
c 8x8   -> 0.81
c 16x16 -> 0.62
c 32x32 -> 0.51
c 64x64 -> 0.42
c 96x96 -> 0.40
      do bj = 1, nrows, tile
      do bi = 1, ncols, tile
         do j = bj, min(nrows,bj+tile-1)
         do i = bi, min(ncols,bi+tile-1)
            dst(i,j) = src(j,i)
         end do
         end do
      end do
      end do
      return
      end

