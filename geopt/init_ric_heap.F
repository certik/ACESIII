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

c This routine allocates a heap for processing RICs and creates pointers
c into the heap for various arrays.

      subroutine init_ric_heap
      implicit none

#include "machsp.com"
#include "ric_heap.com"
#include "mxatms.par"
#include "coord.com"

      integer iNeed 
      integer*8 ndx
      integer s_DerBMat, s_BMat, s_GMat, s_BTGInv
      integer ierr
      integer i,sheap_flag
      integer totrednco

c   o don't allocate memory more than once
      if (heapptr.ne.0) return

      iNeed = 0
      CALL IGETREC(20, 'JOBARC', 'REDNCORD', 1, TOTREDNCO)

c   o DerBMat: The derivative of the B-matrix for transforming the
c              Hessian between RICs and Cartesians.
C      s_DerBMat = 9*mxatms*mxatms*maxredunco
      s_DerBMat = 9*natoms*natoms*totrednco
      iNeed = iNeed + s_DerBMat

c   o BMat: The B-matrix for RIC/Cartesian transformations.
      s_BMat = 9*totrednco*natoms
c      s_BMat = 9*mxatms*mxatms
      iNeed = iNeed + s_BMat

c   o GMat: The G-matrix for RIC/Cartesian gradient transformations.
      s_GMat = totrednco*max(totrednco,3*natoms)
C      s_GMat = 9*mxatms*mxatms
      iNeed = iNeed + s_GMat

c   o BTGInv: inv(trans(B)*G) required for transforming the Hessian.
      s_BTGInv = 3*natoms*totrednco
C      s_BTGInv = 3*mxatms*maxredunco
      iNeed = iNeed + s_BTGInv

c      ndx = iintfp*iNeed
c      print *,'INIT_RIC_HEAP: iintfp, iNeed, ndx ',iintfp, iNeed, ndx
c      call aces_malloc(ndx,dRICHeap,z_RICHeap, heapptr)
c      if (heapptr .eq. 0) then
c         print *, '@INIT_RIC_HEAP: Failed to allocate memory.'
c         print *, '                need ',iNeed/1024/1024,' MB'
c         call errex
c      end if

c      print *,'AFTER ACES_MALLOC: heapptr = ',heapptr,' z_RICHeap ',
c     *   z_RICHeap
c      z_RICHeap = (z_RICHeap+iintfp-1)/iintfp

#ifdef ALTIX
      sheap_flag = 1
#else
      sheap_flag = 0
#endif

      call mem_alloc(dRICHeap, iNeed, ifltln,
     *               z_RICHeap, .true., ierr)
      if (ierr .ne. 0) then
         print *,'ERROR in INIT_RIC_HEAP: Cannot alloc ',iNeed,
     *      ' words of double precision data'
         call errex() 
      endif

      ndx = z_RICHeap
#ifdef _DEBUG
      call c_memset(dRICHeap(z_RICHeap),255,ifltln*iNeed)
#endif

c   o DerBMat location
      z_DerBMat = ndx
      ndx = ndx + s_DerBMat

c   o BMat location
      z_BMat = ndx
      ndx = ndx + s_BMat

c   o GMat location
      z_GMat = ndx
      ndx = ndx + s_GMat

c   o BTGInv location
      z_BTGInv = ndx
      ndx = ndx + s_BTGInv

c      print *, '@INIT_RIC_HEAP: RIC HEAP INDICES'
c      print *, '                z_RICHeap = ',z_RICHeap
c      print *, '                z_DerBMat = ',z_DerBMat
c      print *, '                z_BMat    = ',z_BMat
c      print *, '                z_GMat    = ',z_GMat
c      print *, '                z_BTGInv  = ',z_BTGInv

      return
c     end subroutine init_ric_heap
      end

