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

      subroutine det_segsizes(auto_seg_switch, nproc)
      implicit none
      integer auto_seg_switch, nproc
      integer oa_segsize, ob_segsize, va_segsize, vb_segsize

#include "int_gen_parms.h"

      integer iatom, nelec, naocc, nbocc, iShell, ierr

#ifdef ACES64
      integer*8 arg_1, arg_2, arg_3, iFlags(100), ivCharge(100)
#else
      integer   arg_1, arg_2, arg_3, iFlags(100), ivCharge(100)
#endif

      integer nCarts, nSphrs, l
      nCarts(l) = (l+1)*(l+2)/2
      nSphrs(l) = 2*l+1

c   o count the number of AOs per shell
      if (iSpherical.ne.0) then
         do iShell = 1, nShells
            nFpS(iShell) = nSphrs(ivAngMom(iShell))*nCFpS(iShell)
         end do
      else
         do iShell = 1, nShells
            nFpS(iShell) = nCarts(ivAngMom(iShell))*nCFpS(iShell)
         end do
      end if

      naocc = nalpha_occupied 
      nbocc = nbeta_occupied 

      if (auto_seg_switch .eq. 1) then
         call det_segsizes_kern(naocc,nbocc,sip_mx_occ_segsize,
     &                       nalpha_virtual, nbeta_virtual,
     &                       sip_mx_virt_segsize,
     &                       nFpS,nShells,sip_mx_ao_segsize
     &                       )
      else if (auto_seg_switch .eq. 2) then
         call find_seg(naocc, nalpha_virtual, naocc+nalpha_virtual, 
     *                 nproc, oa_segsize, va_segsize,
     *                 sip_mx_ao_segsize)
         call find_seg(nbocc, nbeta_virtual, naocc+nalpha_virtual, 
     *                 nproc, ob_segsize, vb_segsize,
     *                 sip_mx_ao_segsize)
         sip_mx_occ_segsize = max(oa_segsize, ob_segsize)
         sip_mx_virt_segsize = max(va_segsize, vb_segsize)
      else 
         print *,'Error: Invalid value for auto_seg_switch ',
     *      'in det_segsizes: ',auto_seg_switch
         call abort_job()
      endif

      return
      end

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------

      subroutine det_segsizes_kern(nocca,noccb,seg_occ,
     &                             nvrta,nvrtb,seg_vrt,
     &                             nFpS,nShells,seg_ao
     &                             )
      implicit none

c ARGUMENTS
      integer nocca, noccb, nvrta, nvrtb, nShells
      integer nFpS(nShells), seg_occ, seg_vrt, seg_ao

c PARAMETERS
      integer mn_seg, mx_seg
      parameter (mn_seg=4)  ! 2kB block
c      parameter (mx_seg=16) ! 0.5MB block
c      parameter (mx_seg=19) ! 1MB block
c      parameter (mx_seg=22) ! 2MB block
c      parameter (mx_seg=26) ! 4MB block
      parameter (mx_seg=32) ! 8MB block
c      parameter (mx_seg=35) ! 12MB block
c      parameter (mx_seg=38) ! 16MB block
c      parameter (mx_seg=40) ! 20MB block
c      parameter (mx_seg=42) ! 24MB block
c      parameter (mx_seg=43) ! 28MB block
c      parameter (mx_seg=45) ! 32MB block
c      parameter (mx_seg=50) ! 48MB block
c      parameter (mx_seg=53) ! 64MB block

c INTERNAL VARIABLES
      double precision d1, d2, d3, d4, dCost
      integer mx_occ, mx_vrt, mx_shl, nBasis
      integer iFacts(100), nFacts, iScan
      integer iNeed, iHave, iDiff
      integer i, j, k, nSegs, nOccSegs, nVrtSegs, nAOSegs
      integer iSpin1, iSpin2, pop(2,3)
      character*5 szSpin(2)
      data szSpin/'alpha','beta '/
      character*2 czu(6)
      data czu/'B ','kB','MB','GB','TB','PB'/

      print *, '******************** SEGMENT SIZE ********************'

c ----------------------------------------------------------------------

      mx_occ = max(nocca,noccb)
      mx_vrt = max(nvrta,nvrtb)
      mx_shl = nFpS(1)
      nBasis = nFpS(1)
      do i = 2, nShells
         mx_shl = max(mx_shl,nFpS(i))
         nBasis = nBasis + nFpS(i)
      end do
      print *, 'constraints: ',mx_occ,' occ; ',
     &                         mx_vrt,' vrt; ',
     &                    '<=',mx_shl,' aos/shl (',nBasis,' AOs)'
      if (mx_occ.eq.0.or.mx_vrt.eq.0.or.nBasis.eq.0) then
         print *, '@DET_SEGSIZES: ERROR: 0 dimension encountered'
         call abort_job()
      end if

      pop(1,1) = nocca
      pop(2,1) = noccb
      pop(1,2) = nvrta
      pop(2,2) = nvrtb
      pop(1,3) = nBasis
      pop(2,3) = nBasis

c ----------------------------------------------------------------------

c   o start at mx_vrt and keep going up until a factor is found
      dCost = 1.d100
      iScan = max(mx_vrt,mx_shl)
      if (seg_vrt.eq.0.and.seg_ao.ne.0) seg_vrt=seg_ao
      do while (seg_vrt.eq.0)
         call fact_dim(iScan,iFacts,nFacts)
         print *, 'vrt spectrum: ', !'(pad=',iScan-mx_vrt,')',
     &            (iFacts(i),i=1,nFacts)
         do i = nFacts, 1, -1 ! prefer larger segments to smaller ones
            if (mx_shl.le.iFacts(i).and.iFacts(i).le.mx_seg) then
c            o vrt waste
               nSegs = (mx_vrt+iFacts(i)-1)/iFacts(i)
               iHave = nSegs*iFacts(i)
               d3 = 1.d0*iHave
               d3 = 8.d0*d3**4
               d4 = 1.d0*mx_vrt
               d4 = 8.d0*d4**4
               d1 = d3-d4
               j = 1
               do while (d1.gt.1000.d0)
                  j = j+1
                  d1 = d1/1024.d0
               end do
c            o AO waste
               call binpack_cost(nFpS,nShells,iFacts(i),nAOSegs,d2)
               if (d2.lt.0.99*dCost) then
                  seg_vrt = iFacts(i)
                  dCost = d2
               end if
               k = 1
               do while (d2.gt.1000.d0)
                  k = k+1
                  d2 = d2/1024.d0
               end do
               print '(4x,i2,a,f5.1,x,a,a,f4.1,a,f5.1,x,a,a,i2,a)',
     &                  iFacts(i),': ',d1,czu(j),' wasted (',
     &                  100.d0*(d3-d4)/d3,'%), AO cost = ',d2,czu(k),
     &                  ' (',nAOSegs,' segs)'
            end if
         end do
         iScan = iScan + 1
      end do
      nVrtSegs = (mx_vrt+seg_vrt-1)/seg_vrt

      if (seg_ao.eq.0) seg_ao=seg_vrt
      call binpack_cost(nFpS,nShells,seg_ao,nAOSegs,dCost)

      iNeed = (mx_occ+nVrtSegs-1)/nVrtSegs
      iDiff = mx_occ+1
      iScan = mx_occ
      if (seg_occ.eq.0.and.mx_occ.le.mn_seg) seg_occ = mx_occ
      do while (seg_occ.eq.0.and.iScan.le.mx_occ+3)
         call fact_dim(iScan,iFacts,nFacts)
c         print *, 'occ spectrum: ',(iFacts(i),i=1,nFacts)
         print *, 'occ spectrum: ', !'(pad=',iScan-mx_occ,')',
     &            (iFacts(i),i=1,nFacts)
         do i = 1, nFacts
            if (mn_seg.le.iFacts(i).and.iFacts(i).le.mx_seg) then
               j = abs(iFacts(i)-iNeed)
               if (j.lt.iDiff) then
                  seg_occ = iFacts(i)
                  iDiff = j
               end if
            end if
         end do
         iScan = iScan + 1
      end do
      nOccSegs = (mx_occ+seg_occ-1)/seg_occ

      if (seg_occ.lt.1.or.seg_vrt.lt.1.or.seg_ao.lt.1) then
         print *, 'Error determining segment sizes: (O,V,N) ',
     &            seg_occ, seg_vrt, seg_ao
         call abort_job()
      end if

      print *, 'Summary (leng,numb): (',
     &         seg_occ,',',nOccSegs,') occ, (',
     &         seg_vrt,',',nVrtSegs,') vrt, (',
     &         seg_ao, ',',nAOSegs,') ao'
      print *, '******************************************************'

      return
      end

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------

      subroutine fact_dim(iMult, iList, nFact)
      implicit none

      integer iMult, i, j, k, nFact, iFact1, iFact2, iList(*)
      logical sqr

      i = sqrt(1.*(iMult+0.1))
      sqr = (i*i.eq.iMult)
      nFact = 0
      do iFact1 = i, 1, -1
         iFact2 = iMult/iFact1
         if (iFact1*iFact2.eq.iMult) then
            iList(nFact+1) = iFact1
            iList(nFact+2) = iFact2
            nFact = nFact + 2
            if (iFact1.eq.iFact2) nFact = nFact - 1
         end if
      end do
      k = 3
      if (sqr) k = 2
      do while (k.lt.nFact)
         iFact1 = iList(k)
         do j = k-1, 1, -1
            iList(j+1) = iList(j)
         end do
         iList(1) = iFact1
         k = k + 2
      end do

      return
      end

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------

c This program implements a worst fit bin packing algorithm.

      subroutine binpack_cost(mass,nShells,iSegSize,nSegs,dCost)
      implicit none

      double precision dCost
      integer nShells, mass(*), iSegSize, nSegs

      integer TotMass, lrg, sml

      integer nGrps, pGrp0, SmlGrp

      integer m(2000), GrpMass(2000)
      integer FinGrps

      double precision dTmp, d1, d2
      logical bValid
      logical bAccept

      integer i, j

c ----------------------------------------------------------------------

      if (nShells.lt.1.or.2000.lt.nShells) then
         print *, 'Error in binpack_cost: nshells = ',nshells
         call abort_job()
      end if

      TotMass = 0
      lrg = 0
      sml = 2**30
      do i = 1, nShells
         TotMass = TotMass + mass(i)
         lrg = max(lrg,mass(i))
         sml = min(sml,mass(i))
         m(i) = mass(i)
      end do

      if (sml.lt.1.or.iSegSize.lt.lrg) then
         print *,'Error in binpack_cost: Shells must have a mass > 0,'
         print *,'and the group limit (',iSegSize,
     *           ') must be >= the largest shell size (',lrg,').'
         call abort_job()
      end if

c   o sort the masses in descending order
      call isort(m,0,nShells,-1)

      nGrps = (TotMass+iSegSize-1)/iSegSize
      FinGrps = 0
      do while (FinGrps.eq.0.and.nGrps.le.nShells)

c      o fill in the largest shells
c        (pGrp0 points to the first group to scan from)
         pGrp0 = 1
         do i = 1, nGrps
            GrpMass(i) = m(i)
            if (m(i).gt.iSegSize-sml) pGrp0 = pGrp0 + 1
         end do
         bValid=(pGrp0.le.nGrps.or.nGrps.eq.nShells)

c      o assign the remaining shells
         do i = nGrps+1, nShells
c         o find the smallest group
            SmlGrp = pGrp0
            do j = pGrp0+1, nGrps
               if (GrpMass(j).lt.GrpMass(SmlGrp)) SmlGrp = j
            end do
            GrpMass(SmlGrp) = GrpMass(SmlGrp) + m(i)
            if (GrpMass(SmlGrp).gt.iSegSize) bValid = .false.
         end do

         if (bValid) FinGrps = nGrps
         nGrps = nGrps + 1
      end do

      nSegs = FinGrps

      dCost = 0.d0
      if (TotMass.ne.FinGrps*iSegSize) then
         d1 = 1.d0*TotMass
         d2 = 1.d0*FinGrps*iSegSize
         d1 = 8.d0*d1**4
         d2 = 8.d0*d2**4
         dCost = (d2-d1) !/1024/1024/1024
      end if

      return
      end

