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
      subroutine compute_erd3C(a1,a2,b1,b2,c1,c2,scr,maxblk,
     *                        iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                        out, nsend)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1, a2, b1, b2, c1, c2  
      integer aa1,aa2,bb1,bb2,cc1,cc2
      integer adim, bdim, cdim   
      integer m1, m2, n1, n2, r1, r2 
      integer i, j, n, m, r 
      integer a,b,c

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst
      integer me, ierr

      integer imin, zmin, iblk, zblk

      logical skip
      logical mn_symmetry
      logical rs_symmetry
      logical mn_rs_symmetry
      logical*8 l8true, l8spherical
      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision x3,y3,z3

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision out(a1:a2,b1:b2,c1:c2)
      double precision scr(*)   
      integer iscr(*)

      integer ccbeg(*), ccend(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      call mpi_comm_rank(mpi_comm_world, me, ierr)

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1

      l8true      = .true.
      spherical   = (ispherical .eq. 1)
      l8spherical = spherical
  
      nsend = adim*bdim*cdim
      if (nsend .lt. 0) then
         print *,'ERROR IN ERD3INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim = ',adim,bdim,cdim
         call abort_job()
      endif

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do c = c1,c2
      do b = b1,b2
      do a = a1,a2
         out(a,b,c) = 0.d0
      enddo
      enddo
      enddo

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

         call lookup_shell(end_nfps, nshells, a1, m1)
         call lookup_shell(end_nfps, nshells, a2, m2)
         call lookup_shell(end_nfps, nshells, b1, n1)
         call lookup_shell(end_nfps, nshells, b2, n2)
         call lookup_shell(end_nfps, nshells, c1, r1)
         call lookup_shell(end_nfps, nshells, c2, r2)

         do m = m1, m2
            if (m .eq. 1) then
               aa1 = 1
            else
               aa1 = end_nfps(m-1) + 1
            endif
            aa2 = end_nfps(m)

            x1 = coords(1,m)
            y1 = coords(2,m)
            z1 = coords(3,m)
         do n = n1, n2
            if (n .eq. 1) then
               bb1 = 1
            else
               bb1 = end_nfps(n-1) + 1
            endif
            bb2 = end_nfps(n)

            x2 = coords(1,n)
            y2 = coords(2,n)
            z2 = coords(3,n)
         do r = r1, r2
            if (r .eq. 1) then
               cc1 = 1
            else
               cc1 = end_nfps(r-1) + 1
            endif
            cc2 = end_nfps(r)

            x3 = coords(1,r)
            y3 = coords(2,r)
            z3 = coords(3,r)

            call pack_coeffs_ovl3c(alphas, ixalpha, coeffs, ixpcoef, 
     *                       ncfps, npfps, m, n, r,  
     *                       alpha_pack, nalpha_pack, 
     *                       pcoeff_pack, npcoeff_pack, 
     *                       ccbeg, ccend, indx_cc,
     *                       ccbeg_pack, ccend_pack)

c---------------------------------------------------------------------------
c   Calling sequence for ERD_3C.
c---------------------------------------------------------------------------

               ncsum = ncfps(m) + ncfps(n) + ncfps(r) 

c                 call OED__GENER_OVL3C_BATCH(intmax, zmax,
c    *                nalpha_pack, npcoeff_pack, ncsum, 
c    *                ncfps(m),ncfps(n), ncfps(r), 
c    *                npfps(m),npfps(n), npfps(r), 
c    *                ivangmom(m), ivangmom(n), ivangmom(r), 
c    *                atom(m), atom(n), atom(r), x1,y1,z1,
c    *                x2,y2,z2,x3,y3,z3, alpha_pack,
c    *                pcoeff_pack, ccbeg_pack, ccend_pack,
c    *                spherical, .true., iscr, nints, 
c    *                nfirst, scr)    

c---------------------------------------------------------------------------
c   Move the integrals into the output block. 
c---------------------------------------------------------------------------

            if (nints .gt. 0) then
               
               call move_OVL3C(out,a1,a2,b1,b2,c1,c2, 
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2) 
            endif

         enddo   ! r

         enddo   ! n
         enddo   ! m

      return
      end

      
