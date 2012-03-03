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
      subroutine compute_der_batch(a1,a2,b1,b2,c1,c2,d1,d2,scr,maxblk,
     *                 iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 out, nsend)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer a,b,c,d

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
      logical spherical

      double precision x1,y1,z1
      double precision x2,y2,z2
      double precision x3,y3,z3
      double precision x4,y4,z4

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision out(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision scr(*)   
      integer iscr(*)

      integer ccbeg(*), ccend(*)
      integer der_flags(12), der_save(12)   
      integer myder(4,3), mySder(4,3) 

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)

      common /d2int_com/jatom, jx, kcenter
      integer jatom, jx, jcenter, dcoord, wder, k, ncder  
      integer katom, kx, kcenter  

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      call mpi_comm_rank(mpi_comm_world, me, ierr)
c      print *,'Task ',me,' computing integrals for ',a1,a2,b1,b2,
c     *     c1,c2,d1,d2
c      call c_flush_stdout()

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1 
      spherical = (ispherical .eq. 1)
  
      nsend = adim*bdim*cdim*ddim
      if (nsend .le. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim,ddim = ',adim,bdim,cdim,ddim
         call abort_job()
      endif

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do d = d1,d2
      do c = c1,c2
      do b = b1,b2
      do a = a1,a2
         out(a,b,c,d) = 0.d0
      enddo
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
         call lookup_shell(end_nfps, nshells, d1, s1)
         call lookup_shell(end_nfps, nshells, d2, s2)

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
         do s = s1, s2
         
               x4 = coords(1,s)
               y4 = coords(2,s)
               z4 = coords(3,s)
               call pack_coeffs(alphas, ixalpha, coeffs, ixpcoef,
     *                          ncfps, npfps, m, n,
     *                          r, s, alpha_pack, nalpha_pack,
     *                          pcoeff_pack, npcoeff_pack,
     *                          ccbeg, ccend, indx_cc,
     *                          ccbeg_pack, ccend_pack)
c
c Loop over the 4-centers and set derivative flags.  
c ------------------------------------------------- 

c     do jcenter = 1, 4 

      do a = 1, 12
         der_flags(a) = 0 
      enddo

c     write(6,*) ' Computing derivative integral batch defined by ' 
c     write(6,*) ' jatom = ', jatom 
c     write(6,*) ' jx    = ', jx 
c     write(6,*) ' jcenter = ', jcenter 
c     write(6,*) ' flags = ', dcoord  
c     write(6,*) ' A B C D = ', a1,a2,':',b1,b2,':',c1,c2,':',d1,d2   
c     write(6,*) ' NSEND = ', nsend 

c----------------------------------------------------------------------------
c   Check if the derivative center matches the atom. 
c----------------------------------------------------------------------------

      if ((atom(m) .ne. jatom) .and. (atom(n) .ne. jatom) .and.
     *    (atom(r) .ne. jatom) .and. (atom(s) .ne. jatom)) go to 777 

      if (atom(m) .eq. jatom) der_flags(jx)   = 1  
      if (atom(n) .eq. jatom) der_flags(jx+3) = 1  
      if (atom(r) .eq. jatom) der_flags(jx+6) = 1  
      if (atom(s) .eq. jatom) der_flags(jx+9) = 1  

c----------------------------------------------------------------------------
c   Save the der_array.
c----------------------------------------------------------------------------

      do a = 1, 12
         der_save(a) = der_flags(a)
      enddo

               if (atom(m) .eq. atom(n)) then
                  if (der_flags(1) .eq. 1) der_flags(4) = 1
                  if (der_flags(4) .eq. 1) der_flags(1) = 1
                  if (der_flags(2) .eq. 1) der_flags(5) = 1
                  if (der_flags(5) .eq. 1) der_flags(2) = 1
                  if (der_flags(3) .eq. 1) der_flags(6) = 1
                  if (der_flags(6) .eq. 1) der_flags(3) = 1
               endif

               if (atom(m) .eq. atom(r)) then
                  if (der_flags(1) .eq. 1) der_flags(7) = 1
                  if (der_flags(7) .eq. 1) der_flags(1) = 1
                  if (der_flags(2) .eq. 1) der_flags(8) = 1
                  if (der_flags(8) .eq. 1) der_flags(2) = 1
                  if (der_flags(3) .eq. 1) der_flags(9) = 1
                  if (der_flags(9) .eq. 1) der_flags(3) = 1
               endif

               if (atom(m) .eq. atom(s)) then
                  if (der_flags(1) .eq. 1) der_flags(10) = 1
                  if (der_flags(10).eq. 1) der_flags(1)  = 1
                  if (der_flags(2) .eq. 1) der_flags(11) = 1
                  if (der_flags(11).eq. 1) der_flags(2)  = 1
                  if (der_flags(3 ).eq. 1) der_flags(12) = 1
                  if (der_flags(12).eq. 1) der_flags(3)  = 1
               endif

               if (atom(n) .eq. atom(r)) then
                  if (der_flags(4) .eq. 1) der_flags(7) = 1
                  if (der_flags(7) .eq. 1) der_flags(4) = 1
                  if (der_flags(5) .eq. 1) der_flags(8) = 1
                  if (der_flags(8) .eq. 1) der_flags(5) = 1
                  if (der_flags(6) .eq. 1) der_flags(9) = 1
                  if (der_flags(9) .eq. 1) der_flags(6) = 1
               endif

               if (atom(n) .eq. atom(s)) then
                  if (der_flags(4) .eq. 1) der_flags(10) = 1
                  if (der_flags(10) .eq. 1) der_flags(4) = 1
                  if (der_flags(5) .eq. 1) der_flags(11) = 1
                  if (der_flags(11) .eq. 1) der_flags(5) = 1
                  if (der_flags(6) .eq. 1) der_flags(12) = 1
                  if (der_flags(12) .eq. 1) der_flags(6) = 1
               endif

               if (atom(r) .eq. atom(s)) then
                  if (der_flags(7) .eq. 1) der_flags(10) = 1
                  if (der_flags(10) .eq. 1) der_flags(7) = 1
                  if (der_flags(8) .eq. 1) der_flags(11) = 1
                  if (der_flags(11) .eq. 1) der_flags(8) = 1
                  if (der_flags(9) .eq. 1) der_flags(12) = 1
                  if (der_flags(12) .eq. 1) der_flags(9) = 1
               endif

c---------------------------------------------------------------------------
c   Check if you want to calculate the integral batch.
c---------------------------------------------------------------------------

               ncder   = 0
               k       = 0
               wder    = 0
               nfirst  = 0
               do i = 1, 4
                  do j = 1, 3
                     k = k + 1
                     myder(i,j)  = der_flags(k)
                     mySder(i,j) = der_save(k)
                     if (der_flags(k) .eq. 1) wder  = j
                     if (der_flags(k) .eq. 1) ncder = ncder + 1
                  enddo
               enddo

               if (ncder .eq. 2) then
                  do i = 1, 4
                     if (myder(i,wder) .ne. 0) then
                        nfirst = i
                        go to 10
                     endif
                  enddo
10                continue
                  if (mySder(nfirst,wder) .ne. 1) go to 777
               endif ! ncder .eq. 2

               if (ncder .eq. 4) then
                  do i = 1, 4
                     if (myder(i,wder) .ne. 0) then
                        nfirst = i
                        go to 12
                     endif
                  enddo
12                continue
                  if (mySder(nfirst,wder) .ne. 1) go to 777
               endif ! ncder .eq. 4

               if (ncder .eq. 3) then
                  do i = 1, 4
                     if (myder(i,wder) .ne. 0) then
                        nfirst = i
                        go to 11
                     endif
                  enddo
11                continue
                  if (mySder(nfirst,wder) .ne. 1) go to 777
               endif ! ncder .eq. 3

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------

               ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)

                  call ERD__GENER_ERI_DERV_BATCH(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     *                der_flags(1), der_flags(2), der_flags(3),
     *                der_flags(4), der_flags(5), der_flags(6),
     *                der_flags(7), der_flags(8), der_flags(9),
     *                der_flags(10), der_flags(11), der_flags(12),
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints,
     *                nfirst, scr)

c---------------------------------------------------------------------------
c   Move the integrals into the output block. 
c---------------------------------------------------------------------------

            if (nints .gt. 0) then

               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

               call add_integrals(out, a1,a2,b1,b2,c1,c2,d1,d2,
     *              scr(nfirst),aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)

            endif

777         continue

c        enddo ! jcenter = 1, 4 

         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m
 
      return
      end

