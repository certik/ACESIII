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
      subroutine compute_der2_batch(a1,a2,b1,b2,c1,c2,d1,d2,
     *                 scr,maxblk,
     *                 iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 in, nsend)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'parallel_info.h'
      include 'hess.h'

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer der_flags(12)
      integer der_save(12)
      integer der_save1(12)
      integer der_test(12)
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, k, l, n, m, r, s
      integer a,b,c,d
      integer wder, ncder, nder_first  

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst

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
      double precision x4,y4,z4

      double precision x1p,y1p,z1p
      double precision x2p,y2p,z2p
      double precision x3p,y3p,z3p
      double precision x4p,y4p,z4p

      integer ihess, jhess, icomponent, jcomponent 
      integer iatom, jatom, iflag, jflag   
      integer katom, kx, latom, lx  

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision in(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision scr(*), eps, fact, y
      integer iscr(*)
      integer myder(4,3), mySder(4,3) 

      integer ccbeg(*), ccend(*)
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)

      common /d4int_com/katom, kx, latom, lx  
      save alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1 
      l8true = .false. ! screening .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do d = d1,d2
      do c = c1,c2
      do b = b1,b2
      do a = a1,a2
         in(a,b,c,d) = 0.d0
      enddo
      enddo
      enddo
      enddo

      nsend = adim*bdim*cdim*ddim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim,ddim = ',adim,bdim,cdim,ddim
         call abort_job()
      endif

c----------------------------------------------------------------------------
c   Set the flags array for the FIRST derivative.
c----------------------------------------------------------------------------

      do iflag = 1, 12 

         do jflag = 1, 12 
            der_flags(jflag) = 0 
            der_save(jflag)  = 0 
         enddo 

         der_flags(iflag) = 1 
         der_save(iflag)  = 1 

c----------------------------------------------------------------------------
c   Save the der_array.
c----------------------------------------------------------------------------

         do a = 1, 12 
            der_save(a) = der_flags(a) 
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

c----------------------------------------------------------------------------
c   Make sure the der_flags are set correctly for equivalent centers.
c----------------------------------------------------------------------------

               do a = 1, 12 
                  der_flags(a) = der_save(a) 
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

c----------------------------------------------------------------------------
c   Set the Hessian component for the first derivative.  
c----------------------------------------------------------------------------

               iatom = 0
               icomponent = 0

               do a = 1, 3
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(m)
                     icomponent = a
c                    if ((iatom.eq.katom).and.(icomponent.eq.kx)) then 
c                       jatom      = latom 
c                       jcomponent = lx 
c                       go to 33
c                    endif 
                     if ((iatom.eq.latom).and.(icomponent.eq.lx)) then 
                        jatom      = katom 
                        jcomponent = kx 
                        go to 33
                     endif 
                  endif
               enddo

               do a = 4, 6
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(n)
                     icomponent = a-3
c                    if ((iatom.eq.katom).and.(icomponent.eq.kx)) then 
c                       jatom      = latom 
c                       jcomponent = lx 
c                       go to 33
c                    endif 
                     if ((iatom.eq.latom).and.(icomponent.eq.lx)) then 
                        jatom      = katom 
                        jcomponent = kx 
                        go to 33
                     endif 
                  endif
               enddo 

               do a = 7, 9
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(r)
                     icomponent = a-6
c                    if ((iatom.eq.katom).and.(icomponent.eq.kx)) then 
c                       jatom      = latom 
c                       jcomponent = lx 
c                       go to 33
c                    endif 
                     if ((iatom.eq.latom).and.(icomponent.eq.lx)) then 
                        jatom      = katom 
                        jcomponent = kx 
                        go to 33
                     endif 
                  endif
               enddo 

               do a = 10, 12
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(s)
                     icomponent = a-9
c                    if ((iatom.eq.katom).and.(icomponent.eq.kx)) then 
c                       jatom      = latom 
c                       jcomponent = lx 
c                       go to 33
c                    endif 
                     if ((iatom.eq.latom).and.(icomponent.eq.lx)) then 
                        jatom      = katom 
                        jcomponent = kx 
                        go to 33
                     endif 
                  endif
               enddo 
               go to 666 

33             continue 

c--------------------------------------------------------------------------
c   Set the perturbed arguments.
c--------------------------------------------------------------------------


               if ((jatom .ne. atom(m)) .and. (jatom .ne. atom(n))
     *         .and. (jatom .ne. atom(r)) .and. (jatom .ne. atom(s)))
     *         go to 555

               do jflag = 1, 12
                  der_test(jflag)  = der_save(jflag)
               enddo 

c-------------------------------------------------------------------------
c
c XX integrals first.
c
c-------------------------------------------------------------------------
c
             if (jcomponent .eq. 1) then

                if (atom(m) .eq. jatom) then
                   der_test(1) = der_test(1) + 1
                endif

                if (atom(n) .eq. jatom) then
                   der_test(4) = der_test(4) + 1
                endif

                if (atom(r) .eq. jatom) then
                   der_test(7) = der_test(7) + 1
                endif

                if (atom(s) .eq. jatom) then
                   der_test(10) = der_test(10) + 1 
                endif

                go to 111

             endif
c
c-------------------------------------------------------------------------
c
c YY integrals first.
c
c-------------------------------------------------------------------------
c
             if (jcomponent .eq. 2) then

                if (atom(m) .eq. jatom) then
                   der_test(2) = der_test(2) + 1
                endif

                if (atom(n) .eq. jatom) then
                   der_test(5) = der_test(5) + 1
                endif

                if (atom(r) .eq. jatom) then
                   der_test(8) = der_test(8) + 1
                endif

                if (atom(s) .eq. jatom) then
                   der_test(11) = der_test(11) + 1
                endif

                go to 111

             endif
c
c-------------------------------------------------------------------------
c
c ZZ integrals next.
c
c-------------------------------------------------------------------------
c
             if (jcomponent .eq. 3) then

                if (atom(m) .eq. jatom) then
                   der_test(3) = der_test(3) + 1
                endif

                if (atom(n) .eq. jatom) then
                   der_test(6) = der_test(6) + 1
                endif

                if (atom(r) .eq. jatom) then
                   der_test(9) = der_test(9) + 1 
                endif

                if (atom(s) .eq. jatom) then
                   der_test(12) = der_test(12) + 1
                endif

                go to 111

             endif

c-------------------------------------------------------------------------
c
c If you get to this point then no second-derivative can be taken so 
c skip computation. 
c
             go to 777 
c
111          continue

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

c----------------------------------------------------------------------------
c   Make sure the der_flags are set correctly for equivalent centers.
c----------------------------------------------------------------------------

               if (atom(m) .eq. atom(n)) then
                  if (der_test(1) .ne. 0) der_test(4) = der_test(1)
                  if (der_test(4) .ne. 0) der_test(1) = der_test(4)
                  if (der_test(2) .ne. 0) der_test(5) = der_test(2)
                  if (der_test(5) .ne. 0) der_test(2) = der_test(5)
                  if (der_test(3) .ne. 0) der_test(6) = der_test(3)
                  if (der_test(6) .ne. 0) der_test(3) = der_test(6)
               endif

               if (atom(m) .eq. atom(r)) then
                  if (der_test(1) .ne. 0) der_test(7) = der_test(1)
                  if (der_test(7) .ne. 0) der_test(1) = der_test(7)
                  if (der_test(2) .ne. 0) der_test(8) = der_test(2)
                  if (der_test(8) .ne. 0) der_test(2) = der_test(8)
                  if (der_test(3) .ne. 0) der_test(9) = der_test(3)
                  if (der_test(9) .ne. 0) der_test(3) = der_test(9)
               endif

               if (atom(m) .eq. atom(s)) then
                  if (der_test(1) .ne. 0) der_test(10) = der_test(1)
                  if (der_test(10).ne. 0) der_test(1)  = der_test(10)
                  if (der_test(2) .ne. 0) der_test(11) = der_test(2)
                  if (der_test(11).ne. 0) der_test(2)  = der_test(11)
                  if (der_test(3 ).ne. 0) der_test(12) = der_test(3)
                  if (der_test(12).ne. 0) der_test(3)  = der_test(12)
               endif

               if (atom(n) .eq. atom(r)) then
                  if (der_test(4) .ne. 0) der_test(7) = der_test(4)
                  if (der_test(7) .ne. 0) der_test(4) = der_test(7)
                  if (der_test(5) .ne. 0) der_test(8) = der_test(5)
                  if (der_test(8) .ne. 0) der_test(5) = der_test(8)
                  if (der_test(6) .ne. 0) der_test(9) = der_test(6)
                  if (der_test(9) .ne. 0) der_test(6) = der_test(9)
               endif

               if (atom(n) .eq. atom(s)) then
                  if (der_test(4)  .ne. 0) der_test(10) = der_test(4)
                  if (der_test(10) .ne. 0) der_test(4)  = der_test(10)
                  if (der_test(5)  .ne. 0) der_test(11) = der_test(5)
                  if (der_test(11) .ne. 0) der_test(5)  = der_test(11)
                  if (der_test(6)  .ne. 0) der_test(12) = der_test(6)
                  if (der_test(12) .ne. 0) der_test(6)  = der_test(12)
               endif

               if (atom(r) .eq. atom(s)) then
                  if (der_test(7)  .ne. 0) der_test(10) = der_test(7)
                  if (der_test(10) .ne. 0) der_test(7)  = der_test(10)
                  if (der_test(8)  .ne. 0) der_test(11) = der_test(8)
                  if (der_test(11) .ne. 0) der_test(8)  = der_test(11)
                  if (der_test(9)  .ne. 0) der_test(12) = der_test(9)
                  if (der_test(12) .ne. 0) der_test(9)  = der_test(12)
               endif

               ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)

                  call ERD__GENER_ERI_DERV_BATCH(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum,
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n),
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4,
     *                der_test(1), der_test(2), der_test(3),
     *                der_test(4), der_test(5), der_test(6),
     *                der_test(7), der_test(8), der_test(9),
     *                der_test(10), der_test(11), der_test(12),
     *                alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints,
     *                nfirst, scr)

               if (nints .gt. 0) then

               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

c---------------------------------------------------------------------------
c   Move the integrals to the output block.
c---------------------------------------------------------------------------

               call add_integrals(in, a1,a2,b1,b2,c1,c2,d1,d2,
     *              scr(nfirst),aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,1.0d0)

c---------------------------------------------------------------------------

            endif

777         continue 

555         continue

666         continue 
            enddo  ! iflag = 1, 12 

         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m

      return
      end

      
