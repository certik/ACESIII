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
      subroutine tester_MRT2(a1,a2,b1,b2,c1,c2,d1,d2,
     *                 der_flags, scr,maxblk,
     *                 iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 in, hess)
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

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision out(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision temp(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision in(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision scr(*), eps, fact, y
      double precision hess(3*ncenters,3*ncenters)
      integer iscr(*)
      integer myder(4,3), mySder(4,3) 

      integer ccbeg(*), ccend(*)
      integer ccbeg_pack(2000), ccend_pack(2000)
      double precision alpha_pack(2000), pcoeff_pack(2000)
      save alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1 
      l8true = .false. ! screening .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical
  
      nsend = adim*bdim*cdim*ddim
      if (nsend .lt. 0) then
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

c              write(6,*) 'SATOM :', atom(m),atom(n),atom(r),atom(s)

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

c--------------------------------------------------------------------------
c   Set the hessian components of the first-derivative.
c--------------------------------------------------------------------------

               iatom = 0
               icomponent = 0

               do a = 1, 3
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(m)
                     icomponent = a
                     go to 33
                  endif
               enddo

               do a = 4, 6
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(n)
                     icomponent = a-3
                     go to 33
                  endif
               enddo 

               do a = 7, 9
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(r)
                     icomponent = a-6
                     go to 33
                  endif
               enddo 

               do a = 10, 12
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(s)
                     icomponent = a-9
                     go to 33
                  endif
               enddo 

33             continue 

               eps = 0.0000001d0 

c--------------------------------------------------------------------------
c   Set the perturbed arguments.
c--------------------------------------------------------------------------

               do jatom = 1, ncenters 

                  if ((jatom .ne. atom(m)) .and. (jatom .ne. atom(n))
     *            .and. (jatom .ne. atom(r)) .and. (jatom .ne. atom(s)))
     *            go to 555

               do jcomponent = 1, 3 

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

                x1p = x1
                y1p = y1
                z1p = z1

                x2p = x2
                y2p = y2
                z2p = z2

                x3p = x3
                y3p = y3
                z3p = z3

                x4p = x4
                y4p = y4
                z4p = z4

                if (atom(m) .eq. jatom) then

                   x1p = x1 + eps
                   der_test(1) = der_test(1) + 1

                endif

                if (atom(n) .eq. jatom) then

                   x2p = x2 + eps
                   der_test(4) = der_test(4) + 1

                endif

                if (atom(r) .eq. jatom) then

                   x3p = x3 + eps
                   der_test(7) = der_test(7) + 1

                endif

                if (atom(s) .eq. jatom) then

                   x4p = x4 + eps
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

                x1p = x1
                y1p = y1
                z1p = z1

                x2p = x2
                y2p = y2
                z2p = z2

                x3p = x3
                y3p = y3
                z3p = z3

                x4p = x4
                y4p = y4
                z4p = z4

                if (atom(m) .eq. jatom) then

                   y1p = y1 + eps
                   der_test(2) = der_test(2) + 1

                endif

                if (atom(n) .eq. jatom) then

                   y2p = y2 + eps
                   der_test(5) = der_test(5) + 1

                endif

                if (atom(r) .eq. jatom) then

                   y3p = y3 + eps
                   der_test(8) = der_test(8) + 1

                endif

                if (atom(s) .eq. jatom) then

                   y4p = y4 + eps
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

                x1p = x1
                y1p = y1
                z1p = z1

                x2p = x2
                y2p = y2
                z2p = z2

                x3p = x3
                y3p = y3
                z3p = z3

                x4p = x4
                y4p = y4
                z4p = z4

                if (atom(m) .eq. jatom) then

                   z1p = z1 + eps
                   der_test(3) = der_test(3) + 1

                endif

                if (atom(n) .eq. jatom) then

                   z2p = z2 + eps
                   der_test(6) = der_test(6) + 1

                endif

                if (atom(r) .eq. jatom) then

                   z3p = z3 + eps
                   der_test(9) = der_test(9) + 1 

                endif

                if (atom(s) .eq. jatom) then

                   z4p = z4 + eps
                   der_test(12) = der_test(12) + 1

                endif

                go to 111

             endif

                  do jflag = 1, 12
                     der_save1(jflag) = der_save(jflag)
                  enddo 
c
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
               
c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2 stationary integrals.
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
 
               call move_integrals(out, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2)
             else 
                 go to 777 
            endif
               
c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2 perturbed integrals.
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
c   Add the integrals into the output block.  
c---------------------------------------------------------------------------

            if (nints .gt. 0) then
               
               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)
 
               fact = -1.0d0 
               call add_integrals(out, a1,a2,b1,b2,c1,c2,d1,d2,
     *                            scr(nfirst), 
     *                            aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,fact)
c              write(6,*) 'SATOM :', atom(m),atom(n),atom(r),atom(s)
c              write(6,*) ' ATOM :', iatom,jatom
c              write(6,*) ' COMP :', icomponent, jcomponent
c              write(6,*) ' FLAGS:', (der_test(i),i=1,12)

c---------------------------------------------------------------------------
c   Contract with Intermediate block.
c---------------------------------------------------------------------------

c              y = 0.0d0
c              do l = dd1, dd2
c              do k = cc1, cc2
c              do j = bb1, bb2
c              do i = aa1, aa2
c                 write(6,*) i,j,k,l,out(i,j,k,l)/eps 
c                 y = y + in(i,j,k,l)*out(i,j,k,l)
c              enddo
c              enddo
c              enddo
c              enddo 

c---------------------------------------------------------------------------
c   Sum into the hessian.
c---------------------------------------------------------------------------

c              ihess = (iatom-1)*3 + icomponent
c              jhess = (jatom-1)*3 + jcomponent

c              hess(jhess,ihess) = hess(jhess,ihess) + y/eps 

            endif

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do d = d1,d2
      do c = c1,c2
      do b = b1,b2
      do a = a1,a2
         temp(a,b,c,d) = 0.d0
      enddo
      enddo
      enddo
      enddo

c           go to 777 

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

c -------------------------------------------------------------------------
c   Set the derivative flag arguments.
c--------------------------------------------------------------------------

c              write(6,*) '       ATOM :', iatom,jatom
c              write(6,*) '       COMP :', icomponent, jcomponent
c              write(6,*) '       FLAGS:', (der_test(a),a=1,12)

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

               if (nints .gt. 0) then

               if (s .eq. 1) then
                  dd1 = 1
               else
                  dd1 = end_nfps(s-1) + 1
               endif
               dd2 = end_nfps(s)

               call move_integrals(temp, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst),
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2)

c              go to 987 

               do l = dd1, dd2
               do k = cc1, cc2
               do j = bb1, bb2
               do i = aa1, aa2

               if (dabs(temp(i,j,k,l)-out(i,j,k,l)/eps) .gt. 
     *            0.000001d0) then 

               write(6,*) ' MISMATCH IN SECOND DERIVATIVE' 
               write(6,*) 'SATOM :', atom(m),atom(n),atom(r),atom(s)
               write(6,*) ' ATOM :', iatom,jatom
               write(6,*) ' COMP :', icomponent, jcomponent
               write(6,*) ' FLAGS:', (der_test(a),a=1,12)
               write(6,*) i, j, k, l, out(i,j,k,l)/eps, temp(i,j,k,l)

               endif 

               enddo
               enddo
               enddo
               enddo

987            continue 

c---------------------------------------------------------------------------
c   Contract with Intermediate block.
c---------------------------------------------------------------------------

               y = 0.0d0
               do l = dd1, dd2
               do k = cc1, cc2
               do j = bb1, bb2
               do i = aa1, aa2
                  y = y + in(i,j,k,l)*temp(i,j,k,l)
               enddo
               enddo
               enddo
               enddo

c---------------------------------------------------------------------------
c   Sum into the hessian.
c---------------------------------------------------------------------------

               ihess = (iatom-1)*3 + icomponent
               jhess = (jatom-1)*3 + jcomponent

               hess(jhess,ihess) = hess(jhess,ihess) + y

            endif

777         continue 

            enddo ! jatom = 1, ncenters 
555         continue
            enddo ! jcomponent = 1, 3 

            enddo  ! iflag = 1, 12 

         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m

      return
      end

      
