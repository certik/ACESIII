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
      subroutine compute_erd_2e_Sderv_integrals(a1,a2,b1,b2,c1,c2,d1,d2,
     *                 der_flags,
     *                 scr,maxblk,
     *                 iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 in, hess)
c---------------------------------------------------------------------------
c   "Work" routine for the integral worker task.
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'hess.h'

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer der_flags(12)
      integer der_save(12)
      integer iatom, jatom, icomponent, jcomponent, ihess, jhess   
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, k, l, n, m, r, s
      integer a,b,c,d
      integer ix, jx 
      integer wder, ncder, nder_first  
      integer nflags,iflags(256,12) 

      integer num_to_do, nsend
      integer nints, maxblk
      integer nalpha_pack, npcoeff_pack
      integer ncsum, next, nfirst, ifirst, jfirst, kfirst 
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
      double precision x4,y4,z4

      double precision coords(3,*), coeffs(*), alphas(*)
      double precision out(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision in(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision scr(*), y, fact
      double precision hess(3*ncenters,3*ncenters)
      integer iscr(*)
      integer myder(4,3), mySder(4,3), satom(4)  
      integer ncder1, ncder2, wder1, wder2, nfirst1, nfirst2   

      integer ccbeg(*), ccend(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack

      call mpi_comm_rank(mpi_comm_world, me, ierr)
c      print *,'Task ',me,' computing integrals for ',a1,a2,b1,b2,
c     *     c1,c2,d1,d2

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1 
      l8true = .true.
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
            if (intpkg .eq. flocke_package) then
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
c   Check if the integral is atomic. 
c----------------------------------------------------------------------------

           if ((atom(m) .eq. atom(n)) .and. 
     *         (atom(m) .eq. atom(r)) .and.    
     *         (atom(m) .eq. atom(s))) go to 777   

c----------------------------------------------------------------------------
c   Determine the flags for this quad. 
c----------------------------------------------------------------------------

            satom(1) = atom(m) 
            satom(2) = atom(n) 
            satom(3) = atom(r) 
            satom(4) = atom(s) 

            call setflag(m,n,r,s,satom,nflags,iflags) 
            if (nflags .eq. 0) go to 777 

            do ix = 1, nflags     

               do a = 1, 12 
                  der_save(a)  = iflags(ix,a)   
                  der_flags(a) = iflags(ix,a)   
               enddo  

c----------------------------------------------------------------------------
c   Make sure the der_flags are set correctly for equivalent centers.
c----------------------------------------------------------------------------

               iatom = 0 
               jatom = 0 
               icomponent = 0 
               jcomponent = 0 

               do a = 1, 3  
                  if (der_flags(a) .eq. 2) then  
                     iatom = atom(m)  
                     jatom = atom(m)  
                     icomponent = a 
                     jcomponent = a 
                     go to 33 
                  endif 
               enddo  

               do a = 4, 6  
                  if (der_flags(a) .eq. 2) then  
                     iatom = atom(n)  
                     jatom = atom(n)  
                     icomponent = a-3  
                     jcomponent = a-3  
                     go to 33 
                  endif 
               enddo  

               do a = 7, 9  
                  if (der_flags(a) .eq. 2) then  
                     iatom = atom(r)  
                     jatom = atom(r)  
                     icomponent = a-6  
                     jcomponent = a-6  
                     go to 33 
                  endif 
               enddo  

               do a = 10, 12  
                  if (der_flags(a) .eq. 2) then  
                     iatom = atom(s)  
                     jatom = atom(s)  
                     icomponent = a-9  
                     jcomponent = a-9  
                     go to 33 
                  endif 
               enddo  

               do a = 1, 3
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(m)
                     icomponent = a
                     kfirst = a  
                     go to 66
                  endif
               enddo

               do a = 4, 6
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(n)
                     icomponent = a-3 
                     kfirst = a  
                     go to 66
                  endif
               enddo

               do a = 7, 9  
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(r)
                     icomponent = a-6  
                     kfirst = a  
                     go to 66
                  endif
               enddo

               do a = 10, 12  
                  if (der_flags(a) .eq. 1) then
                     iatom = atom(s)
                     icomponent = a-9  
                     kfirst = a  
                     go to 66
                  endif
               enddo

66             continue

               do a = 1, 3
                  if ((der_flags(a) .eq. 1).and.(a .gt. kfirst)) then
                     jatom = atom(m)
                     jcomponent = a
                     go to 33
                  endif
               enddo

               do a = 4, 6
                  if ((der_flags(a) .eq. 1).and.(a .gt. kfirst)) then
                     jatom = atom(n)
                     jcomponent = a-3 
                     go to 33
                  endif
               enddo

               do a = 7, 9  
                  if ((der_flags(a) .eq. 1).and.(a .gt. kfirst)) then
                     jatom = atom(r)
                     jcomponent = a-6  
                     go to 33
                  endif
               enddo

               do a = 10, 12  
                  if ((der_flags(a) .eq. 1).and.(a .gt. kfirst)) then
                     jatom = atom(s)
                     jcomponent = a-9  
                     go to 33
                  endif
               enddo

33             continue 

               kfirst = 0
               do k = 1, 12
                  if(der_flags(k) .eq. 2) kfirst = 2
                  if(der_flags(k) .eq. 1) kfirst = 1
               enddo

               if (kfirst .eq. 0) then 
                  write(6,*) ' NO DER_FLAGS SET! ' 
                  stop 
               endif  

               if (atom(m) .eq. atom(n)) then 
                 if (der_flags(1) .ne. 0) der_flags(4) = der_flags(1)  
                 if (der_flags(4) .ne. 0) der_flags(1) = der_flags(4)  
                 if (der_flags(2) .ne. 0) der_flags(5) = der_flags(2)  
                 if (der_flags(5) .ne. 0) der_flags(2) = der_flags(5) 
                 if (der_flags(3) .ne. 0) der_flags(6) = der_flags(3) 
                 if (der_flags(6) .ne. 0) der_flags(3) = der_flags(6) 
               endif 

               if (atom(m) .eq. atom(r)) then 
                 if (der_flags(1) .ne. 0) der_flags(7) = der_flags(1) 
                 if (der_flags(7) .ne. 0) der_flags(1) = der_flags(7) 
                 if (der_flags(2) .ne. 0) der_flags(8) = der_flags(2)  
                 if (der_flags(8) .ne. 0) der_flags(2) = der_flags(8) 
                 if (der_flags(3) .ne. 0) der_flags(9) = der_flags(3) 
                 if (der_flags(9) .ne. 0) der_flags(3) = der_flags(9) 
               endif 

               if (atom(m) .eq. atom(s)) then 
                 if (der_flags(1) .ne. 0) der_flags(10) = der_flags(1) 
                 if (der_flags(10).ne. 0) der_flags(1)  = der_flags(10) 
                 if (der_flags(2) .ne. 0) der_flags(11) = der_flags(2) 
                 if (der_flags(11).ne. 0) der_flags(2)  = der_flags(11) 
                 if (der_flags(3) .ne. 0) der_flags(12) = der_flags(3) 
                 if (der_flags(12).ne. 0) der_flags(3)  = der_flags(12) 
               endif 

               if (atom(n) .eq. atom(r)) then 
                 if (der_flags(4) .ne. 0) der_flags(7) = der_flags(4) 
                 if (der_flags(7) .ne. 0) der_flags(4) = der_flags(7) 
                 if (der_flags(5) .ne. 0) der_flags(8) = der_flags(5) 
                 if (der_flags(8) .ne. 0) der_flags(5) = der_flags(8) 
                 if (der_flags(6) .ne. 0) der_flags(9) = der_flags(6) 
                 if (der_flags(9) .ne. 0) der_flags(6) = der_flags(9) 
               endif 

               if (atom(n) .eq. atom(s)) then 
                 if (der_flags(4)  .ne. 0) der_flags(10) = der_flags(4) 
                 if (der_flags(10) .ne. 0) der_flags(4)  = der_flags(10) 
                 if (der_flags(5)  .ne. 0) der_flags(11) = der_flags(5) 
                 if (der_flags(11) .ne. 0) der_flags(5)  = der_flags(11) 
                 if (der_flags(6)  .ne. 0) der_flags(12) = der_flags(6) 
                 if (der_flags(12) .ne. 0) der_flags(6)  = der_flags(12) 
               endif 

               if (atom(r) .eq. atom(s)) then 
                 if (der_flags(7)  .ne. 0) der_flags(10) = der_flags(7) 
                 if (der_flags(10) .ne. 0) der_flags(7)  = der_flags(10) 
                 if (der_flags(8)  .ne. 0) der_flags(11) = der_flags(8) 
                 if (der_flags(11) .ne. 0) der_flags(8)  = der_flags(11) 
                 if (der_flags(9)  .ne. 0) der_flags(12) = der_flags(9) 
                 if (der_flags(12) .ne. 0) der_flags(9)  = der_flags(12) 
               endif 

c----------------------------------------------------------------------------
c              Clear the output array.
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
 
               call move_integrals(out, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2)

c---------------------------------------------------------------------------
c   Contract with Intemediate block.
c---------------------------------------------------------------------------

               y = 0.0d0  
               do i = aa1, aa2 
               do j = bb1, bb2 
               do k = cc1, cc2 
               do l = dd1, dd2 
                  y = y + in(i,j,k,l)*out(i,j,k,l) 
               enddo  
               enddo  
               enddo  
               enddo  

c---------------------------------------------------------------------------
c   Sum into the hessian.
c---------------------------------------------------------------------------

              ihess = (iatom-1)*3 + icomponent
              jhess = (jatom-1)*3 + jcomponent

              if (ihess .eq. jhess)
     *        hess(jhess,ihess) = hess(jhess,ihess) + y

c---------------------------------------------------------------------------
c   Symmetrize the hessian if necassary??? Why???
c---------------------------------------------------------------------------
 
              if (ihess .ne. jhess) then 
                 y = 0.5d0*y  
                 hess(ihess,jhess) = hess(ihess,jhess) + y
                 hess(jhess,ihess) = hess(jhess,ihess) + y
              endif 

            endif

            enddo !ix = 1, 12 

777         continue 

            endif ! (intpkg .eq. flocke_package) then
         enddo   ! s
         enddo   ! r

         enddo   ! n
         enddo   ! m

      return
      end

      
      subroutine setflag(m,n,r,s,satom,nflags,iflags) 

      implicit none
      integer i,j,k,kmin,nmin,wmin,map(4),satom(4),xtype  
      integer ix,jx,ic,jc,jf,nflags,iflags(256,12) 
      integer m,n,r,s 

c------------------------------------------------------------------------------
c Determine which flags to set.  
c------------------------------------------------------------------------------

      nflags = 0 

c------------------------------------------------------------------------------
c (AB|CD) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .ne. satom(2)) .and. (satom(1) .ne. satom(3)) .and.  
     *    (satom(1) .ne. satom(4)) .and. (satom(2) .ne. satom(3)) .and.  
     *    (satom(2) .ne. satom(4)) .and. (satom(3) .ne. satom(4))) then   

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = ic, 4  
            do jx = ix, 3  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AA|BB) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(2)) .and. (satom(3) .eq. satom(4)) .and. 
     *    (satom(1) .ne. satom(3))) then 

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4   
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 11  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 1
     *                                      .and. ic .ne. 2) go to 11  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

11             continue 

            enddo 
            enddo 

         enddo 
         enddo 

c        endif 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AA|BC) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(2)) .and. (satom(3) .ne. satom(4)) .and. 
     *    (satom(1) .ne. satom(3)) .and. (satom(1) .ne. satom(4))) then 

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 12  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 1
     *                                      .and. ic .ne. 2) go to 12  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

12             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (BC|AA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(3) .eq. satom(4)) .and. (satom(1) .ne. satom(2)) .and. 
     *    (satom(3) .ne. satom(1)) .and. (satom(3) .ne. satom(2))) then 

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 13  
               if (satom(ic) .ne. satom(jc) .and. jc .ne. 3
     *                                      .and. jc .ne. 4) go to 13  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

13             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AB|AB) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(3)) .and. (satom(2) .eq. satom(4))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 14  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 1
     *                                      .and. ic .ne. 3) go to 14  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

14             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AB|BA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(4)) .and. (satom(2) .eq. satom(3))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 15  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 1
     *                                      .and. ic .ne. 4) go to 15  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

15             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AB|AC) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(3)) .and. (satom(2) .ne. satom(4)) .and.  
     *    (satom(1) .ne. satom(2)) .and. (satom(1) .ne. satom(4))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 24  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 1
     *                                      .and. ic .ne. 3) go to 24   

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

24             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AB|CA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(4)) .and. (satom(1) .ne. satom(2)) .and.    
     *    (satom(1) .ne. satom(3)) .and. (satom(2) .ne. satom(3))) then    

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 34  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 1
     *                                      .and. ic .ne. 4) go to 34   

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

34             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (BA|CA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(2) .eq. satom(4)) .and. (satom(1) .ne. satom(2)) .and.    
     *    (satom(1) .ne. satom(3)) .and. (satom(2) .ne. satom(3))) then    

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 44  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 2
     *                                      .and. ic .ne. 4) go to 44   

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

44             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (BA|AC) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(2) .eq. satom(3)) .and. (satom(1) .ne. satom(2)) .and.    
     *    (satom(1) .ne. satom(4)) .and. (satom(2) .ne. satom(4))) then    

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 54  
               if (satom(ic) .ne. satom(jc) .and. ic .ne. 2
     *                                      .and. ic .ne. 3) go to 54   

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

54             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c     go to 300 
200   continue 

c------------------------------------------------------------------------------
c (AA|AB) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(2)) .and. (satom(1) .eq. satom(3)) .and. 
     *    (satom(1) .ne. satom(4))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4   
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

c Omit duplicate loops. 
c --------------------- 

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 16  

               if ((satom(ic) .eq. satom(jc)) .and. (ix .eq. jx) .and. 
     *             ((ic .ne. 4)  .and. (ic .ne. 1))) go to 16  

               if ((satom(ic)  .ne. satom(jc))   .and. 
     *            ((ic .ne. 1) .or. (jc .ne. 4)) .and. 
     *            ((ic .ne. 4) .or. (jc .ne. 1))) go to 16  

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

16             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AA|BA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(2)) .and. (satom(1) .eq. satom(4)) .and. 
     *    (satom(1) .ne. satom(3))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

c Omit duplicate loops. 
c --------------------- 

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 17  

               if ((satom(ic) .eq. satom(jc)) .and. (ix .eq. jx) .and. 
     *             (ic .ne. 1)  .and. (ic .ne. 3)) go to 17  

               if ((ix .eq. jx) .and. (satom(ic) .ne. satom(jc)) .and. 
     *            ((ic .ne. 1) .or. (jc .ne. 3)) .and.  
     *            ((ic .ne. 3) .or. (jc .ne. 1))) go to 17 

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

17             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (AB|AA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(1) .eq. satom(3)) .and. (satom(1) .eq. satom(4)) .and. 
     *    (satom(1) .ne. satom(2))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

c Omit duplicate loops. 
c --------------------- 

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 18  

               if ((satom(ic) .eq. satom(jc)) .and. (ix .eq. jx) .and. 
     *             (ic .ne. 1)  .and. (ic .ne. 2)) go to 18  

               if ((ix .eq. jx) .and. (satom(ic) .ne. satom(jc)) .and. 
     *            ((ic .ne. 2) .or. (jc .ne. 1)) .and.  
     *            ((ic .ne. 1) .or. (jc .ne. 2))) go to 18 

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

18             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

c------------------------------------------------------------------------------
c (BA|AA) TYPE. 
c------------------------------------------------------------------------------

      if ((satom(2) .eq. satom(3)) .and. (satom(2) .eq. satom(4)) .and. 
     *    (satom(1) .ne. satom(2))) then  

         do ix = 1, 256 
         do jx = 1, 12 
            iflags(ix,jx) = 0 
         enddo 
         enddo 

         jf = 0 

         do ic = 1, 4 
         do ix = 1, 3 

            i = 3*(ic-1) + ix 

            do jc = 1, 4  
            do jx = 1, 3  

c Omit duplicate loops. 
c --------------------- 

               if (satom(ic) .eq. satom(jc) .and. ix .eq. jx
     *                                      .and. ic .ne. jc) go to 19  

               if ((satom(ic) .eq. satom(jc)) .and. (ix .eq. jx) .and. 
     *             (ic .ne. 1)  .and. (ic .ne. 4)) go to 19  

               if ((ix .eq. jx) .and. (satom(ic) .ne. satom(jc)) .and. 
     *            ((ic .ne. 4) .or. (jc .ne. 1)) .and.  
     *            ((ic .ne. 1) .or. (jc .ne. 4))) go to 19 

               j = 3*(jc-1) + jx 

               jf = jf + 1 

               iflags(jf,i) = 1 
               iflags(jf,j) = iflags(jf,j) + 1 

19             continue 

            enddo 
            enddo 

         enddo 
         enddo 

         nflags = jf 
         go to 300 
            
      endif 

300   continue 

c     if (nflags .eq. 0) then 
c        write(6,*) ' No flags set in compute_erd_2e_Sderv_integrals '
c        STOP 
c     endif  

      return 
      end 
 
