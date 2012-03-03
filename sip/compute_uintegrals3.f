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
      subroutine compute_uintegrals3(a1,a2,b1,b2,c1,c2,d1,d2,scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 out, nsend, ca, cb, nc1,nc2, nd1, nd2,
     *                 fa,fb) 
c---------------------------------------------------------------------------
c   The block of integrals (a1:a2,b1:b2,c1:c2,d1:d2) is computed for the 
c   following 'types' of integrals based on atomic labels.
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
      integer nc1, nc2, nd1, nd2 

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
      double precision ca(nc1:nc2,nd1:nd2)
      double precision cb(nc1:nc2,nd1:nd2)
      double precision fa(nc1:nc2,nc1:nc2)
      double precision fb(nc1:nc2,nc1:nc2)
      double precision scr(*)   
      integer iscr(*)

      integer ccbeg(*), ccend(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      integer*8 ccbeg_pack64(max_dim_coeff), ccend_pack64(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      integer*8 arg64(25)

      common /Imax_com/sz_max(max_shells,max_shells), delta 
      double precision sz_max, delta
      double precision itol, bmax, dtemp, emax    

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter 

      save me,alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *     ccbeg_pack64, ccend_pack64

      call mpi_comm_rank(mpi_comm_world, me, ierr)
c      print *,'Task ',me,' computing integrals for ',a1,a2,b1,b2,
c     *     c1,c2,d1,d2
c      call c_flush_stdout()

      adim = a2-a1+1
      bdim = b2-b1+1
      cdim = c2-c1+1
      ddim = d2-d1+1 
      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

c Set the integral tolerance 

      call set_itol(delta,itol)

c     write(6,*) ' DELTA = ', delta, itol  
c     itol = delta ! 1.0d-9 
  
      nsend = adim*bdim*cdim*ddim
      if (nsend .lt. 0) then
         print *,'ERROR IN INTEGRAL WORKER ',me,' nsend = ',nsend
         print *,'adim,bdim,cdim,ddim = ',adim,bdim,cdim,ddim
         call mpi_abort(mpi_comm_world, ierr)
      endif

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
            if (s .eq. 1) then
               dd1 = 1
            else
               dd1 = end_nfps(s-1) + 1
            endif
            dd2 = end_nfps(s)

c-----------------------------------------------------------------------
c   Determine the largest density element.
c-----------------------------------------------------------------------

            emax = 1.0 
c           do d = aa1, aa2  
c           do b = cc1, cc2  
c              dtemp = 0.0 
c              do a = 1, nalpha_occupied ! nd1,nd2
c                  dtemp = dtemp + ca(b,a)*ca(d,a) 
c              enddo 
c              dtemp = dabs(dtemp) 
c              if (dtemp .gt. emax) emax = dtemp  
c           enddo
c           enddo
c           do d = bb1, bb2  
c           do b = dd1, dd2  
c              dtemp = 0.0 
c              do a = 1, nalpha_occupied ! nd1,nd2
c                  dtemp = dtemp + ca(b,a)*ca(d,a) 
c              enddo 
c              dtemp = dabs(dtemp) 
c              if (dtemp .gt. emax) emax = dtemp  
c           enddo
c           enddo

            bmax = sz_max(m,n)*sz_max(r,s) 
            bmax = dsqrt(bmax)*emax  
c           if ((sz_max(m,r) .gt. 0.0) .or. 
c    &           sz_max(n,s) .gt. 0.0)  
c    &      write(6,*) m,n,r,s,sz_max(m,r),sz_max(n,s),bmax 
c           if (bmax .lt. itol) go to 30 
            if (intpkg .eq. flocke_package) then
               x4 = coords(1,s)
               y4 = coords(2,s)
               z4 = coords(3,s)
               call pack_coeffs(alphas, ixalpha, coeffs, ixpcoef, 
     *                          ncfps, npfps, m, n, 
     *                          r, s, alpha_pack, nalpha_pack, 
     *                          pcoeff_pack, npcoeff_pack, 
     *                          ccbeg, ccend, indx_cc,
     *                          ccbeg_pack, ccend_pack, 
     *                          ccbeg_pack64, ccend_pack64)

c---------------------------------------------------------------------------
c   Calling sequence for ERD version 2.
c---------------------------------------------------------------------------

               ncsum = ncfps(m) + ncfps(n) + ncfps(r) + ncfps(s)

               if (aces64) then
                  arg64(1) = intmax
                  arg64(2) = zmax
                  arg64(3) = nalpha_pack
                  arg64(4) = npcoeff_pack
                  arg64(5) = ncsum
                  arg64(6) = ncfps(m)
                  arg64(7) = ncfps(n)
                  arg64(8) = ncfps(r)
                  arg64(9) = ncfps(s)
                  arg64(10) = npfps(m)
                  arg64(11) = npfps(n)
                  arg64(12) = npfps(r)
                  arg64(13) = npfps(s)
                  arg64(18) = ivangmom(m)
                  arg64(19) = ivangmom(n)
                  arg64(20) = ivangmom(r)
                  arg64(21) = ivangmom(s)
                  call ERD__GENER_ERI_BATCH(arg64(1), arg64(2),
     *                arg64(3), arg64(4), arg64(5),
     *                arg64(6), arg64(7), arg64(8), arg64(9),
     *                arg64(10), arg64(11), arg64(12), arg64(13),
     *                arg64(18), arg64(19), arg64(20), arg64(21),
     *                x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, alpha_pack,
     *                pcoeff_pack, ccbeg_pack64, ccend_pack64,
     *                l8spherical, l8true, iscr, arg64(22), 
     *                arg64(23), scr)    
                  nints = arg64(22)
                  nfirst = arg64(23)
               else
                  call ERD__GENER_ERI_BATCH(intmax, zmax,
     *                nalpha_pack, npcoeff_pack, ncsum, 
     *                ncfps(m),ncfps(n), ncfps(r), ncfps(s),
     *                npfps(m),npfps(n), npfps(r), npfps(s),
     *                ivangmom(m), ivangmom(n), 
     *                ivangmom(r), ivangmom(s), x1,y1,z1,
     *                x2,y2,z2,x3,y3,z3,x4,y4,z4, alpha_pack,
     *                pcoeff_pack, ccbeg_pack, ccend_pack,
     *                spherical, .true., iscr, nints, 
     *                nfirst, scr)    
               endif
            else if (intpkg .eq. gamess_package .or.
     *               intpkg .eq. gamess_derivative_package) then
            else 
               print *,'Error: Invalid integral package: ',intpkg
               call abort_job()
            endif

c---------------------------------------------------------------------------
c   Move the integrals into the output block.  (For the GAMESS integral 
c   package, this has been already performed.)
c---------------------------------------------------------------------------

           if (intpkg .ne. gamess_package .and.
     *         intpkg .ne. gamess_derivative_package .and.
     *         nints .gt. 0) then
 
              call uscf_tp3a(nc1,nc2,nd1,nd2,ca,cb,a1,a2,b1,b2,
     *                             c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                             Fa,Fb)

            endif

30       continue 
         enddo   ! s
20       continue 
         enddo   ! r

10       continue 
         enddo   ! n
100      continue 
         enddo   ! m

c     write(6,*) ' Fock_a inside' 
c     write(6,*) ' -------------' 
c     do b = nc1, nc2
c     do a = nc1, nc2
c        write(6,*) a, b, Fa(a,b) 
c     enddo
c     enddo 

      return
      end
 
c ------------------------------------------------------------------- 
c  Used if 
c ------------------------------------------------------------------- 
c 
      subroutine uscf_tp3a(nc1,nc2,nd1,nd2,ca,cb,va1,va2,vb1,vb2,
     &                vc1,vc2,
     &                vd1,vd2,intblk, a1, a2, b1, b2, c1, c2, d1, d2,
     &                Fa,Fb)
      implicit none
      include 'int_gen_parms.h'
      integer va1, va2, vb1,vb2, vc1, vc2, vd1, vd2
      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer a,b,c,d,sym 
      integer nc1,nc2,nd1,nd2,p  
      integer drange1, drange2
      integer crange1, crange2
      integer brange1, brange2
      integer arange1, arange2
      integer delta_a,itype  

      double precision ca(nc1:nc2,nd1:nd2)
      double precision cb(nc1:nc2,nd1:nd2)
      double precision Fa(nc1:nc2,nc1:nc2)
      double precision Fb(nc1:nc2,nc1:nc2)
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision dtemp, etemp, dcrit, factor  

      dcrit = 0.0 
      drange1 = max(vd1, d1)
      drange2 = min(vd2, d2)
      crange1 = max(vc1, c1)
      crange2 = min(vc2, c2)
      brange1 = max(vb1, b1)
      brange2 = min(vb2, b2)
      arange1 = max(va1, a1)
      arange2 = min(va2, a2)

c ----------------------------------------------------------------- 

      do d = drange1, drange2
      do c = crange1, crange2
         dtemp = ca(c,d) + cb(c,d)  
         dtemp = 2.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do b = brange1, brange2
            do a = arange1, arange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fa(a,b) = Fa(a,b) + etemp ! intblk(a,b,c,d)*dtemp 
               Fa(b,a) = Fa(b,a) + etemp ! intblk(a,b,c,d)*dtemp 

               Fb(a,b) = Fb(a,b) + etemp ! intblk(a,b,c,d)*dtemp 
               Fb(b,a) = Fb(b,a) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do a = arange1, arange2
      do b = brange1, brange2
         dtemp = ca(a,b) + cb(a,b)  
         dtemp = 2.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do c = crange1, crange2
            do d = drange1, drange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fa(c,d) = Fa(c,d) + etemp ! intblk(a,b,c,d)*dtemp 
               Fa(d,c) = Fa(d,c) + etemp ! intblk(a,b,c,d)*dtemp 

               Fb(c,d) = Fb(c,d) + etemp ! intblk(a,b,c,d)*dtemp 
               Fb(d,c) = Fb(d,c) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do b = brange1, brange2
      do d = drange1, drange2
         dtemp = ca(b,d) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do a = arange1, arange2
            do c = crange1, crange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fa(a,c) = Fa(a,c) + etemp ! intblk(a,b,c,d)*dtemp 
               Fa(c,a) = Fa(c,a) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do b = brange1, brange2
      do c = crange1, crange2
         dtemp = ca(b,c) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do a = arange1, arange2
            do d = drange1, drange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fa(a,d) = Fa(a,d) + etemp ! intblk(a,b,c,d)*dtemp 
               Fa(d,a) = Fa(d,a) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do a = arange1, arange2
      do d = drange1, drange2
         dtemp = ca(a,d) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do b = brange1, brange2
            do c = crange1, crange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fa(b,c) = Fa(b,c) + etemp ! intblk(a,b,c,d)*dtemp 
               Fa(c,b) = Fa(c,b) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do a = arange1, arange2
      do c = crange1, crange2
         dtemp = ca(a,c) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do b = brange1, brange2
            do d = drange1, drange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fa(b,d) = Fa(b,d) + etemp ! intblk(a,b,c,d)*dtemp 
               Fa(d,b) = Fa(d,b) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do b = brange1, brange2
      do d = drange1, drange2
         dtemp = cb(b,d) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do a = arange1, arange2
            do c = crange1, crange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fb(a,c) = Fb(a,c) + etemp ! intblk(a,b,c,d)*dtemp 
               Fb(c,a) = Fb(c,a) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do b = brange1, brange2
      do c = crange1, crange2
         dtemp = cb(b,c) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do a = arange1, arange2
            do d = drange1, drange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fb(a,d) = Fb(a,d) + etemp ! intblk(a,b,c,d)*dtemp 
               Fb(d,a) = Fb(d,a) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do a = arange1, arange2
      do d = drange1, drange2
         dtemp = cb(a,d) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do b = brange1, brange2
            do c = crange1, crange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fb(b,c) = Fb(b,c) + etemp ! intblk(a,b,c,d)*dtemp 
               Fb(c,b) = Fb(c,b) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      do a = arange1, arange2
      do c = crange1, crange2
         dtemp = cb(a,c) 
         dtemp = -1.0d0*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            do b = brange1, brange2
            do d = drange1, drange2
               etemp = dtemp*intblk(a,b,c,d) 
               Fb(b,d) = Fb(b,d) + etemp ! intblk(a,b,c,d)*dtemp 
               Fb(d,b) = Fb(d,b) + etemp ! intblk(a,b,c,d)*dtemp 
            enddo
            enddo
         endif 
      enddo
      enddo

      return
      end
c
