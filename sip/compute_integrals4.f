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
      subroutine compute_integrals4(a1,a2,b1,b2,c1,c2,d1,d2,scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 out, nsend, ca, nc1,nc2, nd1, nd2,
     *                 fa) 
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
      double precision fa(nc1:nc2,nc1:nc2)
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
            if (bmax .lt. itol) go to 30 
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
     *          intpkg .ne. gamess_derivative_package .and.
     *          nints .gt. 0) then
 
            if (jx .eq. 1) then 
               call scf_tp4(nc1,nc2,nd1,nd2,ca, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                             Fa)
            else 
               call scf_tp4a(nc1,nc2,nd1,nd2,ca,a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                             Fa)
            endif 

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
c
c ------------------------------------------------------------------- 
c  Used if 
c               IF nu == lambda
c               IF mu < nu
c               IF lambda < sigma
c               IF mu < sigma 
c ------------------------------------------------------------------- 
c 
      subroutine scf_tp4(nc1,nc2,nd1,nd2,ca,va1,va2,vb1,vb2,vc1,vc2,
     &                vd1,vd2,intblk, a1, a2, b1, b2, c1, c2, d1, d2,
     &                Fa)
      implicit none
      include 'int_gen_parms.h'
      integer va1, va2, vb1,vb2, vc1, vc2, vd1, vd2
      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer a,b,c,d
      integer nc1,nc2,nd1,nd2,p  
      integer drange1, drange2
      integer crange1, crange2
      integer brange1, brange2
      integer arange1, arange2

      double precision ca(nc1:nc2,nd1:nd2)
      double precision Fa(nc1:nc2,nc1:nc2)
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision dtemp  

      drange1 = max(vd1, d1)
      drange2 = min(vd2, d2)
      crange1 = max(vc1, c1)
      crange2 = min(vc2, c2)
      brange1 = max(vb1, b1)
      brange2 = min(vb2, b2)
      arange1 = max(va1, a1)
      arange2 = min(va2, a2)

      do d = drange1, drange2
      do c = crange1, crange2
         dtemp = 0.0 
         do p = 1, nalpha_occupied ! nd1,nd2
            dtemp = dtemp + ca(c,p)*ca(d,p) 
         enddo 
         dtemp = 4.0d0*dtemp  
         do b = brange1, brange2
         do a = arange1, arange2
            Fa(a,b) = Fa(a,b) + intblk(a,b,c,d)*dtemp 
            Fa(b,a) = Fa(b,a) + intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      do b = brange1, brange2
      do a = arange1, arange2
         dtemp = 0.0 
         do p = 1, nalpha_occupied ! nd1,nd2
            dtemp = dtemp + ca(a,p)*ca(b,p) 
         enddo 
         dtemp = 4.0d0*dtemp  
         do d = drange1, drange2
         do c = crange1, crange2
            Fa(c,d) = Fa(c,d) + intblk(a,b,c,d)*dtemp 
            Fa(d,c) = Fa(d,c) + intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      do d = drange1, drange2
      do b = brange1, brange2
         dtemp = 0.0 
         do p = 1, nalpha_occupied ! nd1,nd2
            dtemp = dtemp + ca(b,p)*ca(d,p) 
         enddo 
         dtemp = -1.0*dtemp 
         do c = crange1, crange2
         do a = arange1, arange2
            Fa(a,c) = Fa(a,c) + intblk(a,b,c,d)*dtemp 
            Fa(c,a) = Fa(c,a) + intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      do c = crange1, crange2
      do b = brange1, brange2
         dtemp = 0.0 
         do p = 1, nalpha_occupied ! nd1,nd2
            dtemp = dtemp + ca(b,p)*ca(c,p) 
         enddo 
         dtemp = -1.0*dtemp 
         do d = drange1, drange2
         do a = arange1, arange2
            Fa(a,d) = Fa(a,d) + intblk(a,b,c,d)*dtemp 
            Fa(d,a) = Fa(d,a) + intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      do d = drange1, drange2
      do a = arange1, arange2
         dtemp = 0.0 
         do p = 1, nalpha_occupied ! nd1,nd2
            dtemp = dtemp + ca(a,p)*ca(d,p) 
         enddo 
         dtemp = -1.0*dtemp 
         do c = crange1, crange2
         do b = brange1, brange2
            Fa(b,c) = Fa(b,c) + intblk(a,b,c,d)*dtemp 
            Fa(c,b) = Fa(c,b) + intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      do c = crange1, crange2
      do a = arange1, arange2
         dtemp = 0.0 
         do p = 1, nalpha_occupied ! nd1,nd2
            dtemp = dtemp + ca(a,p)*ca(c,p) 
         enddo 
         dtemp = -1.0*dtemp 
         do d = drange1, drange2
         do b = brange1, brange2
            Fa(b,d) = Fa(b,d) + intblk(a,b,c,d)*dtemp 
            Fa(d,b) = Fa(d,b) + intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      return
      end
c
c ------------------------------------------------------------------- 
c  Used if 
c               IF nu == lambda
c               IF mu < nu
c               IF lambda < sigma
c               IF mu < sigma 
c ------------------------------------------------------------------- 
c 
      subroutine scf_tp4a(nc1,nc2,nd1,nd2,ca,va1,va2,vb1,vb2,vc1,vc2,
     &                vd1,vd2,intblk, a1, a2, b1, b2, c1, c2, d1, d2,
     &                Fa)
      implicit none
      include 'int_gen_parms.h'
      integer va1, va2, vb1,vb2, vc1, vc2, vd1, vd2
      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer a,b,c,d
      integer nc1,nc2,nd1,nd2,p  
      integer drange1, drange2
      integer crange1, crange2
      integer brange1, brange2
      integer arange1, arange2
      integer itype, sym

      double precision ca(nc1:nc2,nd1:nd2)
      double precision Fa(nc1:nc2,nc1:nc2)
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision dtemp, factor, dcrit   

      drange1 = max(vd1, d1)
      drange2 = min(vd2, d2)
      crange1 = max(vc1, c1)
      crange2 = min(vc2, c2)
      brange1 = max(vb1, b1)
      brange2 = min(vb2, b2)
      arange1 = max(va1, a1)
      arange2 = min(va2, a2)

      dcrit = 0.0d0 

      factor = 4.0
      itype  = 1
      call fdmult1a(arange1,arange2,brange1,brange2,
     *              crange1,crange2,drange1,drange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype)

      factor = 4.0
      itype  = 2
      call fdmult1a(crange1,crange2,drange1,drange2,
     *              arange1,arange2,brange1,brange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

      itype  = 3
      factor =-1.0
      call fdmult1a(arange1,arange2,crange1,crange2,
     *              brange1,brange2,drange1,drange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

      itype = 4 
      factor =-1.0 
      call fdmult1a(arange1,arange2,drange1,drange2,
     *              crange1,crange2,brange1,brange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

      itype = 5
      factor =-1.0
      call fdmult1a(crange1,crange2,brange1,brange2,
     *              arange1,arange2,drange1,drange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype)

      itype = 6 
      factor =-1.0 
      call fdmult1a(brange1,brange2,drange1,drange2,
     *              arange1,arange2,crange1,crange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do d = drange1, drange2
c     do c = crange1, crange2
c        dtemp = ca(c,d) 
c        dtemp = 4.0d0*dtemp  
c        do b = brange1, brange2
c        do a = arange1, arange2
c           Fa(a,b) = Fa(a,b) + intblk(a,b,c,d)*dtemp 
c           Fa(b,a) = Fa(b,a) + intblk(a,b,c,d)*dtemp 
c        enddo
c        enddo
c     enddo
c     enddo

c     do b = brange1, brange2
c     do a = arange1, arange2
c        dtemp = ca(a,b) 
c        dtemp = 4.0d0*dtemp  
c        do d = drange1, drange2
c        do c = crange1, crange2
c           Fa(c,d) = Fa(c,d) + intblk(a,b,c,d)*dtemp 
c           Fa(d,c) = Fa(d,c) + intblk(a,b,c,d)*dtemp 
c        enddo
c        enddo
c     enddo
c     enddo

c     do d = drange1, drange2
c     do b = brange1, brange2
c        dtemp = ca(b,d) 
c        dtemp = -1.0*dtemp 
c        do c = crange1, crange2
c        do a = arange1, arange2
c           Fa(a,c) = Fa(a,c) + intblk(a,b,c,d)*dtemp 
c           Fa(c,a) = Fa(c,a) + intblk(a,b,c,d)*dtemp 
c        enddo
c        enddo
c     enddo
c     enddo

c     do c = crange1, crange2
c     do b = brange1, brange2
c        dtemp = ca(b,c) 
c        dtemp = -1.0*dtemp 
c        do d = drange1, drange2
c        do a = arange1, arange2
c           Fa(a,d) = Fa(a,d) + intblk(a,b,c,d)*dtemp 
c           Fa(d,a) = Fa(d,a) + intblk(a,b,c,d)*dtemp 
c        enddo
c        enddo
c     enddo
c     enddo

c     do d = drange1, drange2
c     do a = arange1, arange2
c        dtemp = ca(a,d) 
c        dtemp = -1.0*dtemp 
c        do c = crange1, crange2
c        do b = brange1, brange2
c           Fa(b,c) = Fa(b,c) + intblk(a,b,c,d)*dtemp 
c           Fa(c,b) = Fa(c,b) + intblk(a,b,c,d)*dtemp 
c        enddo
c        enddo
c     enddo
c     enddo

c     do c = crange1, crange2
c     do a = arange1, arange2
c        dtemp = ca(a,c) 
c        dtemp = -1.0*dtemp 
c        do d = drange1, drange2
c        do b = brange1, brange2
c           Fa(b,d) = Fa(b,d) + intblk(a,b,c,d)*dtemp 
c           Fa(d,b) = Fa(d,b) + intblk(a,b,c,d)*dtemp 
c        enddo
c        enddo
c     enddo
c     enddo

      return
      end
 
