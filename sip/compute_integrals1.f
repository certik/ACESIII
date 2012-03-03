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
      subroutine compute_integrals1(a1,a2,b1,b2,c1,c2,d1,d2,scr,
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
               call scf_tp1(nc1,nc2,nd1,nd2,ca, a1,a2,b1,b2,c1,c2,d1,d2,
     *                             scr(nfirst), 
     *                             aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                             Fa)
             else 
               call scf_tp1a(nc1,nc2,nd1,nd2,ca,a1,a2,b1,b2,c1,c2,d1,d2,
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
 
      subroutine set_itol(delta,itol) 

      implicit none
      double precision delta, itol  

      itol = 1.0d-5 
      if (delta .le. 0.01)   itol = 1.0d-6  
      if (delta .le. 0.001)  itol = 1.0d-7  
      if (delta .le. 0.0001) itol = 1.0d-8 
      if (delta .le. 0.00001) itol = 1.0d-9 
c
c The screening tolerance is set at 10^-12 for now. VFL 
c 
      itol = 1.0d-12 

      return 
      end 
c
c ------------------------------------------------------------------- 
c  Used if 
c                IF m  < n
c                IF l  < s
c                IF m  < l
c                IF n  != s
c                IF n  != l
c                IF m  != s  
c ------------------------------------------------------------------- 
c 
      subroutine scf_tp1a(nc1,nc2,nd1,nd2,ca,va1,va2,vb1,vb2,vc1,vc2,
     &                vd1,vd2,intblk, a1, a2, b1, b2, c1, c2, d1, d2,
     &                Fa)
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
      double precision Fa(nc1:nc2,nc1:nc2)
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

      factor = 4.0 
      itype  = 1 
      call fdmult1a(arange1,arange2,brange1,brange2, 
     *              crange1,crange2,drange1,drange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,  
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do d = drange1, drange2
c     do c = crange1, crange2
c        dtemp = ca(c,d)  
c        dtemp = 4.0d0*dtemp  
c        if (dabs(dtemp) .gt. dcrit) then 
c           do b = brange1, brange2
c           do a = arange1, arange2
c              etemp = dtemp*intblk(a,b,c,d) 
c              Fa(a,b) = Fa(a,b) + etemp ! intblk(a,b,c,d)*dtemp 
c              Fa(b,a) = Fa(b,a) + etemp ! intblk(a,b,c,d)*dtemp 
c           enddo
c           enddo
c        endif 
c     enddo
c     enddo

      factor = 4.0 
      itype  = 2 
      call fdmult1a(crange1,crange2,drange1,drange2, 
     *              arange1,arange2,brange1,brange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,  
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do b = brange1, brange2
c     do a = arange1, arange2
c        dtemp = ca(a,b) 
c        dtemp = 4.0d0*dtemp  
c        if (dabs(dtemp) .gt. dcrit) then 
c           do d = drange1, drange2
c           do c = crange1, crange2
c              etemp = intblk(a,b,c,d)*dtemp 
c              Fa(c,d) = Fa(c,d) + etemp ! intblk(a,b,c,d)*dtemp 
c              Fa(d,c) = Fa(d,c) + etemp ! intblk(a,b,c,d)*dtemp 
c           enddo
c           enddo
c        endif 
c     enddo
c     enddo

      itype  = 3 
      factor =-1.0 
      call fdmult1a(arange1,arange2,crange1,crange2, 
     *              brange1,brange2,drange1,drange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,  
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do d = drange1, drange2
c     do b = brange1, brange2
c        dtemp = ca(b,d) 
c        if (dabs(dtemp) .gt. dcrit) then 
c           dtemp = -1.0*dtemp 
c           do c = crange1, crange2
c           do a = arange1, arange2
c              etemp = intblk(a,b,c,d)*dtemp 
c              Fa(a,c) = Fa(a,c) + etemp ! intblk(a,b,c,d)*dtemp 
c              Fa(c,a) = Fa(c,a) + etemp ! intblk(a,b,c,d)*dtemp 
c           enddo
c           enddo
c        endif 
c     enddo
c     enddo

      itype = 4 
      factor =-1.0 
      call fdmult1a(arange1,arange2,drange1,drange2, 
     *              crange1,crange2,brange1,brange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,  
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do c = crange1, crange2
c     do b = brange1, brange2
c        dtemp = ca(b,c) 
c        if (dabs(dtemp) .gt. dcrit) then 
c           dtemp = -1.0*dtemp 
c           do d = drange1, drange2
c           do a = arange1, arange2
c              etemp = intblk(a,b,c,d)*dtemp 
c              Fa(a,d) = Fa(a,d) + etemp ! intblk(a,b,c,d)*dtemp 
c              Fa(d,a) = Fa(d,a) + etemp ! intblk(a,b,c,d)*dtemp 
c           enddo
c           enddo
c        endif 
c     enddo
c     enddo

      itype = 5 
      factor =-1.0 
      call fdmult1a(crange1,crange2,brange1,brange2, 
     *              arange1,arange2,drange1,drange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,  
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do d = drange1, drange2
c     do a = arange1, arange2
c        dtemp = ca(a,d) 
c        if (dabs(dtemp) .gt. dcrit) then 
c           dtemp = -1.0d0*dtemp  
c           do b = brange1, brange2
c           do c = crange1, crange2
c              etemp = intblk(a,b,c,d)*dtemp 
c              Fa(b,c) = Fa(b,c) + etemp ! intblk(a,b,c,d)*dtemp 
c              Fa(c,b) = Fa(c,b) + etemp ! intblk(a,b,c,d)*dtemp 
c           enddo
c           enddo
c        endif 
c     enddo
c     enddo

      itype = 6 
      factor =-1.0 
      call fdmult1a(brange1,brange2,drange1,drange2, 
     *              arange1,arange2,crange1,crange2,
     *              a1,a2,b1,b2,c1,c2,d1,d2,intblk,factor,  
     *              sym,nc1,nc2,nd1,nd2,ca,fa,dcrit,itype) 

c     do c = crange1, crange2
c     do a = arange1, arange2
c        dtemp = ca(a,c) 
c        if (dabs(dtemp) .gt. dcrit) then 
c           dtemp = -1.0d0*dtemp  
c           do b = brange1, brange2
c           do d = drange1, drange2
c              etemp = intblk(a,b,c,d)*dtemp 
c              Fa(b,d) = Fa(b,d) + etemp ! intblk(a,b,c,d)*dtemp 
c              Fa(d,b) = Fa(d,b) + etemp ! intblk(a,b,c,d)*dtemp 
c           enddo
c           enddo
c        endif 
c     enddo
c     enddo

      return
      end
c
c ------------------------------------------------------------------- 
c  Used if 
c                IF m  < n
c                IF l  < s
c                IF m  < l
c                IF n  != s
c                IF n  != l
c                IF m  != s  
c ------------------------------------------------------------------- 
c 
      subroutine scf_tp1(nc1,nc2,nd1,nd2,ca,va1,va2,vb1,vb2,vc1,vc2,
     &                vd1,vd2,intblk, a1, a2, b1, b2, c1, c2, d1, d2,
     &                Fa)
      implicit none
      include 'int_gen_parms.h'
      include 'mpif.h'
      integer va1, va2, vb1,vb2, vc1, vc2, vd1, vd2
      integer a1, a2, b1, b2, c1, c2, d1, d2
      integer a,b,c,d
      integer nc1,nc2,nd1,nd2,p, istride, sym   
      integer drange1, drange2
      integer crange1, crange2
      integer brange1, brange2
      integer arange1, arange2

      double precision t1, t2, t3  
      integer me 

      double precision ca(nc1:nc2,nd1:nd2)
      double precision Fa(nc1:nc2,nc1:nc2)
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision dtemp, etemp, etemp1, etemp2, etemp3, etemp4, 
     *                 etemp5, dcrit  

      dcrit = 1.0d-12 
      istride = 5 
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
         dtemp = 4.0*dtemp  
         do b = brange1, brange2
         do a = arange1, arange2, 1 
                etemp1 = dtemp*intblk(a,b,c,d) 
                Fa(a,b) = Fa(a,b) + etemp1 
                Fa(b,a) = Fa(b,a) + etemp1
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
            etemp = intblk(a,b,c,d)*dtemp 
            Fa(c,d) = Fa(c,d) + etemp ! intblk(a,b,c,d)*dtemp 
            Fa(d,c) = Fa(d,c) + etemp ! intblk(a,b,c,d)*dtemp 
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
            etemp = intblk(a,b,c,d)*dtemp 
            Fa(a,c) = Fa(a,c) + etemp ! intblk(a,b,c,d)*dtemp 
            Fa(c,a) = Fa(c,a) + etemp ! intblk(a,b,c,d)*dtemp 
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
            etemp = intblk(a,b,c,d)*dtemp 
            Fa(a,d) = Fa(a,d) + etemp ! intblk(a,b,c,d)*dtemp 
            Fa(d,a) = Fa(d,a) + etemp ! intblk(a,b,c,d)*dtemp 
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
         dtemp = -1.0d0*dtemp  
         do b = brange1, brange2
         do c = crange1, crange2
            etemp = intblk(a,b,c,d)*dtemp 
            Fa(b,c) = Fa(b,c) + etemp ! intblk(a,b,c,d)*dtemp 
            Fa(c,b) = Fa(c,b) + etemp ! intblk(a,b,c,d)*dtemp 
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
         dtemp = -1.0d0*dtemp  
         do b = brange1, brange2
         do d = drange1, drange2
            etemp = intblk(a,b,c,d)*dtemp 
            Fa(b,d) = Fa(b,d) + etemp ! intblk(a,b,c,d)*dtemp 
            Fa(d,b) = Fa(d,b) + etemp ! intblk(a,b,c,d)*dtemp 
         enddo
         enddo
      enddo
      enddo

      return
      end

      subroutine fdmult1a(arange1,arange2,brange1,brange2, 
     *                    crange1,crange2,drange1,drange2,
     *                    a1,a2,b1,b2,c1,c2,d1,d2,intblk, 
     *                    factor,sym,nc1,nc2,nd1,nd2,
     *                    ca,fa,dcrit,itype) 
      implicit none
      include 'int_gen_parms.h'
      integer arange1,arange2,brange1,brange2,nc1,nc2,sym
      integer crange1,crange2,drange1,drange2,nd1,nd2
      integer a, b, c, d, p, delta_a, itype  
      integer a1,a2,b1,b2,c1,c2,d1,d2 
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2) 
      double precision ca(nc1:nc2,nd1:nd2)
      double precision fa(nc1:nc2,nc1:nc2)
      double precision Xa(nc1:nc2,nc1:nc2)
      double precision factor, dcrit, dtemp, etemp, etemp1, 
     *                 etemp2, etemp3, etemp4, etemp5  

      delta_a = arange2 - arange1 

      do d = drange1, drange2
      do c = crange1, crange2
         dtemp = ca(c,d)  
         dtemp = factor*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            if (itype .eq. 1) then 
               do b = brange1, brange2
               do p = arange1, arange2, 5 
                  if (arange2-p .gt. 4) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                      Xa(p+3,b) = intblk(p+3,b,c,d) 
                      Xa(p+4,b) = intblk(p+4,b,c,d) 
                      go to 1 
                  endif 
                  if (arange2-p .eq. 4) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                      Xa(p+3,b) = intblk(p+3,b,c,d) 
                      Xa(p+4,b) = intblk(p+4,b,c,d) 
                  endif 
                  if (arange2-p .eq. 3) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                      Xa(p+3,b) = intblk(p+3,b,c,d) 
                  endif 
                  if (arange2-p .eq. 2) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                  endif 
                  if (arange2-p .eq. 1) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                  endif 
                  if (arange2-p .eq. 0) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                  endif 
1              continue 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 2) then 
               do b = brange1, brange2
               do p = arange1, arange2 
                  Xa(p,b) = intblk(c,d,p,b) 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 3) then 
               do b = brange1, brange2
               do p = arange1, arange2, 5  
                  if (arange2-p .ge. 4) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                     Xa(p+2,b) = intblk(p+2,c,b,d) 
                     Xa(p+3,b) = intblk(p+3,c,b,d) 
                     Xa(p+4,b) = intblk(p+4,c,b,d) 
                     go to 2 
                  endif 
                  if (arange2-p .eq. 3) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                     Xa(p+2,b) = intblk(p+2,c,b,d) 
                     Xa(p+3,b) = intblk(p+3,c,b,d) 
                  endif 
                  if (arange2-p .eq. 2) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                     Xa(p+2,b) = intblk(p+2,c,b,d) 
                  endif 
                  if (arange2-p .eq. 1) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                  endif 
                  if (arange2-p .eq. 0) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                  endif 
2              continue 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 4) then 
               do b = brange1, brange2
               do p = arange1, arange2, 5  
                  if (arange2-p .ge. 4) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     Xa(p+2,b) = intblk(p+2,d,c,b) 
                     Xa(p+3,b) = intblk(p+3,d,c,b) 
                     Xa(p+4,b) = intblk(p+4,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 3) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     Xa(p+2,b) = intblk(p+2,d,c,b) 
                     Xa(p+3,b) = intblk(p+3,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 2) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     Xa(p+2,b) = intblk(p+2,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 1) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 0) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     go to 3 
                  endif 
3              continue 
               enddo  
               enddo  
            endif 
            if (itype .eq. 5) then 
               do b = brange1, brange2
               do p = arange1, arange2 
                  Xa(p,b) = intblk(c,b,p,d) 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 6) then 
               do b = brange1, brange2
               do p = arange1, arange2 
                  Xa(p,b) = intblk(c,p,d,b) 
               enddo  
               enddo  
            go to 44 
            endif 

44          continue 

            do b = brange1, brange2

               a = arange1 

               if (delta_a .eq. 0) then  
                     etemp1  = dtemp*Xa(a,b) 
                     Fa(a,b) = Fa(a,b) + etemp1 
                     Fa(b,a) = Fa(b,a) + etemp1
                  go to 33 
               endif  

               if (delta_a .eq. 1) then  
                     etemp   = dtemp*Xa(a,b) 
                     etemp1  = dtemp*Xa(a+1,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(b,a) = Fa(b,a) + etemp 

                     Fa(a+1,b) = Fa(a+1,b) + etemp1 
                     Fa(b,a+1) = Fa(b,a+1) + etemp1
                  go to 33 
               endif  

               if (delta_a .eq. 2) then  
                     etemp  = dtemp*Xa(a,b) 
                     etemp1 = dtemp*Xa(a+1,b) 
                     etemp2 = dtemp*Xa(a+2,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(b,a) = Fa(b,a) + etemp 

                     Fa(a+1,b) = Fa(a+1,b) + etemp1  
                     Fa(b,a+1) = Fa(b,a+1) + etemp1 

                     Fa(a+2,b) = Fa(a+2,b) + etemp2 
                     Fa(b,a+2) = Fa(b,a+2) + etemp2
                  go to 33 
               endif  

               if (delta_a .eq. 3) then  
                     etemp = dtemp*Xa(a,b) 
                     etemp1 = dtemp*Xa(a+1,b) 
                     etemp2 = dtemp*Xa(a+2,b) 
                     etemp3 = dtemp*Xa(a+3,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(b,a) = Fa(b,a) + etemp 

                     Fa(a+1,b) = Fa(a+1,b) + etemp1 
                     Fa(b,a+1) = Fa(b,a+1) + etemp1

                     Fa(a+2,b) = Fa(a+2,b) + etemp2 
                     Fa(b,a+2) = Fa(b,a+2) + etemp2

                     Fa(a+3,b) = Fa(a+3,b) + etemp3 
                     Fa(b,a+3) = Fa(b,a+3) + etemp3
                  go to 33 
               endif  

               if (delta_a .eq. 4) then  
                     etemp = dtemp*Xa(a,b) 
                     etemp1 = dtemp*Xa(a+1,b) 
                     etemp2 = dtemp*Xa(a+2,b) 
                     etemp3 = dtemp*Xa(a+3,b) 
                     etemp4 = dtemp*Xa(a+4,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(b,a) = Fa(b,a) + etemp 

                     Fa(a+1,b) = Fa(a+1,b) + etemp1 
                     Fa(b,a+1) = Fa(b,a+1) + etemp1

                     Fa(a+2,b) = Fa(a+2,b) + etemp2 
                     Fa(b,a+2) = Fa(b,a+2) + etemp2

                     Fa(a+3,b) = Fa(a+3,b) + etemp3 
                     Fa(b,a+3) = Fa(b,a+3) + etemp3

                     Fa(a+4,b) = Fa(a+4,b) + etemp4 
                     Fa(b,a+4) = Fa(b,a+4) + etemp4
                  go to 33 
               endif  

               if (delta_a .gt. 4) then  
99                continue 

                     if (arange2-a .gt. 4) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 
                        etemp3 = dtemp*Xa(a+3,b) 
                        etemp4 = dtemp*Xa(a+4,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(b,a) = Fa(b,a) + etemp 

                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(b,a+1) = Fa(b,a+1) + etemp1

                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(b,a+2) = Fa(b,a+2) + etemp2

                        Fa(a+3,b) = Fa(a+3,b) + etemp3 
                        Fa(b,a+3) = Fa(b,a+3) + etemp3

                        Fa(a+4,b) = Fa(a+4,b) + etemp4 
                        Fa(b,a+4) = Fa(b,a+4) + etemp4
                        go to 22 
                     endif 

                     if (arange2-a .eq. 4) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 
                        etemp3 = dtemp*Xa(a+3,b) 
                        etemp4 = dtemp*Xa(a+4,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(b,a) = Fa(b,a) + etemp 

                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(b,a+1) = Fa(b,a+1) + etemp1

                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(b,a+2) = Fa(b,a+2) + etemp2

                        Fa(a+3,b) = Fa(a+3,b) + etemp3 
                        Fa(b,a+3) = Fa(b,a+3) + etemp3

                        Fa(a+4,b) = Fa(a+4,b) + etemp4 
                        Fa(b,a+4) = Fa(b,a+4) + etemp4
                        go to 33 
                     endif 
                     if (arange2-a .eq. 3) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 
                        etemp3 = dtemp*Xa(a+3,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(b,a) = Fa(b,a) + etemp 

                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(b,a+1) = Fa(b,a+1) + etemp1

                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(b,a+2) = Fa(b,a+2) + etemp2

                        Fa(a+3,b) = Fa(a+3,b) + etemp3 
                        Fa(b,a+3) = Fa(b,a+3) + etemp3
                        go to 33 
                     endif 

                     if (arange2-a .eq. 2) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 


                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(b,a) = Fa(b,a) + etemp 

                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(b,a+1) = Fa(b,a+1) + etemp1

                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(b,a+2) = Fa(b,a+2) + etemp2
                        go to 33 
                     endif 

                     if (arange2-a .eq. 1) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(b,a) = Fa(b,a) + etemp 

                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(b,a+1) = Fa(b,a+1) + etemp1
                        go to 33 
                     endif 

                     if (arange2-a .eq. 0) then 
                        etemp = dtemp*Xa(a,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(b,a) = Fa(b,a) + etemp 
                        go to 33 
                     endif 

22                continue   
                  a = a + 5 
                  go to 99 

               endif 
33             continue 

            enddo ! b 
         endif 
      enddo ! c 
      enddo ! d 

      return 
      end 


      subroutine fdmult1b(arange1,arange2,brange1,brange2, 
     *                    crange1,crange2,drange1,drange2,
     *                    a1,a2,b1,b2,c1,c2,d1,d2,intblk, 
     *                    factor,sym,nc1,nc2,nd1,nd2,
     *                    ca,fa,dcrit,itype) 
      implicit none
      include 'int_gen_parms.h'
      integer arange1,arange2,brange1,brange2,nc1,nc2,sym
      integer crange1,crange2,drange1,drange2,nd1,nd2
      integer a, b, c, d, p, delta_a, itype  
      integer a1,a2,b1,b2,c1,c2,d1,d2 
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2) 
      double precision ca(nc1:nc2,nd1:nd2)
      double precision fa(nc1:nc2,nc1:nc2)
      double precision Xa(nc1:nc2,nc1:nc2)
      double precision factor, dcrit, dtemp, etemp, etemp1, 
     *                 etemp2, etemp3, etemp4, etemp5  

      delta_a = arange2 - arange1 

      do d = drange1, drange2
      do c = crange1, crange2
         dtemp = ca(c,d)  
         dtemp = factor*dtemp  
         if (dabs(dtemp) .gt. dcrit) then 
            if (itype .eq. 1) then 
               do b = brange1, brange2
               do p = arange1, arange2, 5 
                  if (arange2-p .gt. 4) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                      Xa(p+3,b) = intblk(p+3,b,c,d) 
                      Xa(p+4,b) = intblk(p+4,b,c,d) 
                      go to 1 
                  endif 
                  if (arange2-p .eq. 4) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                      Xa(p+3,b) = intblk(p+3,b,c,d) 
                      Xa(p+4,b) = intblk(p+4,b,c,d) 
                  endif 
                  if (arange2-p .eq. 3) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                      Xa(p+3,b) = intblk(p+3,b,c,d) 
                  endif 
                  if (arange2-p .eq. 2) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                      Xa(p+2,b) = intblk(p+2,b,c,d) 
                  endif 
                  if (arange2-p .eq. 1) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                      Xa(p+1,b) = intblk(p+1,b,c,d) 
                  endif 
                  if (arange2-p .eq. 0) then  
                      Xa(p,  b) = intblk(p,b,c,d) 
                  endif 
1              continue 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 2) then 
               do b = brange1, brange2
               do p = arange1, arange2 
                  Xa(p,b) = intblk(c,d,p,b) 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 3) then 
               do b = brange1, brange2
               do p = arange1, arange2, 5  
                  if (arange2-p .ge. 4) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                     Xa(p+2,b) = intblk(p+2,c,b,d) 
                     Xa(p+3,b) = intblk(p+3,c,b,d) 
                     Xa(p+4,b) = intblk(p+4,c,b,d) 
                     go to 2 
                  endif 
                  if (arange2-p .eq. 3) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                     Xa(p+2,b) = intblk(p+2,c,b,d) 
                     Xa(p+3,b) = intblk(p+3,c,b,d) 
                  endif 
                  if (arange2-p .eq. 2) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                     Xa(p+2,b) = intblk(p+2,c,b,d) 
                  endif 
                  if (arange2-p .eq. 1) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                     Xa(p+1,b) = intblk(p+1,c,b,d) 
                  endif 
                  if (arange2-p .eq. 0) then  
                     Xa(p,  b) = intblk(p,c,b,d) 
                  endif 
2              continue 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 4) then 
               do b = brange1, brange2
               do p = arange1, arange2, 5  
                  if (arange2-p .ge. 4) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     Xa(p+2,b) = intblk(p+2,d,c,b) 
                     Xa(p+3,b) = intblk(p+3,d,c,b) 
                     Xa(p+4,b) = intblk(p+4,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 3) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     Xa(p+2,b) = intblk(p+2,d,c,b) 
                     Xa(p+3,b) = intblk(p+3,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 2) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     Xa(p+2,b) = intblk(p+2,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 1) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     Xa(p+1,b) = intblk(p+1,d,c,b) 
                     go to 3 
                  endif 
                  if (arange2-p .eq. 0) then  
                     Xa(p,b)   = intblk(p,d,c,b) 
                     go to 3 
                  endif 
3              continue 
               enddo  
               enddo  
            endif 
            if (itype .eq. 5) then 
               do b = brange1, brange2
               do p = arange1, arange2 
                  Xa(p,b) = intblk(c,b,p,d) 
               enddo  
               enddo  
            go to 44 
            endif 
            if (itype .eq. 6) then 
               do b = brange1, brange2
               do p = arange1, arange2 
                  Xa(p,b) = intblk(c,p,d,b) 
               enddo  
               enddo  
            go to 44 
            endif 

44          continue 

            do b = brange1, brange2

               a = arange1 

               if (delta_a .eq. 0) then  
                     etemp1  = dtemp*Xa(a,b) 
                     Fa(a,b) = Fa(a,b) + etemp1 
                  go to 33 
               endif  

               if (delta_a .eq. 1) then  
                     etemp   = dtemp*Xa(a,b) 
                     etemp1  = dtemp*Xa(a+1,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(a+1,b) = Fa(a+1,b) + etemp1 
                  go to 33 
               endif  

               if (delta_a .eq. 2) then  
                     etemp  = dtemp*Xa(a,b) 
                     etemp1 = dtemp*Xa(a+1,b) 
                     etemp2 = dtemp*Xa(a+2,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(a+1,b) = Fa(a+1,b) + etemp1  
                     Fa(a+2,b) = Fa(a+2,b) + etemp2 
                  go to 33 
               endif  

               if (delta_a .eq. 3) then  
                     etemp = dtemp*Xa(a,b) 
                     etemp1 = dtemp*Xa(a+1,b) 
                     etemp2 = dtemp*Xa(a+2,b) 
                     etemp3 = dtemp*Xa(a+3,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(a+1,b) = Fa(a+1,b) + etemp1 
                     Fa(a+2,b) = Fa(a+2,b) + etemp2 
                     Fa(a+3,b) = Fa(a+3,b) + etemp3 
                  go to 33 
               endif  

               if (delta_a .eq. 4) then  
                     etemp = dtemp*Xa(a,b) 
                     etemp1 = dtemp*Xa(a+1,b) 
                     etemp2 = dtemp*Xa(a+2,b) 
                     etemp3 = dtemp*Xa(a+3,b) 
                     etemp4 = dtemp*Xa(a+4,b) 

                     Fa(a,b) = Fa(a,b) + etemp  
                     Fa(a+1,b) = Fa(a+1,b) + etemp1 
                     Fa(a+2,b) = Fa(a+2,b) + etemp2 
                     Fa(a+3,b) = Fa(a+3,b) + etemp3 
                     Fa(a+4,b) = Fa(a+4,b) + etemp4 
                  go to 33 
               endif  

               if (delta_a .gt. 4) then  
99                continue 

                     if (arange2-a .gt. 4) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 
                        etemp3 = dtemp*Xa(a+3,b) 
                        etemp4 = dtemp*Xa(a+4,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(a+3,b) = Fa(a+3,b) + etemp3 
                        Fa(a+4,b) = Fa(a+4,b) + etemp4 
                        go to 22 
                     endif 

                     if (arange2-a .eq. 4) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 
                        etemp3 = dtemp*Xa(a+3,b) 
                        etemp4 = dtemp*Xa(a+4,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(a+3,b) = Fa(a+3,b) + etemp3 
                        Fa(a+4,b) = Fa(a+4,b) + etemp4 
                        go to 33 
                     endif 
                     if (arange2-a .eq. 3) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 
                        etemp3 = dtemp*Xa(a+3,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        Fa(a+3,b) = Fa(a+3,b) + etemp3 
                        go to 33 
                     endif 

                     if (arange2-a .eq. 2) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 
                        etemp2 = dtemp*Xa(a+2,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        Fa(a+2,b) = Fa(a+2,b) + etemp2 
                        go to 33 
                     endif 

                     if (arange2-a .eq. 1) then 
                        etemp  = dtemp*Xa(a,b) 
                        etemp1 = dtemp*Xa(a+1,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        Fa(a+1,b) = Fa(a+1,b) + etemp1 
                        go to 33 
                     endif 

                     if (arange2-a .eq. 0) then 
                        etemp = dtemp*Xa(a,b) 

                        Fa(a,b) = Fa(a,b) + etemp  
                        go to 33 
                     endif 

22                continue   
                  a = a + 5 
                  go to 99 

               endif 
33             continue 

            enddo ! b 
         endif 
      enddo ! c 
      enddo ! d 

      return 
      end 


