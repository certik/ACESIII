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
      subroutine compute_uintegrals8(a1,a2,b1,b2,c1,c2,d1,d2,scr,
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
     *         intpkg .ne. gamess_derivative_package .and.
     *         nints .gt. 0) then
 
              call uscf_tp8a(jx,nc1,nc2,nd1,nd2,ca,cb,a1,a2,b1,b2,
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
      subroutine uscf_tp8a(jx,nc1,nc2,nd1,nd2,ca,cb,va1,va2,vb1,vb2,
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
      integer dstart, dend  
      integer delta_a, ainc, amax, jx, itype  

      double precision ca(nc1:nc2,nd1:nd2)
      double precision cb(nc1:nc2,nd1:nd2)

      double precision Xa(nc1:nc2,nd1:nd2)
      double precision Xb(nc1:nc2,nd1:nd2)

      double precision Fa(nc1:nc2,nc1:nc2)
      double precision Fb(nc1:nc2,nc1:nc2)
      double precision ca_d(a1:a2)
      double precision cb_d(a1:a2)
      double precision ca_c(a1:a2)
      double precision cb_c(a1:a2) 
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2)
      double precision dtemp, etemp, itemp(5), dcrit, factor  
      double precision dtemp1, etemp1 
      double precision dtemp2, etemp2 
      double precision dtemp3, etemp3 
      double precision dtemp4, etemp4 
      double precision dtemp5, etemp5 
      double precision dtemp6, etemp6 
      double precision dtemp7, etemp7 
      double precision dtemp8, etemp8 
      double precision etempa_bd, etempa_bc, etempb_bd, etempb_bc 

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


      dstart = min(arange1,brange1,crange1,drange1)
      dend = max(arange2,brange2,crange2,drange2)

      do d = dstart, dend ! nd1, nd2 
      do c = dstart, dend ! nc1, nc2 

      if (jx .eq. 1) then

             dtemp = 0.0
             do p = 1, nalpha_occupied ! nd1,nd2
                dtemp = dtemp + ca(c,p)*ca(d,p)
             enddo
             xa(c,d) = dtemp

             dtemp = 0.0
             do p = 1, nbeta_occupied ! nd1,nd2
                dtemp = dtemp + cb(c,p)*cb(d,p)
             enddo
             xb(c,d) = dtemp

      else

             xa(c,d) = ca(c,d)
             xb(c,d) = cb(c,d)

      endif

      enddo
      enddo

      do d = drange1, drange2
         do a = arange1, arange2
            ca_d(a) = xa(a,d)
            cb_d(a) = xb(a,d)
         enddo 
      do c = crange1, crange2
         dtemp = 2.0*(xa(c,d) + xb(c,d))  
         do a = arange1, arange2
            ca_c(a) = xa(a,c)
            cb_c(a) = xb(a,c)
         enddo
            do b = brange1, brange2
               dtemp1 = -1.0*xa(b,d)   
               dtemp2 = -1.0*xa(b,c)   
               dtemp5 = -1.0*xb(b,d)   
               dtemp6 = -1.0*xb(b,c)   
               etempa_bd = 0.0 
               etempa_bc = 0.0 
               etempb_bd = 0.0 
               etempb_bc = 0.0 
            do a = arange1, arange2, 5 
                  amax = min(arange2-a,4)
                  ainc = 0

                                   itemp(1)  = intblk(a  ,b,c,d)
                  if (amax .ge. 1) itemp(2)  = intblk(a+1,b,c,d)
                  if (amax .ge. 2) itemp(3)  = intblk(a+2,b,c,d)
                  if (amax .ge. 3) itemp(4)  = intblk(a+3,b,c,d)
                  if (amax .ge. 4) itemp(5)  = intblk(a+4,b,c,d) 
33       continue

c              itemp = intblk(a+ainc,b,c,d) 

               etemp  = dtemp*itemp(ainc+1)  
               etemp1 = dtemp1*itemp(ainc+1) 
               etemp2 = dtemp2*itemp(ainc+1) 
               etemp5 = dtemp5*itemp(ainc+1) 
               etemp6 = dtemp6*itemp(ainc+1) 

               Fa(a+ainc,b) = Fa(a+ainc,b) + etemp 
               Fa(b,a+ainc) = Fa(b,a+ainc) + etemp 

               Fb(a+ainc,b) = Fb(a+ainc,b) + etemp 
               Fb(b,a+ainc) = Fb(b,a+ainc) + etemp 

               Fa(a+ainc,c) = Fa(a+ainc,c) + etemp1 

               Fa(a+ainc,d) = Fa(a+ainc,d) + etemp2 

               Fb(a+ainc,c) = Fb(a+ainc,c) + etemp5 

               Fb(a+ainc,d) = Fb(a+ainc,d) + etemp6 

               etempa_bc = etempa_bc + ca_d(a+ainc)*itemp(ainc+1)  
               etempb_bc = etempb_bc + cb_d(a+ainc)*itemp(ainc+1)  

               etempa_bd = etempa_bd + ca_c(a+ainc)*itemp(ainc+1)  
               etempb_bd = etempb_bd + cb_c(a+ainc)*itemp(ainc+1)  

                  ainc = ainc + 1
                  if (ainc .le. amax) go to 33

            enddo

               Fa(b,c) = Fa(b,c) - etempa_bc  
               Fb(b,c) = Fb(b,c) - etempb_bc 

               Fa(b,d) = Fa(b,d) - etempa_bd 
               Fb(b,d) = Fb(b,d) - etempb_bd 

            enddo

      enddo    
      enddo    


      return
      end
c
