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
      subroutine compute_oed_integrals(type, z)
c---------------------------------------------------------------------------
c  Computes the set of OED integrals matching "type" and returns them in 
c  array Z.
c
c   type = 1 --> NAI integrals
c          2 --> KIN integrals
c          3 --> OVL integrals
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'machine_types.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer type
      double precision z(*)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif
      integer i, m, n, nfirst, nints, nalpha, ncoeff, iatom, ncsum
       integer j
      integer a1,a2,b1,b2,lastao
      logical spherical
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2

      if (type .lt. 1 .or. type .gt. 3) then
         print *,'Error: Invalid type argument in compute_oed_integrals'
         print *,'type = ',type
         call abort_job()
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      spherical = (ispherical .eq. 1)
      lastao    = end_nfps(nshells)

      do i = 1, lastao*lastao
        z(i) = 0.
      enddo

c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = 1, nshells
         iatom = atom(m)
         xn(iatom) = dbuf(icoord+3*(m-1))
         yn(iatom) = dbuf(icoord+3*(m-1)+1)
         zn(iatom) = dbuf(icoord+3*(m-1)+2)
      enddo

c---------------------------------------------------------------------------
c   Calculate each shell block.
c---------------------------------------------------------------------------

      do m = 1, nshells
      do n = 1, nshells
         ncsum = ncfps(m) + ncfps(n) 
         x1 = dbuf(icoord+3*(m-1))
         y1 = dbuf(icoord+3*(m-1)+1)
         z1 = dbuf(icoord+3*(m-1)+2)
 
         x2 = dbuf(icoord+3*(n-1))
         y2 = dbuf(icoord+3*(n-1)+1)
         z2 = dbuf(icoord+3*(n-1)+2)

         call pack_coeffs_oed(dbuf(ialpha), ixalpha, dbuf(ipcoeff), 
     *                       ixpcoef, ncfps, npfps, m, n,
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ibuf(iccbeg), ibuf(iccend), 
     *                       indx_cc,
     *                       ccbeg_pack, ccend_pack,
     *                       max_dim_coeff)

            if (type .eq. 1) then
               call oed__gener_nai_batch(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 ncenters, xn,yn,zn, charge, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), nints, nfirst, 
     *                 dbuf(dpkgscr))     
            else if (type .eq. 2) then
               call oed__gener_kin_batch(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))
            else if (type .eq. 3) then
               call oed__gener_ovl_batch(intmax, zmax, nalpha, ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,spherical,
     *                 .true., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))
            endif

         if (nints .gt. 0) then
            a2 = end_nfps(m)
            if (m .eq. 1) then
               a1 = 1
            else
               a1 = end_nfps(m-1)+1
            endif

            b2 = end_nfps(n)
            if (n .eq. 1) then
               b1 = 1
            else
               b1 = end_nfps(n-1)+1
            endif

            call move_integrals2(z, 1,lastao,1,lastao,
     *             dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2)
         endif
      enddo
      enddo

      return
      end
