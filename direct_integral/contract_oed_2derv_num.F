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
      subroutine contract_oed_2derv_numint(type, density, z,
     *                       da1, da2, db1, db2, hess)
c---------------------------------------------------------------------------
c  Computes the set of OED derivative integrals matching "type" and 
c  and contracts them with density block (da1:da2, db1:db2).
c
c   type = 
c          2 --> KIN integrals
c          3 --> OVL integrals
c---------------------------------------------------------------------------
      implicit none
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'hess.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      integer type
      integer da1, da2, db1, db2
      double precision density(da1:da2, db1:db2)
      double precision z(da1:da2, db1:db2)
      double precision hess(3*ncenters,3*ncenters)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif

      integer i, ix, m, n, nfirst, nints, nalpha, ncoeff, iatom,jatom
      integer component,icomponent,jcomponent,ihess,jhess,newatoms  
      integer a, b, a1,a2,b1,b2,lastao,ncsum,igrad
      integer m1, m2, n1, n2
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer ixderc
      integer ind(2)
      integer flopcount
      logical spherical
      logical*8 l8true, l8spherical
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *       pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2, x, y, fact
      double precision x1p,y1p,z1p,x2p,y2p,z2p, eps 

      if (type .lt. 2 .or. type .gt. 3) then
         print *,'Error: Invalid type in contract_oed_numint'
         print *,'type = ',type
         call abort_job()
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      spherical = (ispherical .eq. 1)
      l8spherical = spherical
      l8true = .true.
      lastao    = end_nfps(nshells)

      ind(1) = 1
      ind(2) = 2

      if (intpkg .ne. flocke_package) then
         print *,'Error: contract_oed_derv_integrals requires ERD ',
     *           'integral package.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Find the shells we need to loop over for this block of density data.
c--------------------------------------------------------------------------

      call lookup_shell(end_nfps, nshells, da1, m1)
      call lookup_shell(end_nfps, nshells, da2, m2)
      call lookup_shell(end_nfps, nshells, db1, n1)
      call lookup_shell(end_nfps, nshells, db2, n2)
      
c-------------------------------------------------------------------------
c   Build an array of coordinates per atom.
c-------------------------------------------------------------------------

      do m = 1, nshells
         iatom = atom(m)
         xn(iatom) = dbuf(icoord+3*(m-1))
         yn(iatom) = dbuf(icoord+3*(m-1)+1)
         zn(iatom) = dbuf(icoord+3*(m-1)+2)
      enddo

c-------------------------------------------------------------------------
c   Calculate the integrals over the necessary shell blocks.
c-------------------------------------------------------------------------

       do m = m1, m2
       do n = n1, n2

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

c loop over derivatives on centers

         do icomponent = 1, 3  ! component is x, y, or z term of gradient
         do ix = 1, 2          ! index of desired center

c loop over perturbations 

         do jatom = 1, ncenters 
         do jcomponent = 1, 3 

c-------------------------------------------------------------------------
c   Clear the integral block.
c-------------------------------------------------------------------------

          do b = db1, db2
          do a = da1, da2
             z(a,b) = 0.
          enddo
          enddo

             der1x = 0
             der1y = 0
             der1z = 0
             der2x = 0
             der2y = 0
             der2z = 0

             ixderc = 0
             dercx  = 0
             dercy  = 0
             dercz  = 0

             if (icomponent .eq. 1) then
                if (ix .eq. 1) then
                   der1x = 1
                   if (atom(m) .eq. atom(n)) der2x = 1
                   iatom = atom(m)
                else
                   der2x = 1
                   if (atom(m) .eq. atom(n)) der1x = 1
                   iatom = atom(n)
                endif
             endif

             if (icomponent .eq. 2) then
                if (ix .eq. 1) then
                   der1y = 1
                   if (atom(m) .eq. atom(n)) der2y = 1
                   iatom = atom(m)
                else
                   der2y = 1
                   if (atom(m) .eq. atom(n)) der1y = 1
                   iatom = atom(n)
                endif
             endif

             if (icomponent .eq. 3) then
                if (ix .eq. 1) then
                   der1z = 1
                   if (atom(m) .eq. atom(n)) der2z = 1
                   iatom = atom(m)
                else
                   der2z = 1
                   if (atom(m) .eq. atom(n)) der1z = 1
                   iatom = atom(n)
                endif
             endif

               eps   = 0.0000001d0
c
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

                if (atom(m) .eq. jatom) then

                   x1p = x1 + eps

                endif

                if (atom(n) .eq. jatom) then

                   x2p = x2 + eps

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

                if (atom(m) .eq. jatom) then

                   y1p = y1 + eps

                endif

                if (atom(n) .eq. jatom) then

                   y2p = y2 + eps

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

                if (atom(m) .eq. jatom) then

                   z1p = z1 + eps

                endif

                if (atom(n) .eq. jatom) then

                   z2p = z2 + eps

                endif

                go to 111

             endif
c
c-------------------------------------------------------------------------
c
111          continue
c
c Stationary integrals.
c --------------------- 
c
               if (type .eq. 2) then
                  call oed__gener_kin_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))
               else if (type .eq. 3) then
                  call oed__gener_ovl_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z, alpha_pack,
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

                call move_integrals2(z, da1,da2,db1,db2,
     *                 dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2)

             endif
c
c Perturbed integrals
c -------------------
c
               if (type .eq. 2) then
                  call oed__gener_kin_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1p,y1p,z1p,x2p,y2p,z2p,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z, alpha_pack,
     *                 pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .true., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))
               else if (type .eq. 3) then
                  call oed__gener_ovl_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m),ivangmom(n),x1p,y1p,z1p,x2p,y2p,z2p,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z, alpha_pack,
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

               fact = -1.0d0 
               call add_integrals2(z, da1,da2,db1,db2,
     *                 dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2,fact)

c--------------------------------------------------------------------------
c   Contract the density block with the block of integrals.
c--------------------------------------------------------------------------

               call cwork0222(y,
     *                        density,da1,da2,db1,db2,
     *                        da1,da2,db1,db2, ind,
     *                        z,da1,da2,db1,db2,
     *                        da1,da2,db1,db2,ind(1), ind(2),
     *                        flopcount)

c---------------------------------------------------------------------------
c   Sum into the gradient.
c---------------------------------------------------------------------------

                 ihess = (iatom-1)*3 + icomponent
                 jhess = (jatom-1)*3 + jcomponent

                 hess(jhess,ihess) = hess(jhess,ihess) + y/eps

             endif
c
         enddo ! jx
         enddo ! jcomponent
         enddo ! ix
         enddo ! icomponent

         enddo   ! n shells
         enddo   ! m shells

      return
      end
