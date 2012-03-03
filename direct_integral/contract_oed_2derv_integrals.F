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
      subroutine contract_oed_2derv_integrals(type, density, z,
     *                       da1, da2, db1, db2, hess)
c---------------------------------------------------------------------------
c  Computes the set of OED second derivative integrals matching "type" and 
c  contracts them with density block (da1:da2, db1:db2). Note that the 
c  density block may be the density(for type = 1,2) or an intermediate
c  array(for type = 3).
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
      include 'parallel_info.h'
      include 'trace.h'

      integer type
      integer da1, da2, db1, db2
      double precision density(da1:da2, db1:db2)
      double precision z(da1:da2, db1:db2)
      double precision hess(3*ncenters, 3*ncenters)
      double precision dbuf(1)
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf, dbuf)
#endif

      integer i, j, m, n, nfirst, nints, nalpha, ncoeff, iatom, jatom 
      integer icomponent, jcomponent, ix, jx, ncsum  
      integer a, b, a1,a2,b1,b2,lastao, ihess, jhess 
      integer m1, m2, n1, n2
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer ixderc, jxderc 
      integer ind(2)
      integer flopcount
      logical spherical
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2, x, y, tess, fact

      if (type .lt. 1 .or. type .gt. 3) then
         print *,'Error: Invalid type in contract_oed_2derv_integrals'
         print *,'type = ',type
         call abort_job()
      endif

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif
      spherical = (ispherical .eq. 1)
      lastao    = end_nfps(nshells)

      ind(1) = 1
      ind(2) = 2

      if (intpkg .ne. flocke_package) then
         print *,'Error: contract_oed_2derv_integrals requires ERD ',
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
     *                       max_dim_coeff )

       if (type .eq. 2 .or. type .eq. 3) then 

c loop over derivatives on centers

         do icomponent = 1, 3   ! i component is x, y, or z term of hessian  
         do jcomponent = 1, 3   ! j component is x, y, or z term of hessian  
         do ix = 1, 2          ! index of desired center of i component 
         do jx = 1, 2          ! index of desired center of j component 

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
c
c Check if 'derivatives' equal 
c
             if((icomponent .eq. jcomponent) .and. (ix .eq. jx)) then 
c 
             if (icomponent .eq. 1) then
                if (ix .eq. 1) then
                   der1x = 2
                   if (atom(m) .eq. atom(n)) der2x = 2
                   iatom = atom(m)
                   jatom = atom(m)
                else
                   der2x = 2
                   if (atom(m) .eq. atom(n)) der1x = 2
                   iatom = atom(n)
                   jatom = atom(n)
                endif
             endif

             if (icomponent .eq. 2) then
                if (ix .eq. 1) then
                   der1y = 2
                   if (atom(m) .eq. atom(n)) der2y = 2
                   iatom = atom(m)
                   jatom = atom(m)
                else
                   der2y = 2
                   if (atom(m) .eq. atom(n)) der1y = 2
                   iatom = atom(n)
                   jatom = atom(n)
                endif
             endif

             if (icomponent .eq. 3) then
                if (ix .eq. 1) then
                   der1z = 2
                   if (atom(m) .eq. atom(n)) der2z = 2
                   iatom = atom(m)
                   jatom = atom(m)
                else
                   der2z = 2
                   if (atom(m) .eq. atom(n)) der1z = 2
                   iatom = atom(n)
                   jatom = atom(n)
                endif
             endif 
             go to 999 
             endif 
c
c set first derivative index 
c 
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
c
c set second derivative index 
c 
             if (jcomponent .eq. 1) then
                if (jx .eq. 1) then
                   der1x = der1x + 1
                   if (atom(m) .eq. atom(n)) der2x = der1x !1
                   jatom = atom(m)
                else
                   der2x = der2x + 1
                   if (atom(m) .eq. atom(n)) der1x = der2x !1
                   jatom = atom(n)
                endif
             endif

             if (jcomponent .eq. 2) then
                if (jx .eq. 1) then
                   der1y = der1y + 1
                   if (atom(m) .eq. atom(n)) der2y = der1y !1
                   jatom = atom(m)
                else
                   der2y = der2y + 1
                   if (atom(m) .eq. atom(n)) der1y = der2y !1
                   jatom = atom(n)
                endif
             endif

             if (jcomponent .eq. 3) then
                if (jx .eq. 1) then
                   der1z = der1z + 1
                   if (atom(m) .eq. atom(n)) der2z = der1z !1
                   jatom = atom(m)
                else
                   der2z = der2z + 1
                   if (atom(m) .eq. atom(n)) der1z = der2z !1
                   jatom = atom(n)
                endif
             endif

999         continue 
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
c   Sum into the hessian.
c---------------------------------------------------------------------------

              ihess = (iatom-1)*3 + icomponent 
              jhess = (jatom-1)*3 + jcomponent 

              hess(jhess,ihess) = hess(jhess,ihess) + y 

             endif

         enddo ! jx
         enddo ! ix
         enddo ! jcomponent
         enddo ! icomponent

         endif ! ktype = 2,3

         enddo   ! n shells
         enddo   ! m shells

      return
      end




