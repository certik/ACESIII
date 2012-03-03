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
      subroutine comp_return_derv_xyz(type,z, da1, da2, db1, db2, nsend)

c---------------------------------------------------------------------------
c  Computes the set of OED derivative moment integrals   
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'
      include 'dbugcom.h'
      include 'parallel_info.h'
#ifdef ALTIX
      include 'sheap.h'
#endif

      double precision dsum
      integer type, icount
      integer nsend 
      integer da1, da2, db1, db2
      integer adim, bdim 
      double precision z(da1:da2, db1:db2)

      double precision dbuf(1), temp
      integer ibuf(1)
#ifdef ALTIX
      pointer (iptr, ibuf)
      pointer (dptr, dbuf)
#else
      common dbuf
      equivalence (ibuf(1), dbuf(1))
#endif

      integer i, ix, m, n, nfirst, nints, nalpha, ncoeff, iatom, ncsum
      integer nalpha_pack, npcoeff_pack

      integer component
      integer a, b, a1,a2,b1,b2,lastao, igrad
      integer m1, m2, n1, n2
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer xmomder, ymomder, zmomder
      integer ind(2)
      integer ierr 

      integer flopcount
      logical spherical
   
      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      double precision xn(max_centers), yn(max_centers), 
     *                 zn(max_centers)
      double precision x1,y1,z1,x2,y2,z2, y, fact 

      common /d2int_com/jatom, jx, jcenter
      integer jatom, jx, jcenter, dcoord

      spherical = (ispherical .eq. 1)
      lastao = end_nfps(nshells)

#ifdef ALTIX
      iptr = ishptr
      dptr = dshptr
#endif

      if (dbg) 
     *  write(6,*) ' Entering comp_return_derv_xyz jatom jx ', 
     *               jatom,jx,jcenter 

c --------------------------------------------------------------------------- 
c  Set up send data. 
c --------------------------------------------------------------------------- 

      adim = da2-da1+1
      bdim = db2-db1+1 

      nsend = adim*bdim
      if (nsend .lt. 0) then
         print *,'ERROR IN COMP_RETURN_DERV_XYZ ',me,' nsend = ',nsend
         print *,'adim,bdim = ',adim,bdim
         call abort_job()
      endif 

c----------------------------------------------------------------------------
c   Clear the output array.
c----------------------------------------------------------------------------

      do b = db1,db2
      do a = da1,da2
         Z(a,b) = 0.D0
      enddo
      enddo

      ind(1) = 1
      ind(2) = 2

      if (intpkg .ne. flocke_package) then
         print *,'Error: comp_return_derv_xyz requires ERD ',
     *           'integral package.'
         call abort_job()
      endif

c-------------------------------------------------------------------------
c   Find the shells we need to loop over for this block.
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

c         if ((atom(m) .eq. jatom) .and. (atom(n) .eq. jatom)) 
c    *       go to 333  

c         if ((atom(m) .ne. jatom) .and. (atom(n) .ne. jatom)) 
c    *       go to 333  

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

             do ix = 1, 2  ! Centers

             der1x = 0
             der1y = 0
             der1z = 0
             der2x = 0
             der2y = 0
             der2z = 0

c            if (atom(m) .eq. jatom) then 
c               if (jx .eq. 1) der1x = 1 
c               if (jx .eq. 2) der1y = 1 
c               if (jx .eq. 3) der1z = 1 
c            endif 

c            if (atom(n) .eq. jatom) then 
c               if (jx .eq. 1) der2x = 1 
c               if (jx .eq. 2) der2y = 1 
c               if (jx .eq. 3) der2z = 1 
c            endif 

             if ((ix .eq. 1) .and. (atom(m) .ne. jatom)) go to 777 
             if ((ix .eq. 2) .and. (atom(n) .ne. jatom)) go to 777 

             if (jx .eq. 1) then
                if (ix .eq. 1) then
                   der1x = 1
C                   if (atom(m) .eq. atom(n)) der2x = 1
                else
                   der2x = 1
C                   if (atom(m) .eq. atom(n)) der1x = 1
                endif
             endif

             if (jx .eq. 2) then
                if (ix .eq. 1) then
                   der1y = 1
C                   if (atom(m) .eq. atom(n)) der2y = 1
                else
                   der2y = 1
C                   if (atom(m) .eq. atom(n)) der1y = 1
                endif
             endif

             if (jx .eq. 3) then
                if (ix .eq. 1) then
                   der1z = 1
C                   if (atom(m) .eq. atom(n)) der2z = 1
                else
                   der2z = 1
C                   if (atom(m) .eq. atom(n)) der1z = 1
                endif
             endif

             xmomder = 0
             ymomder = 0
             zmomder = 0
             if ( type .eq. 1 ) xmomder = 1
             if ( type .eq. 2 ) ymomder = 1
             if ( type .eq. 3 ) zmomder = 1

             call oed__gener_xyz_derv_batch(intmax,zmax,
     *                 nalpha,ncoeff,
     *                 ncsum, ncfps(m), ncfps(n), npfps(m),npfps(n),
     *                 ivangmom(m), ivangmom(n), x1,y1,z1,x2,y2,z2,
     *                 der1x, der1y, der1z,
     *                 der2x, der2y, der2z,
     *                 xmomder, ymomder, zmomder,
     *                 alpha_pack, pcoeff_pack, ccbeg_pack, ccend_pack,
     *                 spherical, .false., ibuf(ipkgscr), nints, nfirst,
     *                 dbuf(dpkgscr))


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

               fact = 1.0d0 

               call add_integrals2(z, da1,da2,db1,db2,
     *                   dbuf(dpkgscr+nfirst-1),a1,a2,b1,b2,fact)

             endif ! (nints .gt. 0) then

777          continue 

         enddo ! ix

333      continue 
         enddo   ! n shells
         enddo   ! m shells

7000  continue 

      return
      end
