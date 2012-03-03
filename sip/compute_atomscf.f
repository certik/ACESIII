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
      subroutine compute_atomscf(watom, scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 nai, kin, ovl,  
     *                 fa,fb) 
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s
      integer a,b,c,d
      integer iatom, n_basis  
      double precision watom 

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
      double precision nai(nc1:nc2,nd1:nd2)
      double precision kin(nc1:nc2,nd1:nd2)
      double precision ovl(nc1:nc2,nd1:nd2)
      double precision H0T(nc1:nc2,nd1:nd2)

      double precision fa(nc1:nc2,nc1:nc2)
      double precision fb(nc1:nc2,nc1:nc2)
      integer map(nc1:nc2) 
      integer umap(nc1:nc2) 
      integer beg_anfps(max_shells)  
      integer end_anfps(max_shells)  
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

c     call mpi_comm_rank(mpi_comm_world, me, ierr)

      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

      iatom = watom 
c      write(6,*) ' Performing an SCF calculation on atom:', iatom 

      call comp_return_h0(H0T, iatom, nc1, nc2, nc1, nc2) 

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

         m1 = 1 
         n1 = 1 
         r1 = 1 
         s1 = 1 

         m2 = (nshells)   
         n2 = (nshells)  
         r2 = (nshells) 
         s2 = (nshells)  

c-----------------------------------------------------------------------
c   Find the number of basis functions and shells in the atom.  
c-----------------------------------------------------------------------

         n_basis = 0 
         do m = m1, m2 
            if(atom(m) .eq. iatom) n_basis = n_basis + 
     *                end_nfps(m) - end_nfps(m-1)   
         enddo 
c         write(6,*) ' The number of basis functions on atom', 
c     *                iatom, '=', n_basis  

c-----------------------------------------------------------------------
c   Find the mapping from atom <--> molecule.  
c-----------------------------------------------------------------------

         do n = nc1, nc2 
            map(n) = 0 
            umap(n) = 0 
         enddo 

         n_basis = 0 
         do m = m1, m2 
            beg_anfps(m) = 0   
            end_anfps(m) = 0   

            if(atom(m) .eq. iatom) then 
               beg_anfps(m) = n_basis + 1  

               if (m .eq. 1) then
                  DO n = 1, end_nfps(m) 
                     n_basis = n_basis + 1 
                     map(n_basis) = n 
                     umap(n) = n_basis  
                  enddo 
               else 
                  DO n = end_nfps(m-1) + 1, end_nfps(m) 
                     n_basis = n_basis + 1 
                     map(n_basis) = n 
                     umap(n) = n_basis  
                  enddo 
               endif 

               end_anfps(m) = n_basis 

            endif 

         enddo 

c        do n = nc1, nc2 
c           write(6,*) 'n umap(n)', n, umap(n)  
c        enddo 
c        do m = m1, m2 
c           write(6,*) ' Mth shell:', m, end_anfps(m) 
c        enddo 
c         write(6,*) ' The number of basis functions on atom', 
c     *                iatom, '=', n_basis  

         call do_atomscf(watom, scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 H0T, nai, kin, ovl,  
     *                 fa,  fb, 
     *                 n_basis, beg_anfps, end_anfps, map, umap) 

      return 
      end 

      subroutine do_atomscf(watom, scr,
     *                 maxblk, iscr, coords,coeffs,alphas, ccbeg, ccend,
     *                 nc1,nc2, nd1, nd2,
     *                 HOT, nai, kin, ovl,  
     *                 fina,  finb, 
     *                 n_basis, beg_anfps, end_anfps, map, umap) 
c---------------------------------------------------------------------------

      implicit none

      include 'mpif.h'
      include 'int_gen_parms.h'
      include 'machine_types.h'

      integer a1, a2, b1, b2, c1, c2, d1, d2 
      integer aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2
      integer adim, bdim, cdim, ddim  
      integer m1, m2, n1, n2, r1, r2, s1, s2
      integer i, j, n, m, r, s, l, mn, rs  
      integer a,b,c,d
      integer iatom, n_basis  
      double precision watom 

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
      double precision nai(nc1:nc2,nd1:nd2)
      double precision kin(nc1:nc2,nd1:nd2)
      double precision ovl(nc1:nc2,nd1:nd2)
      double precision HOT(nc1:nc2,nd1:nd2)

      double precision fina(nc1:nc2,nc1:nc2)
      double precision finb(nc1:nc2,nc1:nc2)

      double precision h0(n_basis,n_basis) 
      double precision aovl(n_basis,n_basis) 
      double precision sos(n_basis,n_basis) 
      double precision Qxx(n_basis,n_basis) 

      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision HFDOLD_A(n_basis,n_basis) 
      double precision HFDOLD_B(n_basis,n_basis) 

      double precision ca(n_basis,n_basis) 
      double precision cb(n_basis,n_basis) 
      double precision cba(n_basis,n_basis) 
      double precision cbb(n_basis,n_basis) 
      double precision FTa(n_basis,n_basis) 
      double precision FTb(n_basis,n_basis) 
      double precision Fa(n_basis,n_basis) 
      double precision Fb(n_basis,n_basis) 
      double precision temp, tempa, tempb   
      integer doit, itemp, p, p1  

      integer nocc_a, nocc_b 
      integer iter, max_iter 

      integer map(nc1:nc2) 
      integer umap(nc1:nc2) 
      integer beg_anfps(max_shells)  
      integer end_anfps(max_shells)  
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

c     call mpi_comm_rank(mpi_comm_world, me, ierr)

      l8true = .true.
      spherical = (ispherical .eq. 1)
      l8spherical = spherical

      iatom = watom 

      do m = 1, max_centers 
         if (m .eq. iatom) then 
            nocc_b = charge(m)/2 
            nocc_a = charge(m) - nocc_b 
         endif 
      enddo 

c      write(6,*) ' Performing an SCF calculation on atom:', iatom, 
c     * 'in a basis of', n_basis, 'functions with', nocc_a, nocc_b, 
c     * 'alpha and beta occupied electrons'   

c-----------------------------------------------------------------------
c   Find the shell blocks for which we shall loop through.
c-----------------------------------------------------------------------

      m1 = 1 
      n1 = 1 
      r1 = 1 
      s1 = 1 

      m2 = (nshells)   
      n2 = (nshells)  
      r2 = (nshells) 
      s2 = (nshells)  

c-----------------------------------------------------------------------
c Sum nai and kin into small array and copy ovl there too. 
c --> initial guess   
c-----------------------------------------------------------------------

      itemp = 0 
      do n = nc1, nc2 
      do m = nc1, nc2 
         if (umap(m).ne.0 .and. umap(n).ne.0) then 
            aovl(umap(m),umap(n)) = ovl(m,n) 
c           h0(umap(m),umap(n))   = nai(m,n) + kin(m,n)  
            h0(umap(m),umap(n))   = HOT(m,n)   
            itemp = itemp + 1 
         endif 
      enddo  
      enddo  
      
      if (itemp .ne. n_basis**2) then 
         write(6,*) ' Something wrong with umap ', itemp, n_basis 
         call abort_job()
      endif 

c-----------------------------------------------------------------------
c Construct the hcore initial guess  
c-----------------------------------------------------------------------

      do n = 1, n_basis  
      do m = 1, n_basis  
         FA(m,n) = h0(m,n) 
         FB(m,n) = h0(m,n) 
      enddo  
      enddo  

c-----------------------------------------------------------------------
c Construct U*S**(-1/2)  
c-----------------------------------------------------------------------

      call diag(aovl,sos,m,n_basis,0,1,1) 
      do m = 1, n_basis  
      do n = 1, n_basis  
         temp = 0.0 
         do l = 1, n_basis  
            temp = temp + sos(m,l)*aovl(l,n) 
         enddo 
         Qxx(m,n) = temp 
       enddo  
       enddo  

c-----------------------------------------------------------------------
c Transpose the Fock matrix -> Construct S^(-1/2) F S^(-1/2)  
c-----------------------------------------------------------------------

       call fock_transpose(FA,FB,Qxx,FTa,FTb,n_basis) 

c-----------------------------------------------------------------------
c Diagonalize the transposed Fock matrix  
c-----------------------------------------------------------------------

       call diag(FTa,ca,m,n_basis,0,0,0) 
       call diag(FTb,cb,m,n_basis,0,0,0) 

c-----------------------------------------------------------------------
c Back transform the coefficient array  
c-----------------------------------------------------------------------

       call c_backtran(Qxx,ca,cb,cba,cbb,n_basis) 

c-----------------------------------------------------------------------
c Compute the HF density  
c-----------------------------------------------------------------------

       call hfdensity(ca,cb,HFD_A,HFD_B,n_basis,nocc_a,nocc_b) 

c-----------------------------------------------------------------------
c Compute the HF energy   
c-----------------------------------------------------------------------

       call hfenergy(HFD_A,HFD_B,FA,FB,h0,n_basis) 

c-----------------------------------------------------------------------
c Copy the HF Density into the old Density. 
c-----------------------------------------------------------------------

      call hfdensity_copy(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis)  

c-----------------------------------------------------------------------
c Start the SCF iterations  
c-----------------------------------------------------------------------

      max_iter = 30  
      DO iter = 1, max_iter 

c-----------------------------------------------------------------------
c       Construct the new Fock matrix  
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       One-electron piece  
c-----------------------------------------------------------------------

        do m = 1, n_basis 
        do n = 1, n_basis 
           FA(m,n) = h0(m,n)  
           FB(m,n) = h0(m,n)  
        enddo 
        enddo 

c-----------------------------------------------------------------------
c       Two-electron piece  
c-----------------------------------------------------------------------

         do m = m1, m2
            if(atom(m) .eq. iatom) then  
            aa1 = beg_anfps(m)
            aa2 = end_anfps(m)

            x1 = coords(1,m)
            y1 = coords(2,m)
            z1 = coords(3,m)
         do n = n1, n2
            if (m .le. n) then 
            if(atom(n) .eq. iatom) then  
            bb1 = beg_anfps(n)
            bb2 = end_anfps(n)

            x2 = coords(1,n)
            y2 = coords(2,n)
            z2 = coords(3,n)
         do r = r1, r2
            if(atom(r) .eq. iatom) then  
            cc1 = beg_anfps(r)
            cc2 = end_anfps(r)

            x3 = coords(1,r)
            y3 = coords(2,r)
            z3 = coords(3,r)
         do s = s1, s2
            if (r .le. s) then 
            if(atom(s) .eq. iatom) then  
            dd1 = beg_anfps(s)
            dd2 = end_anfps(s)

            doit = 0 

         if (( r .lt. s) .and.  
     *       ( n .ne. s) .and.  
     *       ( n .ne. r) .and.  
     *       ( m .lt. n) .and.  
     *       ( m .lt. r) .and.  
     *       ( m .ne. s)) doit = 1  

         if ((r .lt. s) .and.  
     *       (m .eq. n)) doit = 1  

         if (( r .lt. s) .and.   
     *    ( n .lt. s) .and.  
     *    ( m .lt. n) .and.  
     *    ( m .eq. r)) doit = 1   

         if (( r .lt. s) .and.  
     *    ( n .eq. r) .and.  
     *    ( m .lt.  n) .and.  
     *    ( m .lt. s)) doit = 1  

         if (( r .lt. s) .and.  
     *    ( n .eq. s) .and.  
     *    ( m .lt. n) .and.  
     *    ( m .lt. r)) doit = 1  
c
         if (( n .eq. m ) .and.  
     *    ( s .eq. m ) .and.  
     *    ( r .eq. m )) doit = 1  
        
         if (( r .eq. s) .and.  
     *    ( m .eq. n) .and.  
     *    ( m .lt. r)) doit = 1  
 
         if (( r .lt. s ) .and.  
     *    ( n .eq. s) .and.  
     *    ( m .lt. n) .and.  
     *    ( m .eq. r)) doit = 1  
 
c
c-----------------------------------------------------------------------
c   Determine the largest density element.
c-----------------------------------------------------------------------

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

c              if (doit .eq. 1) then 

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

c               endif 

c---------------------------------------------------------------------------
c   Move the integrals into the output block.  
c---------------------------------------------------------------------------

           if (nints .gt. 0) then

               call form_ss1fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss2fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss3fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss4fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss5fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss6fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss7fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

               call form_ss8fock(scr(nfirst), n_basis,
     *                        aa1,aa2,bb1,bb2,cc1,cc2,dd1,dd2,
     *                        HFD_A,HFD_B,FA,FB)

           endif
c 
            endif 
            endif 
         enddo   ! s
            endif 
         enddo   ! r
            endif 
            endif 
         enddo   ! n
            endif 
         enddo   ! m

c-----------------------------------------------------------------------
c       Done computing the Fock matrix   
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       Compute the new HF energy  
c-----------------------------------------------------------------------

        call hfenergy(HFD_A,HFD_B,FA,FB,h0,n_basis) 

c-----------------------------------------------------------------------
c       Transpose the new Fock Matrix   
c-----------------------------------------------------------------------

        call fock_transpose(Fa,FB,Qxx,FTa,FTb,n_basis) 

c-----------------------------------------------------------------------
c       Diagonalize the new Transposed Fock Matrix   
c-----------------------------------------------------------------------

        call diag(FTa,ca,m,n_basis,0,0,0) 
        call diag(FTb,cb,m,n_basis,0,0,0) 

c-----------------------------------------------------------------------
c       Back Transform the coefficient array  
c-----------------------------------------------------------------------

        call c_backtran(Qxx,ca,cb,cba,cbb,n_basis) 

c-----------------------------------------------------------------------
c       Check on convergence and replace the old density with the new   
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       Compute the new HF density  
c-----------------------------------------------------------------------

        call hfdensity(ca,cb,HFD_A,HFD_B,n_basis,nocc_a,nocc_b) 

c-----------------------------------------------------------------------
c       Check for convergence    
c-----------------------------------------------------------------------

        call check_conv(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis,doit)  
        if (doit .eq. 1) go to 100 

c-----------------------------------------------------------------------
c Copy the HF Density into the old Density. 
c-----------------------------------------------------------------------

        call hfdensity_copy(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis)  


      ENDDO ! iter = 1, max_iter 
100   continue 

c      write(6,*) ' Alpha Orbital energies ' 
c      write(6,*) ' -----------------------------------------------' 
c      do m = 1, nocc_a 
c         write(6,*) '  ', m, FTA(m,m), FTA(m,m)*27.21138386 
c      enddo  
c      write(6,*) ' -----------------------------------------------' 
c      do m = nocc_a + 1, n_basis  
c         write(6,*) '  ', m, FTA(m,m), FTA(m,m)*27.21138386  
c      enddo  
c      write(6,*) ' -----------------------------------------------' 
c      write(6,*) ' ' 

c      write(6,*) ' Beta Orbital energies ' 
c      write(6,*) ' -----------------------------------------------' 
c      do m = 1, nocc_b 
c         write(6,*) '  ', m, FTB(m,m), FTB(m,m)*27.21138386 
c      enddo  
c      write(6,*) ' -----------------------------------------------' 
c      do m = nocc_b + 1, n_basis  
c         write(6,*) '  ', m, FTB(m,m), FTB(m,m)*27.21138386  
c      enddo  
c      write(6,*) ' -----------------------------------------------' 
c      write(6,*) ' ' 

c         write(6,*) ' Done SCF calculation of ATOM :', iatom 

c-----------------------------------------------------------------------
c     Put the density matrix into the full matrix   
c-----------------------------------------------------------------------

      do n = 1, nbasis
      do m = 1, nbasis
         if ((map(m) .ne. 0) .and. (map(n) .ne. 0)) then 
         Fina(map(m),map(n)) = HFD_a(m,n)  
         Finb(map(m),map(n)) = HFD_b(m,n)   
         endif 
      enddo 
      enddo 

      return
      end
 

      subroutine form_fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                     dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s 
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      do s = dd1, dd2 
      do r = cc1, cc2 
         temp_rs = HFD_A(r,s) + HFD_B(r,s) 
      do n = bb1, bb2 
      do m = aa1, aa2 
           temp_mn  = int_block(m,n,r,s)*temp_rs 
           F_A(m,n) = F_A(m,n) + temp_mn 
           F_B(m,n) = F_B(m,n) + temp_mn 
      enddo 
      enddo 
      enddo 
      enddo 

      do s = dd1, dd2 
      do r = cc1, cc2 
      do n = bb1, bb2 
         tempa = HFD_A(n,s) 
         tempb = HFD_B(n,s) 
      do m = aa1, aa2 
           F_A(m,r) = F_A(m,r) - int_block(m,n,r,s)*tempa ! HFD_A(n,s)  
           F_B(m,r) = F_B(m,r) - int_block(m,n,r,s)*tempb ! HFD_B(n,s)  
      enddo 
      enddo 
      enddo 
      enddo 

      return 
      end 

      subroutine form_sfock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s 
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

c-----------------------------------------------------------------------
c m=n=r=s  
c-----------------------------------------------------------------------

c     if (m.eq.n .and. r.eq.s .and. n.eq.r) then 

         do s = dd1, dd2 
         do r = cc1, cc2 
            temp_rs = HFD_A(r,s) + HFD_B(r,s) 
         do n = bb1, bb2 
         do m = aa1, aa2 
              temp_mn  = int_block(m,n,r,s)*temp_rs 
              F_A(m,n) = F_A(m,n) + temp_mn 
              F_B(m,n) = F_B(m,n) + temp_mn 
         enddo 
         enddo 
         enddo 
         enddo 

         do s = dd1, dd2 
         do r = cc1, cc2 
            do m = aa1, aa2 
               F1A(m) = 0.0 
               F1B(m) = 0.0 
            enddo 
         do n = bb1, bb2 
            tempa = HFD_A(n,s) 
            tempb = HFD_B(n,s) 
         do m = aa1, aa2 
              F1A(m) = F1A(m) - int_block(m,n,r,s)*tempa ! HFD_A(n,s)  
              F1B(m) = F1B(m) - int_block(m,n,r,s)*tempb ! HFD_B(n,s)  
         enddo 
         enddo 
            do m = aa1, aa2 
               F_A(m,r) = F_A(m,r) + F1A(m)  
               F_B(m,r) = F_A(m,r) + F1B(m)  
            enddo 
         enddo 
         enddo 

c     endif 

c-----------------------------------------------------------------------
c  m.lt.n .and. r.lt.s .and. m.lt.r .and. n.ne.s .and. n.ne.r .and. m.ne.s   
c-----------------------------------------------------------------------

      if(m.lt.n .and. r.lt.s .and. m.lt.r .and. n.ne.s .and. 
     *   n.ne.r .and. m.ne.s) then  

      endif 

      return 
      end 

      subroutine fock_transpose(Fa,FB,Qxx,FTa,FTb,n_basis) 
      implicit none 
      integer n_basis, m, n, l, s  
      double precision Fa(n_basis,n_basis)  
      double precision Fb(n_basis,n_basis)  
      double precision FTa(n_basis,n_basis)  
      double precision FTb(n_basis,n_basis)  
      double precision Qxx(n_basis,n_basis)  
      double precision tempa, tempb 

      do m = 1, n_basis 
      do n = 1, n_basis 
         FTa(m,n) = 0.0  
         FTb(m,n) = 0.0  
      enddo 
      enddo 

      do m = 1, n_basis  
      do s = 1, n_basis  
         tempa = 0.0 
         tempb = 0.0 
         do l = 1, n_basis  
            tempa = tempa + Qxx(l,m)*FA(l,s) 
            tempb = tempb + Qxx(l,m)*FB(l,s) 
         enddo 

         do n = 1, n_basis  
            FTa(m,n) = FTa(m,n) + tempa*Qxx(s,n) 
            FTb(m,n) = FTb(m,n) + tempb*Qxx(s,n) 
         enddo 
      enddo 
      enddo 

      return 
      end 

      subroutine c_backtran(Qxx,ca,cb,cba,cbb,n_basis) 
      implicit none 
      integer n_basis 
      integer m, p, n  
      double precision Qxx(n_basis,n_basis)  
      double precision ca(n_basis,n_basis)  
      double precision cb(n_basis,n_basis)  
      double precision cba(n_basis,n_basis)  
      double precision cbb(n_basis,n_basis) 
      double precision tempa, tempb  

      do p = 1, n_basis  
      do m = 1, n_basis  
         tempa = 0.0  
         tempb = 0.0  
         do n = 1, n_basis  
            tempa = tempa + Qxx(m,n)*ca(n,p) 
            tempb = tempb + Qxx(m,n)*cb(n,p) 
         enddo 
         cba(m,p) = tempa
         cbb(m,p) = tempb
      enddo 
      enddo 

      do p = 1, n_basis  
      do m = 1, n_basis  
         ca(m,p) = cba(m,p) 
         cb(m,p) = cbb(m,p) 
      enddo 
      enddo 

      return 
      end 


      subroutine hfdensity(ca,cb,HFD_A,HFD_B,n_basis,nocc_a,nocc_b) 
      implicit none 
      integer n_basis, nocc_a, nocc_b 
      integer m, n, i, j 
      double precision ca(n_basis,n_basis) 
      double precision cb(n_basis,n_basis) 
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision tempa, tempb  

      do m = 1, n_basis  
      do n = 1, n_basis  
         tempa = 0.0 
         tempb = 0.0 
         do i = 1, nocc_a 
            tempa = tempa + ca(m,i)*ca(n,i) 
         enddo 
         do j = 1, nocc_b 
            tempb = tempb + cb(m,j)*cb(n,j) 
         enddo 
         HFD_A(m,n) = tempa 
         HFD_B(m,n) = tempb 
      enddo 
      enddo 

      return 
      end 

      subroutine hfenergy(HFD_A,HFD_B,FA,FB,h0,n_basis) 
      implicit none 
      integer n_basis 
      integer m, n 
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision FA(n_basis,n_basis) 
      double precision FB(n_basis,n_basis) 
      double precision h0(n_basis,n_basis) 
      double precision ea, eb, etotal 

      ea = 0.0 
      eb = 0.0 
      etotal = 0.0 

      do n = 1, n_basis  
      do m = 1, n_basis  
         ea = ea + (h0(m,n)+FA(m,n))*HFD_A(m,n)   
         eb = eb + (h0(m,n)+FB(m,n))*HFD_B(m,n)   
      enddo  
      enddo  

      etotal = 0.5d0*(ea + eb) 
c      write(6,*) ' Total SCF energy(-NN) = ', etotal 

      return 
      end 
   
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss1fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
      double precision t8xx
      double precision t9xx
      double precision t10xx
      double precision t11xx
      double precision t12xx
      double precision t13xx
      double precision t14xx
      double precision t15xx
      double precision t16xx
      double precision t17xx
      double precision t18xx
      double precision t19xx
      double precision t20xx
      double precision t21xx
      double precision t22xx
c
c
c############ CLASS 6  ############
c SS1 
c#  Do four-center part (mu nu |la si)->(mu nu |si la ),(nu mu|la si),(nu mu|si la)
c#                      (la si |mu nu),(si la |mu nu),(la si |nu mu),(si la |nu mu)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if ( l .lt. s) then 
            T1xx = HFD_A(l,s) + HFD_B(l,s)
            T5xx = 0.0 
      do n = bb1, bb2 
         if ( n .ne. s) then 
         if ( n .ne. l) then 
            T11xx = 0.0 
            T13xx = 0.0 
            T19xx = 0.0 
            T21xx = 0.0 
      do m = aa1, aa2 
c
         if ( m .lt. n) then 
         if ( m .lt. l) then 
         if ( m .ne. s) then 

            wint      = int_block(m,n,l,s) 
            T4xx      = HFD_A(m,n) + HFD_B(m,n)
            T2xx      = 2.0*wint*T1xx
            T5xx      = T5xx + 2.0*wint*T4xx
            T7xx      = wint*HFD_A(n,s)
            T9xx      = wint*HFD_A(n,l)
            T11xx     = T11xx + wint*HFD_A(m,s)
            T13xx     = T13xx + wint*HFD_A(m,l)
            T19xx     = T19xx + wint*HFD_B(m,s)
            T21xx     = T21xx + wint*HFD_b(m,l)

            F_a(m,n)  = F_a(m,n) + T2xx
            F_b(m,n)  = F_b(m,n) + T2xx
            F_a(n,m)  = F_a(n,m) + T2xx
            F_b(n,m)  = F_b(n,m) + T2xx

            F_A(m,l) = F_a(m,l) - T7xx
            F_A(l,m) = F_a(l,m) - T7xx

            F_a(m,s) = F_a(m,s) - T9xx
            F_A(s,m) = F_a(s,m) - T9xx

            T15xx    = wint*HFD_B(n,s)
            F_B(m,l) = F_b(m,l) - T15xx
            F_B(l,m) = F_b(l,m) - T15xx

            T17xx    = wint*HFD_B(n,l)
            F_b(m,s) = F_b(m,s) - T17xx
            F_b(s,m) = F_b(s,m) - T17xx

         endif 
         endif 
         endif 

      enddo ! m  
            F_A(n,l) = F_a(n,l) - T11xx
            F_A(l,n) = F_a(l,n) - T11xx
            F_A(n,s) = F_a(n,s) - T13xx
            F_A(s,n) = F_a(s,n) - T13xx
            F_b(n,l) = F_b(n,l) - T19xx
            F_b(l,n) = F_b(l,n) - T19xx
            F_b(n,s) = F_b(n,s) - T21xx
            F_b(s,n) = F_b(s,n) - T21xx
         endif 
         endif 
      enddo ! n  
            F_a(l,s) = F_a(l,s) + T5xx
            F_b(l,s) = F_b(l,s) + T5xx
            F_a(s,l) = F_a(s,l) + T5xx
            F_b(s,l) = F_b(s,l) + T5xx
         endif 
      enddo ! l  
      enddo ! s  
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss2fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
      double precision t8xx
      double precision t9xx
      double precision t10xx
      double precision t11xx
      double precision t12xx
      double precision t13xx
      double precision t14xx
      double precision t15xx
      double precision t16xx
      double precision t17xx
      double precision t18xx
      double precision t19xx
      double precision t20xx
      double precision t21xx
      double precision t22xx
c
c########### CLASS 3  ############
c SS2 
c  Do two-center part (m m |m n)->(m n |m m ),(m m|n m),(n m|m m)
c  Do three-center part (m m |n l)->(m m |l n ),(n l|m m),(l n|m m)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if (l .lt. s) then 
             T1xx = HFD_a(l,s) + HFD_b(l,s)
             T4xx = 0.0  
      do n = bb1, bb2 
      do m = aa1, aa2 
         if (m .eq. n) then 

             wint     = int_block(m,n,l,s) 
             T3xx     = HFD_a(m,n) + HFD_b(m,n)
             T2xx     = 2.0*wint*T1xx
             T4xx     = T4xx + wint*T3xx
             T6xx     = wint*HFD_a(n,s)
             T8xx     = wint*HFD_a(n,l)
             T10xx    = wint*HFD_b(n,s)
             T12xx    = wint*HFD_b(n,l)

             F_a(m,n) = F_a(m,n) + T2xx
             F_b(m,n) = F_b(m,n) + T2xx

             F_a(m,l) = F_a(m,l) - T6xx
             F_a(l,m) = F_a(l,m) - T6xx

             F_a(m,s) = F_a(m,s) - T8xx
             F_a(s,m) = F_a(s,m) - T8xx

             F_b(m,l) = F_b(m,l) - T10xx
             F_b(l,m) = F_b(l,m) - T10xx

             F_b(m,s) = F_b(m,s) - T12xx
             F_b(s,m) = F_b(s,m) - T12xx

         endif 
      enddo 
      enddo 
             F_a(l,s) = F_a(l,s) + T4xx
             F_b(l,s) = F_b(l,s) + T4xx
             F_a(s,l) = F_a(s,l) + T4xx
             F_b(s,l) = F_b(s,l) + T4xx
         endif 
      enddo 
      enddo 

c ---------------------------------------------------------------------------- 
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss3fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
      double precision t8xx
      double precision t9xx
      double precision t10xx
      double precision t11xx
      double precision t12xx
      double precision t13xx
      double precision t14xx
      double precision t15xx
      double precision t16xx
      double precision t17xx
      double precision t18xx
      double precision t19xx
      double precision t20xx
      double precision t21xx
      double precision t22xx
c
c########### CLASS 5  ############
c##################  CLASS A ##############
c SS3 
c  Do three-center part (m n |m la)->(m n |la m ),(n m|m la),(n m|la m)
c                       (m la |m n),(la m |m n ),(m la |n m),(la m |n m)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if ( l .lt. s) then 
            T1xx = HFD_a(l,s) + HFD_b(l,s)
            T5xx = 0.0 
      do n = bb1, bb2 
         if ( n .lt. s) then 
            T19xx = 0.0 
            T21xx = 0.0 
            T11xx = 0.0 
            T13xx = 0.0 
      do m = aa1, aa2 
         if ( m .lt. n) then 
         if ( m .eq. l) then 

                wint     = int_block(m,n,l,s) 
                T4xx     = HFD_a(m,n) + HFD_b(m,n)
                T2xx     = 2.0*wint*T1xx
                T5xx     = T5xx + 2.0*wint*T4xx
                T7xx     = wint*HFD_a(n,s)
                T9xx     = wint*HFD_a(n,l)
                T11xx    = T11xx + wint*HFD_a(m,s)
                T13xx    = T13xx + wint*HFD_a(m,l)
                T15xx    = wint*HFD_b(n,s)
                T17xx    = wint*HFD_b(n,l)
                T19xx    = T19xx + wint*HFD_b(m,s)
                T21xx    = T21xx + wint*HFD_b(m,l)

                F_a(m,n) = F_a(m,n) + T2xx
                F_b(m,n) = F_b(m,n) + T2xx
                F_a(n,m) = F_a(n,m) + T2xx
                F_b(n,m) = F_b(n,m) + T2xx

                F_a(m,l) = F_a(m,l) - T7xx
                F_a(l,m) = F_a(l,m) - T7xx
                F_a(m,s) = F_a(m,s) - T9xx
                F_a(s,m) = F_a(s,m) - T9xx

                F_b(m,l) = F_b(m,l) - T15xx
                F_b(l,m) = F_b(l,m) - T15xx
                F_b(m,s) = F_b(m,s) - T17xx
                F_b(s,m) = F_b(s,m) - T17xx
c
         endif 
         endif 
      enddo ! m  
            F_b(n,l) = F_b(n,l) - T19xx
            F_b(l,n) = F_b(l,n) - T19xx
            F_b(n,s) = F_b(n,s) - T21xx
            F_b(s,n) = F_b(s,n) - T21xx
            F_a(n,l) = F_a(n,l) - T11xx
            F_a(l,n) = F_a(l,n) - T11xx
            F_a(n,s) = F_a(n,s) - T13xx
            F_a(s,n) = F_a(s,n) - T13xx
         endif 
      enddo ! n  
            F_a(l,s) = F_a(l,s) + T5xx
            F_b(l,s) = F_b(l,s) + T5xx
            F_a(s,l) = F_a(s,l) + T5xx
            F_b(s,l) = F_b(s,l) + T5xx
         endif 
      enddo ! l  
      enddo ! s  
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss4fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
      double precision t8xx
      double precision t9xx
      double precision t10xx
      double precision t11xx
      double precision t12xx
      double precision t13xx
      double precision t14xx
      double precision t15xx
      double precision t16xx
      double precision t17xx
      double precision t18xx
      double precision t19xx
      double precision t20xx
      double precision t21xx
      double precision t22xx
c
c##################  CLASS B ##############
c SS4 
c  Do three-center part (m n |m la)->(m n |la m ),(n m|m la),(n m|la m)
c                       (m la |m n),(la m |m n ),(m la |n m),(la m |n m)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if ( l .lt. s) then 
            T1xx = HFD_a(l,s) + HFD_b(l,s)
            T11xx = 0.0 
            T13xx = 0.0 
            T5xx  = 0.0 
      do n = bb1, bb2 
         if ( n .eq. l) then 
            T19xx = 0.0 
            T21xx = 0.0 
      do m = aa1, aa2 
         if ( m .lt.  n) then 
         if ( m .lt. s) then 

             wint     = int_block(m,n,l,s) 
             T4xx     = HFD_a(m,n) + HFD_b(m,n)
             T2xx     = 2.0*wint*T1xx
             T7xx     = wint*HFD_a(n,s)
             T9xx     = wint*HFD_a(n,l)
             T15xx    = wint*HFD_b(n,s)
             T17xx    = wint*HFD_b(n,l)

             T19xx    = T19xx + wint*HFD_b(m,s)
             T21xx    = T21xx + wint*HFD_b(m,l)

             T5xx     = T5xx + 2.0*wint*T4xx
             T11xx    = T11xx + wint*HFD_a(m,s)
             T13xx    = T13xx + wint*HFD_a(m,l)

             F_a(m,n) = F_a(m,n) + T2xx
             F_b(m,n) = F_b(m,n) + T2xx
             F_a(n,m) = F_a(n,m) + T2xx
             F_b(n,m) = F_b(n,m) + T2xx

             F_a(m,l) = F_a(m,l) - T7xx
             F_a(l,m) = F_a(l,m) - T7xx

             F_a(m,s) = F_a(m,s) - T9xx
             F_a(s,m) = F_a(s,m) - T9xx

             F_b(m,l) = F_b(m,l) - T15xx
             F_b(l,m) = F_b(l,m) - T15xx

             F_b(m,s) = F_b(m,s) - T17xx
             F_b(s,m) = F_b(s,m) - T17xx
c
         endif 
         endif 
      enddo ! m  
             F_a(n,l) = F_a(n,l) - T11xx
             F_a(l,n) = F_a(l,n) - T11xx
             F_a(n,s) = F_a(n,s) - T13xx
             F_a(s,n) = F_a(s,n) - T13xx
             F_b(n,l) = F_b(n,l) - T19xx
             F_b(l,n) = F_b(l,n) - T19xx
             F_b(n,s) = F_b(n,s) - T21xx
             F_b(s,n) = F_b(s,n) - T21xx
         endif 
      enddo ! n  
             F_a(l,s) =  F_a(l,s) + T5xx
             F_b(l,s) =  F_b(l,s) + T5xx
             F_a(s,l) =  F_a(s,l) + T5xx
             F_b(s,l) =  F_b(s,l) + T5xx
         endif 
      enddo ! l  
      enddo ! s  
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss5fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
      double precision t8xx
      double precision t9xx
      double precision t10xx
      double precision t11xx
      double precision t12xx
      double precision t13xx
      double precision t14xx
      double precision t15xx
      double precision t16xx
      double precision t17xx
      double precision t18xx
      double precision t19xx
      double precision t20xx
      double precision t21xx
      double precision t22xx
c
c##################  CLASS C ##############
c SS5 
c  Do three-center part (m n |m la)->(m n |la m ),(n m|m la),(n m|la m)
c                       (m la |m n),(la m |m n ),(m la |n m),(la m |n m)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if ( l .lt. s) then 
            T1xx = HFD_a(l,s) + HFD_b(l,s)
            T5xx = 0.0 
      do n = bb1, bb2 
         if ( n .eq. s) then 
            T11xx = 0.0 
            T13xx = 0.0 
            T19xx = 0.0 
            T21xx = 0.0 
      do m = aa1, aa2 
         if ( m .lt. n) then 
         if ( m .lt. l) then 
c 
             wint     = int_block(m,n,l,s) 
             T4xx     = HFD_a(m,n) + HFD_b(m,n)
             T2xx     = 2.0*wint*T1xx
             T7xx     = wint*HFD_a(n,s)
             T9xx     = wint*HFD_a(n,l)
             T15xx    = wint*HFD_b(n,s)
             T17xx    = wint*HFD_b(n,l)

             T11xx    = T11xx + wint*HFD_a(m,s)
             T13xx    = T13xx + wint*HFD_a(m,l)
             T5xx     = T5xx + 2.0*wint*T4xx
             T19xx    = T19xx + wint*HFD_b(m,s)
             T21xx    = T21xx + wint*HFD_b(m,l)

             F_a(m,n) = F_a(m,n) + T2xx
             F_a(n,m) = F_a(n,m) + T2xx
             F_b(m,n) = F_b(m,n) + T2xx
             F_b(n,m) = F_b(n,m) + T2xx

             F_a(m,l) = F_a(m,l) - T7xx
             F_a(l,m) = F_a(l,m) - T7xx

             F_a(m,s) = F_a(m,s) - T9xx
             F_a(s,m) = F_a(s,m) - T9xx

             F_b(m,l) = F_b(m,l) - T15xx
             F_b(l,m) = F_b(l,m) - T15xx

             F_b(m,s) = F_b(m,s) - T17xx
             F_b(s,m) = F_b(s,m) - T17xx
c
         endif 
         endif 
      enddo ! m  
            F_b(n,l) = F_b(n,l) - T19xx
            F_b(l,n) = F_b(l,n) - T19xx
            F_a(n,l) = F_a(n,l) - T11xx
            F_a(l,n) = F_a(l,n) - T11xx
            F_a(n,s) = F_a(n,s) - T13xx
            F_a(s,n) = F_a(s,n) - T13xx
            F_b(n,s) = F_b(n,s) - T21xx
            F_b(s,n) = F_b(s,n) - T21xx
         endif 
      enddo ! n  
            F_a(l,s) = F_a(l,s) + T5xx
            F_b(l,s) = F_b(l,s) + T5xx
            F_a(s,l) = F_a(s,l) + T5xx
            F_b(s,l) = F_b(s,l) + T5xx
         endif 
      enddo ! l  
      enddo ! s  
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss6fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
c
c########### CLASS 1  ############
c SS6 
c  Do one-center part (m m |m m)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
            T1xx = HFD_a(l,s) + HFD_b(l,s)
      do n = bb1, bb2 
      do m = aa1, aa2 
         if ( n .eq. m ) then 
         if ( s .eq. m ) then 
         if ( l .eq. m ) then 
        
            wint     = int_block(m,n,l,s) 
            T2xx     = wint*T1xx
            T3xx     = wint*HFD_a(n,s)
            T4xx     = wint*HFD_b(n,s)

            F_a(m,n) = F_a(m,n) + T2xx
            F_b(m,n) = F_b(m,n) + T2xx
            F_a(m,l) = F_a(m,l) - T3xx
            F_b(m,l) = F_b(m,l) - T4xx
c
         endif 
         endif 
         endif 
      enddo ! m  
      enddo ! n  
      enddo ! l  
      enddo ! s   
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss7fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
c
c########### CLASS 2  ############
c SS7 
c  Do two-center part (m m |n n)->(n n |m m )
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if ( l .eq. s) then 
            T1xx = HFD_a(l,s) + HFD_b(l,s)
            T4xx = 0.0 
      do n = bb1, bb2 
      do m = aa1, aa2 
         if ( m .eq. n) then 
         if ( m .lt. l) then 
 
             wint     = int_block(m,n,l,s) 
             T3xx     = HFD_a(m,n) + HFD_b(m,n)
             T2xx     = wint*T1xx
             T4xx     = T4xx + wint*T3xx
             T5xx     = wint*HFD_a(n,s)
             T7xx     = wint*HFD_b(n,s)

             F_a(m,n) = F_a(m,n) + T2xx
             F_b(m,n) = F_b(m,n) + T2xx

             F_a(m,l) = F_a(m,l) - T5xx
             F_a(l,m) = F_a(l,m) - T5xx

             F_b(m,l) = F_b(m,l) - T7xx
             F_b(l,m) = F_b(l,m) - T7xx
c
         endif 
         endif 
c
      enddo ! m  
      enddo ! n  
            F_a(l,s) = F_a(l,s) + T4xx
            F_b(l,s) = F_b(l,s) + T4xx
         endif 
      enddo ! l  
      enddo ! s  
c
c########### END CLASS 2  ############
c
      return 
      end 
c
c ---------------------------------------------------------------------------- 

      subroutine form_ss8fock(int_block,n_basis,aa1,aa2,bb1,bb2,cc1,cc2,
     *                      dd1,dd2,HFD_A,HFD_B,F_A,F_B)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s, l  
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision F_A(n_basis,n_basis) 
      double precision F_B(n_basis,n_basis) 
      double precision F1A(n_basis) 
      double precision F1B(n_basis) 
      double precision tempa, tempb, temp_rs, temp_mn   

      double precision wint  
      double precision t1xx
      double precision t2xx
      double precision t3xx
      double precision t4xx
      double precision t5xx
      double precision t6xx
      double precision t7xx
      double precision t8xx
      double precision t9xx
      double precision t10xx
      double precision t11xx
c
c########### CLASS 4  ############
c SS8 
c  Do two-center part (m n |m n)->(m n |n m ),(n m|m n),(n m|n m)
c
      do s = dd1, dd2 
      do l = cc1, cc2 
         if ( l .lt. s ) then 
            T1xx = HFD_a(l,s) + HFD_b(l,s)
      do n = bb1, bb2 
         if ( n .eq. s) then 
            T6xx  = 0.0 
            T7xx  = 0.0 
            T10xx = 0.0 
            T11xx = 0.0 
      do m = aa1, aa2 
         if ( m .lt. n) then 
         if ( m .eq. l) then 
 
             wint     = int_block(m,n,l,s) 
             T2xx     = 2.0*wint*T1xx

             T4xx     = wint*HFD_a(n,s)
             T5xx     = wint*HFD_a(n,l)
             T8xx     = wint*HFD_b(n,s)
             T9xx     = wint*HFD_b(n,l)

             T6xx     = T6xx  + wint*HFD_a(m,s)
             T7xx     = T7xx  + wint*HFD_a(m,l)
             T10xx    = T10xx + wint*HFD_b(m,s)
             T11xx    = T11xx + wint*HFD_b(m,l)
 
             F_a(m,n) = F_a(m,n) + T2xx
             F_b(m,n) = F_b(m,n) + T2xx
             F_a(n,m) = F_a(n,m) + T2xx
             F_b(n,m) = F_b(n,m) + T2xx

             F_a(m,l) = F_a(m,l) - T4xx
             F_a(m,s) = F_a(m,s) - T5xx
             F_b(m,l) = F_b(m,l) - T8xx
             F_b(m,s) = F_b(m,s) - T9xx
c
         endif 
         endif 
      enddo ! m  
            F_a(n,l) = F_a(n,l) - T6xx
            F_a(n,s) = F_a(n,s) - T7xx
            F_b(n,l) = F_b(n,l) - T10xx
            F_b(n,s) = F_b(n,s) - T11xx
         endif 
      enddo ! n  
         endif 
      enddo ! l  
      enddo ! s  
c
c ---------------------------------------------------------------------------- 
c
      return 
      end 

      subroutine form_ss(int_block,m,n,r,s,n_basis,aa1,aa2,bb1,bb2,
     *                   cc1,cc2,dd1,dd2)
      implicit none 
      integer aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2, n_basis 
      integer a, b, c, d 
      integer m, n, r, s   
      double precision int_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  
      double precision temp_block(aa1:aa2,bb1:bb2,cc1:cc2,dd1:dd2)  

      if ((m .lt. n) .and. (r .lt. s)) then 
         do a = aa1, aa2 
         do b = bb1, bb2 
         do c = cc1, cc2 
         do d = dd1, dd2 
            int_block(b,a,c,d) = int_block(a,b,c,d)  
            int_block(a,b,d,c) = int_block(a,b,c,d)  
            int_block(b,a,d,c) = int_block(a,b,c,d)  
         enddo 
         enddo 
         enddo 
         enddo 
      endif 

      if ((m .lt. n) .and. (r .eq. s)) then 
         do a = aa1, aa2 
         do b = bb1, bb2 
         do c = cc1, cc2 
         do d = dd1, dd2 
            int_block(b,a,c,d) = int_block(a,b,c,d)  
            int_block(a,b,d,c) = int_block(a,b,c,d)  
         enddo 
         enddo 
         enddo 
         enddo 
      endif 

      if ((m .eq. n) .and. (r .lt. s)) then 
         do a = aa1, aa2 
         do b = bb1, bb2 
         do c = cc1, cc2 
         do d = dd1, dd2 
            int_block(a,b,d,c) = int_block(a,b,c,d)  
            int_block(b,a,d,c) = int_block(a,b,c,d)  
         enddo 
         enddo 
         enddo 
         enddo 
      endif 

      return 
      end 

      subroutine hfdensity_copy(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis)  
      implicit none 

      integer n_basis, a, b 
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision HFDOLD_A(n_basis,n_basis) 
      double precision HFDOLD_B(n_basis,n_basis) 

      do a = 1, n_basis 
      do b = 1, n_basis 
         HFDOLD_A(b,a) = HFD_A(b,a) 
         HFDOLD_B(b,a) = HFD_B(b,a) 
      enddo 
      enddo 

      return 
      end 


      subroutine check_conv(HFD_A,HFD_B,HFDOLD_A,HFDOLD_B,n_basis,flag)  
      implicit none 

      integer n_basis, a, b, flag  
      double precision HFD_A(n_basis,n_basis) 
      double precision HFD_B(n_basis,n_basis) 
      double precision HFDOLD_A(n_basis,n_basis) 
      double precision HFDOLD_B(n_basis,n_basis) 
      double precision max, diffa, diffb, thresh    

      max = 0.0 
      do a = 1, n_basis 
      do b = 1, n_basis 
         diffa = dabs(HFDOLD_A(b,a) - HFD_A(b,a))  
         diffb = dabs(HFDOLD_B(b,a) - HFD_B(b,a))  
         if (diffa .gt. max) max = diffa 
         if (diffb .gt. max) max = diffb 
      enddo 
      enddo 

      thresh = 1.0d-6 
      flag = 0 
      if (max .lt. thresh) flag = 1 
c      write(6,*) '      Max density difference :', max 

c     if (flag .eq. 0) then 
c        do a = 1, n_basis 
c        do b = 1, n_basis 
c           HFDOLD_A(b,a) = HFD_A(b,a) 
c           HFDOLD_B(b,a) = HFD_B(b,a) 
c        enddo 
c        enddo 
c     endif 

      return 
      end 

