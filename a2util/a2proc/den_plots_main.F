#include "flags.h"
      subroutine Den_Plots_Main(Icore, Icrsiz)
c
      implicit double precision (a-h,o-z)
      parameter (mxcoef=30)
c
#include "machsp.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "baslims.par"
#include "offsets_2getden.com"
C
      character*32 szFile
      logical bExist
      Character*4 Comp_pgrp, Full_pgrp
      Dimension Nocc(16), Atommass(Mxatms), Iatmchrg(Mxatms),  
     &          Coord(3*Mxatms), Norbits_fullG(Mxatms), 
     &          NOrbits_compG(Mxatms),Nbsfns_4irrep(8)
      Dimension Icore(Icrsiz)
C
      Data Ione, Ieight /1, 8/
C
      Iuhf = 1
      If (iflags(h_IFLAGS_ref).eq.0) iuhf = 0 
C
      Maxcor   = Icrsiz
      Mxangmom = Nht
      Length   = 0
c
c Read the JOBARC file for basic data of the molecule. 
c
      
      Call Getrec(-1, 'JOBARC', 'NUMDROPA', Length, Ijunk)       
      If (Length .GT. 0) Then 
         Print*, "Frozen-core is not allowed density plots"
         Call Errex
      Endif
c 
      Call Getrec(20, 'JOBARC', 'NREALATM', Ione, Nreal_atoms)
      Call Getrec(20, 'JOBARC', 'NATOMS  ', Ione, Natoms)
      Call Getrec(20, 'JOBARC', 'FULLNORB', Ione, Iful_unq_atoms)
      Call Getrec(20, 'JOBARC', 'COMPNORB', Ione, Icmp_unq_atoms)
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Coord)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
      Call Getrec(20, 'JOBARC', 'COMPNIRR', Ione, Nirrep)
      Call Getrec(20, 'JOBARC', 'OCCUPYA ', Nirrep, Nocc(1))
      Call Getrec(20, 'JOBARC', 'NBASTOT ', Ione, Nbfns)
      Call Getrec(20, 'JOBARC', 'NAOBASFN', Ione, Naobfns)
      Call Getrec(20, 'JOBARC', 'NUMBASIR', Nirrep, Nbsfns_4irrep)
      Call Getrec(20, 'JOBARC', 'FULLPOPV', Iful_unq_atoms, 
     &            Norbits_fullG)
      Call Getrec(20, 'JOBARC', 'COMPPOPV', Icmp_unq_atoms, 
     &            Norbits_compG)
C
      nOCCa = 0
      Do Irrep = 1, Nirrep
         Nocca = Nocca + Nocc(Irrep)
      End Do
      If (Iuhf .EQ. 1) Then
         Call Getrec(20, 'JOBARC', 'OCCUPYB', Nirrep, Nocc(I+8)) 
         Noccb = 0
         Do Irrep = 1, Nirrep
            Noccb = Noccb + Nocc(8+I) 
         End Do
      Else
         Noccb = Nocca
      Endif 
C
      Call Getcrec(20, 'JOBARC', 'COMPPTGP', 4, Comp_pgrp)
      Call Getcrec(20, 'JOBARC', 'FULLPTGP', 4, Full_pgrp)
C
#ifdef _DEBUG_LVLM2
      print*, "Entering A2_PREP4_MOL_READ.."
      Print*, (Norbits_compG(i), i=1, Icmp_unq_atoms)
      Print*, (Norbits_fullG(i), i=1, Iful_unq_atoms)
#endif 
      Call A2_PREP4_MOL_READ(Icore, Maxcor, Natoms, Nbfns, 
     &                       Naobfns, Icmp_unq_atoms, 
     &                       Iful_unq_atoms, Coord, Iatmchrg,
     &                       Norbits_fullG, Norbits_compG, 
     &                       Mxshel, Mxangmom, Itot_prim, Iunqshl, 
     &                       Ntotshl, Max_prim_4atom, 
     &                       Max_prim_4shell, Inext)

#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "Offsets for various arrays initilized in Pre4_Den_plo."
      Write(*, '(31I5)'), IMEMBC,INUC,NUFCT,NUMOM,NFCT,NANGMOM,
     &                   IATMNAM,NAOUATM,NAOATM,IANGMOMSHL,
     &                   ICONFUNSHL,IPRMFUNSHL,INPRIMSHL,
     &                   IANGMOMTSHL,ICONFUNTSHL,IOFFSETPRM,
     &                   IOFFSETCON,IOFFSETSHL,IMAP_SHL2CNT,
     &                   IPRMFUNTSHL,NMOMAO,IALPHA,ICOEFFA,
     &                   ICOEFFB,IDENSA,IDENSB,IPCOEFFA,
     &                   IPCOEFFB, IPCOEFF,ISCR1,ISCR2
C
      Write(6,*) 
      Print*, "The orbital exponents in main", Itot_prim, ialpha
      call output(icore(ialpha), 1, Itot_prim, 1, 1, 1, Itot_prim,
     &            1,1,1)  
      Print*, "The Contraction Coefs."
      call output(icore(ipcoeff), 1, Itot_prim, 1, Naobfns, Itot_prim,
     &            Naobfns, 1)
#endif
#ifdef _ARCHIVED
C This block is archived. Obtaining the density need to be handled
C somewhat genrelly and now it is done inside generate_12Dgrid.F 
C using a2get_den.F
C
      Inturldens = Inext
      Iscfdens_t = Inturldens + Nbfns*Nbfns*Iintfp
      Iscfdens_d = Iscfdens_t + Naobfns*Naobfns*Iintfp
      Icordens_t = Iscfdens_d + Naobfns*Naobfns*Iintfp 
      Icordens_d = Icordens_t + Naobfns*Naobfns*Iintfp 
      Iscfcoef_a = Icordens_d + Naobfns*Naobfns*Iintfp  
      Iscfcoef_b = Iscfcoef_a + Naobfns*Nbfns*Iintfp 
      Isphr2cart = Iscfcoef_b + Naobfns*Nbfns*Iintfp
      Iscftmp1   = Isphr2cart + Naobfns*Nbfns*Iintfp
      Iscftmp2   = Iscftmp1   + Nbfns*Nbfns*Iintfp 
      Inext      = Iscftmp2   + Nbfns*Nbfns*Iintfp 
c      
#ifdef _DEBUG_LVLM2
      print*, "Entering A2get_den."
#endif 
      Call A2get_den(Icore(Iscfcoef_a), Icore(Inturldens),
     &               Icore(Iscfdens_t), Icore(Icordens_t), 
     &               Icore(Isphr2cart), Icore(IScftmp1), Naobfns,
     &               Nbfns, Nocc, Nbfns_4irrep, Iuhf+1)
     &               
      If (Iuhf .GT. 0) Then

         Call A2get_den(Icore(Iscfcoef_b), Icore(Inturldens),
     &                  Icore(Iscfdens_d), Icore(Icordens_d), 
     &                  Icore(Isphr2cart), Icore(IScftmp1), Naobfns,
     &                  Nbfns, Nocc, Nbsfn_4irrep, Iuhf+1)
     &               
      Endif 
C
      Call Generate_12Dgrid(Icore(Inext), Ileft, Natoms, Itot_prim,
     &                      Nbfns, Naobfns, Lnp1, Norbits_compG,
     &                      Iful_unq_atoms, Coord, Icore(Ialpha),
     &                      Icore(Ipcoeff), Iunqshl,Icore(Inprimshl),
     &                      Icore(Iangmomshl), Icore(Iconfunshl),
     &                      Icore(Iprmfunshl), Icore(Nufct),
     &                      Icore(Naouatm), Icore(Ioffsetprm),
     &                      Icore(Ioffsetcon), Icore(Itmp1),
     &                      Icore(Itmp2), Icore(Itmp3),
     &                      Icore(Iscfdens_t), Icore(Iscfdens_d),
     &                      Icore(Icordens_t), Icore(Icordens_d),
     &                      Icore(Iprductint), Iatmchrg, Iuhf)
c
#endif
C
       Lnp1             = Max_prim_4shell
       Icnterbf         = Inext 
       Iscfcoef_a       = Icnterbf         + Nbfns
       Iscfcoef_b       = Iscfcoef_a       + Naobfns*Nbfns*Iintfp
       Iscfcoef_reord_a = Iscfcoef_b       + Naobfns*Nbfns*Iintfp 
       iSCfcoef_reord_b = Iscfcoef_reord_a + Nbfns*Naobfns*Iintfp
       Icoef_a          = Iscfcoef_reord_b + Nbfns*Naobfns*Iintfp
       Icoef_b          = Icoef_a          + Nbfns*Itot_prim*Iintfp
       Itmp1            = Icoef_b          + Nbfns*Itot_prim*Iintfp
       Itmp2            = Itmp1            + Lnp1*Lnp1*Iintfp
       Itmp3            = Itmp2            + Lnp1*Lnp1*Iintfp
       Iprductint       = Itmp3            + Lnp1*Lnp1*Iintfp
       Inuctr           = Iprductint       + Naobfns*Naobfns*Iintfp
       Inext            = Inuctr           + Natoms
       Ileft            = Maxcor - Inext
c
#ifdef _DEBUG_LVLM2
      print*, "Entering Reord_MOs:"
      Print*, "Input variables", Inext, Max_prim_4atom, Itot_prim, 
     &         icnterbf, Iscfcoef_reord_a, Iscfcoef_reord_b,
     &         Icoef_a, Icoef_b, Nangmom, Nmomao, Naoatm,ipcoeff
      write(*,*)
#endif 
C
       Call Reord_mos(0, Nreal_atmos, Naobfns, Nbfns, Max_prim_4atom, 
     &                Icore(Icnterbf), ICore(Iscfcoef_a),  
     &                ICore(Iscfcoef_reord_a), Icore(Nangmom),
     &                Icore(Nmomao), Icore(Naoatm), .False.)
C
      CALL Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &            Icore(Ipcoeff),Itot_prim,Icore(Iscfcoef_reord_a),
     &            Naobfns,0.0D+00,Icore(Icoef_a),Itot_prim)


      If (Iuhf .ne. 0) Then
         Call Reord_mos(1, Nreal_atoms, Naobfns, Nbfns, 
     &                  Max_prim_4atom, Icore(Icnterbf), 
     &                  ICore(Iscfcoef_b), ICore(Iscfcoef_reord_b), 
     &                  Icore(Nangmom), Icore(Nmomao), Icore(Naoatm),
     &                  .False.)
C
         Call Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &               Icore(Ipcoeff),Itot_prim,
     &               Icore(Iscfcoef_reord_b), Naobfns,0.0D+00,
     &               Icore(Icoef_b),Itot_prim)
      End if

#ifdef _DEBUG_LVLM2
      print*, "Entering Nodummy:"
      write(*,*)
#endif 
      CALL Nodummy(Natoms, Max_prim_4atom, Iatmchrg, Coord,
     &            Icore(Nangmom), Icore(Nmomfct), Icore(Inuctr))
C
      if (Comp_pgrp .Eq. "C1") Iful_unq_atoms = Icmp_unq_atoms
C
#ifdef _DEBUG_LVLM2
      print*, "Entering Generate_12Dgrid:"
      Write(6,*)
      Print*, "The start of Icore", Inext, ialpha, Ileft,
     &         Iful_unq_atoms, Ipcoeff, Iunqshl, Inprimshl,
     &         Iangmomshl, Iconfunshl, Iprmfunshl, Nufct,
     &         Naouatm, Ioffsetprm, Ioffsetcon, Ioffsetprm,
     &         Ioffsetcon, Itmp1, Itmp2, Itmp3, Iprductint
C   
      Print*, "The exponents at entry 2 Gen_12Grid" 
      call output(icore(ialpha), 1, Itot_prim, 1, 1, 1, Itot_prim,
     &            1,1,1)
      Write(6,*)
      print*, "The offsets for primitives and contr. functions"
      Write(6,"(4i5)") (Icore(Nufct -1 +i), i=1, Iful_unq_atoms)
      Write(6,"(4i5)") (Icore(naouatm -1 +i), i=1, Iful_unq_atoms)
      Write(6,*)
      Write(6,"(4i5)") (Icore(Iangmomshl -1 +i), 
     &                 i=1,Iunqshl*Iful_unq_atoms)
      Write(6,"(4i5)") (Icore(Inprimshl -1 +i), 
     &                 i=1,iful_unq_atoms)
      Write(6,"(4i5)") (Icore(Iconfunshl -1 +i), 
     &                 i=1,iful_unq_atoms*Iunqshl)
      Write(6,"(4i5)") (Icore(Iprmfunshl -1 +i), 
     &                 i=1,iful_unq_atoms*Iunqshl)
      Write(6,"(4i5)") (Icore(Ioffsetprm -1 +i), 
     &                 i=1,iful_unq_atoms*Iunqshl)
      Write(6,"(4i5)") (Icore(Ioffsetcon-1 +i), 
     &                 i=1,iful_unq_atoms*Iunqshl)

#endif 
      Call Generate_12Dgrid(Icore(Inext), Ileft, Natoms, Itot_prim,
     &                      Nbfns, Naobfns, Lnp1, Norbits_compG,
     &                      Iful_unq_atoms, Coord, Icore(Ialpha),
     &                      Icore(Ipcoeff), Iunqshl,Icore(Inprimshl), 
     &                      Icore(Iangmomshl), Icore(Iconfunshl), 
     &                      Icore(Iprmfunshl), Icore(Nufct), 
     &                      Icore(Naouatm), Icore(Ioffsetprm),
     &                      Icore(Ioffsetcon), Icore(Itmp1),
     &                      Icore(Itmp2), Icore(Itmp3),
     &                      Icore(Iprductint), Iatmchrg, Iuhf)
      write(6,*)
      Write(6,"(t3,50a)")"@-den_plot_main: The density distribution ",
     &                  "calculation is succesfully completed."
      Return
      End

