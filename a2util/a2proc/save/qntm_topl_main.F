#include "flags.h"
      subroutine Qntm_Topl_Main(Icore, Icrsiz)
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
#ifdef _DEBUG_LVL0
      print*, "Entering a2_prep4_mol_read"
      Print*, (Norbits_compG(i), i=1, Icmp_unq_atoms)
      Print*, (Norbits_fullG(i), i=1, Iful_unq_atoms)
#endif 
      Call a2_prep4_mol_read(Icore, Maxcor, Natoms, Nbfns, 
     &                       Naobfns, Icmp_unq_atoms, 
     &                       Iful_unq_atoms, Coord, Iatmchrg, 
     &                       Norbits_fullG, Norbits_compG, 
     &                       Mxshel, Mxangmom, Itot_prim, Iunqshl, 
     &                       Ntotshl, Max_prim_4atom, 
     &                       Max_prim_4shell, Inext)

#ifdef _DEBUG_LVL0
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
C
      Icordens_a = Inext
      Icordens_b = Icordens_a + Nbfns*Nbfns*Iintfp
      Inatorbs_a = Icordens_b + Nbfns*Nbfns*Iintfp
      Inatorbs_b = Inatorbs_a + Naobfns*Nbfns*Iintfp
      Ieigvecs_a = Inatorbs_b + Naobfns*Nbfns*Iintfp
      Ieigvecs_b = Ieigvecs_a + Nbfns*Nbfns*Iintfp
      Ieigvals_a = Ieigvecs_b + Nbfns*Nbfns*Iintfp
      Ieigvals_b = Ieigvals_a + Nbfns*Iintfp
      Iscftmp1   = Ieigvals_b + Nbfns*Iintfp
      Iscftmp2   = Iscftmp1   + Nbfns*Nbfns*Iintfp 
      Iocca      = Iscftmp2   + Nbfns*Nbfns*Iintfp 
      Ioccb      = Iocca      + Nbfns*Iintfp
      Inext      = Ioccb      + Nbfns*Iintfp
#ifdef _DEBUG_LVL0
      print*, "Entering A2make_wfn."
#endif 
      Call A2make_wfn(Icore(Icordens_a), Icore(Icordens_b),
     &                Icore(Inatorbs_a), Icore(Inatorbs_b),
     &                Icore(Ieigvecs_a), Icore(Ieigvecs_b),
     &                ICore(Ieigvals_a), Icore(Ieigvals_b), 
     &                Icore(Iscftmp1), Icore(Iscftmp2), 
     &                ICore(Iocca), Icore(Ioccb), Nocca, Noccb,
     &                Maxocca, Maxoccb, Naobfns, Nbfns, Iuhf)
C
       Lnp1             = Max_prim_4shell
       Icnterbf         = Inext 
       Inatorbs_reord_a = Icnterbf         + Nbfns
       Inatorbs_reord_b = Inatorbs_reord_a + Naobfns*Nbfns*Iintfp 
       Icoef_a          = Inatorbs_reord_b + Nbfns*Naobfns*Iintfp
       Icoef_b          = Icoef_a          + Nbfns*Itot_prim*Iintfp
       Itype            = Icoef_b          + Nbfns*Itot_prim*Iintfp
       ICent            = Itype            + Itot_prim*Iintfp
       INext            = ICent            + Itot_prim*Iintfp
c
#ifdef _DEBUG_LVLM2
      print*, "Entering Reord_MOs:"
      Print*, "Input variables", Inext, Max_prim_4atom, Itot_prim, 
     &         icnterbf, Iscfcoef_reord_a, Iscfcoef_reord_b,
     &         Icoef_a, Icoef_b, Nangmom, Nmomao, Naoatm,ipcoeff
      write(*,*)
#endif 
C
       Call Reord_mos(0, Nreal_atoms, Naobfns, Nbfns, Max_prim_4atom, 
     &                Icore(Icnterbf), ICore(Inatorbs_a),  
     &                ICore(Inatorbs_reord_a), Icore(Nangmom),
     &                Icore(Nmomao), Icore(Naoatm),.True.)
C
      CALL Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &            Icore(Ipcoeff),Itot_prim,Icore(Inatorbs_reord_a),
     &            Naobfns,0.0D+00,Icore(Icoef_a),Itot_prim)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) "Alpha orbital for Quantum Topology analysis"
      Call output(Icore(Icoef_a), 1, Itot_prim, 1, Nbfns, Itot_prim,
     &            Nbfns, 1)
#endif

      If (Iuhf .ne. 0) Then
         Call Reord_mos(1, Nreal_atoms, Naobfns, Nbfns, 
     &                  Max_prim_4atom, Icore(Icnterbf), 
     &                  ICore(Inatorbs_b), ICore(Inatorbs_reord_b), 
     &                  Icore(Nangmom), Icore(Nmomao), Icore(Naoatm),
     &                  .True.)
C
         Call Xgemm('N','N',Itot_prim,Nbfns,Naobfns,1.0D+00,
     &               Icore(Ipcoeff),Itot_prim,
     &               Icore(Inatorbs_reord_b), Naobfns,0.0D+00,
     &               Icore(Icoef_b),Itot_prim)
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,*) "Beta orbital for Quantum Topology analysis"
      Call output(Icore(Icoef_b), 1, Itot_prim, 1, Nbfns, Itot_prim,
     &            Nbfns, 1)
#endif

      End if

#ifdef _DEBUG_LVLM2
      print*, "Entering Nodummy:"
      write(*,*)
#endif 
      Call Nodummy(Natoms, Max_prim_4atom, Iatmchrg, Coord,
     &            Icore(Nangmom), Icore(Nmomfct), Icore(Inuctr))
C
      if (Comp_pgrp .Eq. "C1") Iful_unq_atoms = Icmp_unq_atoms
C
#ifdef _DEBUG_LVL0            
      Write(6,*)              
      Write(6,*) "The number of contracted function per shell"
      Write(6,*) (ICORE(ICONFUNTSHL - 1 + I), I=1, 7)
      Write(6,*) "The shell angular momentum"
      Write(6,*) (ICORE(IANGMOMTSHL - 1 + I), I=1, 7)
      Write(6,*) "The number of primitive functions per shell"
      Write(6,*) (ICORE(IPRMFUNTSHL - 1 + I), I=1, 7)
#endif
      Call A2write_wfn(Icore(Icoef_a), Icore(Nangmom), Icore(Ialpha),
     &                 Icore(Iprmfunshl), Icore(Iangmomshl), 
     &                 Icore(Inprimshl), Icore(Itype), Icore(Icent),
     &                 ICore(Ieigvals_a), Icore(Iocca),
     &                 Coord, Maxocca, Mxshel, Itot_prim, 
     &                 Nbfns, Nreal_atoms,
     &                 Iatmchrg)
C
      If (Iuhf .ne. 0)  Call A2write_wfn(Icore(Icoef_a), 
     &                                   Icore(Nangmom), 
     &                                   Icore(Ialpha) , 
     &                                   Icore(Iprmfunshl)   , 
     &                                   Icore(Iangmomshl),
     &                                   Icore(Iprmfunshl),
     &                                   Icore(Itype), Icore(Icent),
     &                                   ICore(Ieigvals_b),
     &                                   Icore(Iocca), 
     &                                   Coord, Maxoccb, 
     &                                   Mxshel, Itot_prim, Nbfns, 
     &                                   Nreal_atoms, Iatmchrg)
C
#ifdef _DEBUG_LVLM1
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
      Iunqshl=7
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
      Return
      End

