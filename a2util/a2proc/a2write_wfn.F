#include "flags.h"
      Subroutine a2write_wfn(Wfn, Nangmom, Alpha, Nprim_4shl, 
     &                       Iangmom, IShell_4atm, Itype, Icent, 
     &                       Eigvals, Occnm, Iorder, Scr, Coord,
     &                       Maxocc, Nshls, Itot_prim, Nbfns,
     &                       Nreal_atoms, Iatmchrg, Iunit, Iuhf,
     &                       Iwrite_H, Iwrite_T, Ekin)
      
      Implicit Double Precision (A-H, O-Z)
C
#include "machsp.com"
#include "jodaflags.com"
C
      Character*2 Atmlab(103)
      Character*11 Method(45)
      Character*4 Ref 
      Logical Iwrite_H, Iwrite_T
      
      Character Title(10)
      Integer Angmom

      Dimension Coord(3*Nreal_atoms), Iatmchrg(Nreal_atoms),
     &          Wfn(Itot_prim, Nbfns), Nprim_4shl(Nreal_atoms*NShls),
     &          Alpha(Itot_prim), Nangmom(Nreal_atoms), 
     &          Itype(Itot_prim), Iangmom(Nreal_atoms*Nshls), 
     &          IShell_4atm(Nreal_atoms), Icent(Itot_prim),
     &          Eigvals(Nbfns), Occnm(Nbfns), Iorder(Itot_prim),
     &          Scr(Itot_prim)
C 
      Data (Atmlab(i),i=1,103)
     & /'H ','HE',
     & 'LI','BE','B ','C ','N ','O ','F ','NE',
     & 'NA','MG','AL','SI','P ','S ','CL','AR',
     & 'K ','CA',
     & 'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     & 'GA','GE','AS','SE','BR','KR',
     & 'RB','SR',
     & 'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',
     & 'IN','SN','SB','TE','I ','XE',
     & 'CS','BA','LA',
     & 'CE','PR','ND','PM','SM','EU','GD',
     & 'TB','DY','HO','ER','TM','YB','LU',
     & 'HF','TA','W ','RE','OS','IR','PT','AU','HG',
     & 'TL','PB','BI','PO','AT','RN',
     & 'FR','RA','AC',
     & 'TH','PA','U ','NP','PU','AM','CM',
     & 'BK','CF','ES','FM','MD','NO','LR'/
C
      Data (Method(I), i=1, 45)/
     &  "SCF",          "MBPT(2)",      "MBPT(3)",      "SDQ-MBPT(4)",
     &  "MBPT(4)",      "LCCD",         "LCCSD",        "UCCSD(4)",
     &  "CCD",          "UCC(4)",       "CCSD",         "CCSD[T]",
     &  "CCSD+TQ*",     "CCSDT-1",      "CCSDT-1b",     "CCSDT-2",
     &  "CCSDT-3",      "CCSDT-4",      "CCSDT",        "LCCSDT",
     &  "CCD+ST(CCD)",  "QCISD(T)",     "CCSD(T)",      "QCISD",
     &  "CID",          "CISD",         "QCISD(TQ)",    "CCSD(TQ)",
     &  "CCSD+TQ",      "CCSDT+Q*",     "CCSDT+Q",      "CC5SD(T)",
     &  "CCSD-T",       "CC3",          "CCSDT-T1T2",   "CCSDTQ-1",
     &  "CCSDTQF-1",    "CCSDTQ-2",     "CCSDTQ-3",     "CCSDTQ",
     &  "ACCSD",        "HFDFT",        "ACCSD(T)",     "CCSD(TQf)",
     &  "CCSDT(Qf)" /
C    
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a)") "The number of functions on each atom"
      Write(6, "(4(1x,I3))") (Nangmom(I), I=1, Nreal_atoms)
      Write(6,*)
      Write(6, "(a)") "The number of prim. functions on each atom"
      Do I =1, Nreal_atoms
      Write(6, "(4(1x,I3))") (NPrim_4shl(ioff+J),J=1,IShell_4atm(I))
      ioff = ioff + IShell_4atm(I)
      Enddo 
      Write(6,*) 
      Write(6,*) "The number of shells per atom"
      Write(6, "(4(1x,I3))") (IShell_4atm(I), I=1, Nreal_atoms)
      Write(6,*) 
      Write(6,*) "The shell angular momentum per atom s=1,p=3,d=6,,"
      Ioff = 0 
      Do I= 1, Nreal_atoms
      Write(6, "(4(1x,I3))") (Iangmom(ioff + J), J=1,IShell_4atm(I))
      ioff = ioff + IShell_4atm(I)
      Enddo
      Write(6,*)
      Write(6, "(3F10.5)"), (coord(i), i=1,3*Nreal_atoms)
      Write(*,*)
#endif
      Call Built_Reord_Vec(Nreal_atoms, Nprim_4shl, Iangmom,
     &                     IShell_4atm, Itot_prim, Nshls, Iorder)

      Call reord_vecs(Alpha, Scr, Iorder, Itot_prim, .True.)
      Call reord_vecs(Eigvals, Scr, Iorder, Itot_prim, .True.)
C
      Do Iorbs = 1, Nbfns
         Call reord_vecs(Wfn(1,Iorbs), Scr, Iorder, Itot_prim, .True.)
      Enddo

      Ioff = 0
      Ifun = 0
      Do Iatms = 1, Nreal_atoms
         Do Ishl = 1, IShell_4atm(Iatms)
            Angmom = Iangmom(Ifun+Ishl)
            Nconfs = Nprim_4shl(Ifun+Ishl)
            If (Angmom .Eq. 1)  Then
                Imin = 1
                Imax = 1
            Elseif (Angmom .Eq. 3) Then
                Imin = 2
                Imax = 4
            Elseif (Angmom .Eq. 6) Then
                Imin = 5
                Imax = 10
            Elseif (Angmom .Eq. 10) Then
                Imin = 11 
                Imax = 20
            Elseif (Angmom .Eq. 15) Then
                Imin = 21 
                Imax = 35
            Elseif (Angmom .Eq. 21) Then
                Imin = 36
                Imax = 56
            Elseif (Angmom .Eq. 36) Then
                Imin = 57
                Imax = 84
            Elseif (Angmom .Gt. 36) Then
                Write(6, "(a)") "The WFN file is only available",
     &                          " for up to L angular momentum" 
                Call Errex
            Endif
  
            Do Iang = Imin, Imax
               Do Nfuncs = 1, Nconfs
                  Ioff = Ioff + 1
                  Itype(Ioff) = Iang
                  Icent(Ioff) = Iatms
               Enddo 
            Enddo
         Enddo
         Ifun = Ifun + IShell_4atm(Iatms)
      Enddo
C           
      If (Iwrite_H) Then
         Write(Iunit, "(a,a)") "----Input file for quantum ",
     &                         "topology analysis----"
C
         Write(IUnit, 100) Maxocc, Itot_prim, Nreal_atoms
  100    Format('GAUSSIAN',10x,I5,' MOL ORBITALS',1x,I6,
     &          ' PRIMITIVES', 4x,I5,' NUCLEI')

         Ioff = 0
         Do Iatms = 1, Nreal_atoms
             Write(Iunit, 110)  Atmlab(Iatmchrg(Iatms)), Iatms, 
     &                          Iatms, (Coord(J+ioff), J=1,3), 
     &                          Dble(Iatmchrg(Iatms))
             Ioff = Ioff + 3
         Enddo
C       
  110    Format(a4,i4,4x,'(CENTRE',i3,')',1x,3f12.8, 
     &          '  CHARGE =',f5.1)
C
         Write(Iunit, 120) (Icent(I), I = 1, Itot_prim)
         Write(Iunit, 130) (Itype(I), I = 1, Itot_prim)
         Write(Iunit, 140) (Alpha(I), I = 1, Itot_prim)

  120    Format('CENTRE ASSIGNMENTS',2x,20i3)
  130    Format('TYPE ASSIGNMENTS',4x,20i3)
  140    Format('EXPONENTS',1x,1p,5e14.7)
      Endif
     
      Do Iorb = 1, Maxocc
         Write(Iunit, 150) Iorb, Occnm(Iorb), Eigvals(Iorb)
         Write(Iunit, 160) (Wfn(Iprim, Iorb), Iprim =1, Itot_prim)
      Enddo 

  150 Format('MO  ',i3,19x,'OCC NO =',f13.7,'  ORB. ENERGY = ',f12.7 )
  160 Format(1p,5e16.8)
  
      If (Iwrite_T) Then
         Write(Iunit, "(a)") "END DATA"
C
         Call Getrec(20, "JOBARC", 'SCFENEG ', Iintfp, Escf)
         Call Getrec(20, "JOBARC", 'TOTENERG', Iintfp, Etot)
         Vir     = 0.0D0
C
         Imethod = Iflags(h_IFLAGS_calc) + 1
         Iref    = Iflags(h_IFLAGS_ref)
         If (Iref .Eq. 0) Then
            Ref = "RHF"
         Elseif (Iref .Eq. 1) Then
            Ref = "UHF"
         Elseif (Iref .Eq. 2) Then
            Ref = "ROHF"
         Endif 
         Epot = Etot - Ekin 
         If (Ekin .Ne. 0.0D0) Vir = Dabs(Epot/Ekin)
CSSS         Write(Iunit, 170) Ref, Method(Imethod), Etot, Vir
CSSS  170    Format(a4,1x,a11,' ENERGY =',f15.10,' VIRIAL(-V/T)  =',f13.8)
         Write(Iunit, 170) Ref, Method(1), Escf, Vir
  170 Format(a4,1x,a3,' ENERGY =',f20.12,' THE VIRIAL(-V/T)=',f13.8)
C
         Write(Iunit, "(a,a)") "----End of input file for quantum ",
     &                         "topology analysis----"
      Endif
C
      Return
      End
