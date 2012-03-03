      Subroutine a2write_wfn(Wfn, Nangmom, Alpha, Nprim_4shl, 
     &                       Iangmom, IShell_4atm, Itype, Icent, 
     &                       Eigvals, Occnm, Coord, Maxocc, 
     &                       Nshls, Itot_prim, Nbfns,
     &                       Nreal_atoms, Iatmchrg)
      
      Implicit Double Precision (A-H, O-Z)
      Character*2 Atmlab(103)
      Character Title(10)
      Integer Angmom

      Dimension Coord(Nreal_atoms, 3), Iatmchrg(Nreal_atoms),
     &          Wfn(Itot_prim, Nbfns), Nprim_4shl(Nreal_atoms, NShls),
     &          Alpha(Itot_prim), Nangmom(Nreal_atoms), 
     &          Itype(Itot_prim), Iangmom(Nreal_atoms, Nshls), 
     &          IShell_4atm(Nreal_atoms), Icent(Itot_prim),
     &          Eigvals(Nbfns), Occnm(Nbfns)
 
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
      Data Iunit /10/
      
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6, "(a)") "The number of functions on each atom"
      Write(6, "(4(1x,I3))") (Nangmom(I), I=1, Nreal_atoms)
      Write(6,*)
      Write(6, "(a)") "The number of prim. functions on each atom"
      Do I=1, Nreal_atoms
      Write(6, "(3(1x,I3))") (NPrim_4shl(I,J), J=1, 3)
      Enddo
      Write(6,*) 
      Write(6,*) "The number of shells per atom"
      Write(6, "(4(1x,I3))") (IShell_4atm(I), I=1, Nreal_atoms)
      Write(6,*) 
      Write(6,*) "The shell angular momentum per atom s=1,p=3,d=6,,"
      Do I=1, Nreal_atoms
         Write(6, "(4(1x,I3))") (Iangmom(I, J), J=1, IShell_4atm(I))
      Enddo 
#endif
      Ioff = 0
      Do Iatms = 1, Nreal_atoms
         Do Ishl = 1, IShell_4atm(Iatms)
            Angmom = Iangmom(Iatms, Ishl)
            Nconfs = Nprim_4shl(Iatms, Ishl)
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
      Enddo
           
      Open(unit=Iunit, File="WFNANL", Form="Formatted", 
     &     Status="New")
     
      Write(Iunit, "(a,a)") "Input file for quantum ",
     &             "Topology analysis"
C
      Write(IUnit, 100) Nbfns, Itot_prim, Nreal_atoms
  100 Format('GAUSSIAN',10x,I5,' MOL ORBITALS',1x,I6,' PRIMITIVES',
     &       4x,I5,' NUCLEI')

      Do Iatms = 1, Nreal_atoms
         Write(Iunit, 110)  Atmlab(Iatmchrg(Iatms)), Itams, Iatms, 
     &                      (Coord(Iatms, J), J=1,3), 
     &                      Dble(Iatmchrg(Iatms))
      Enddo
C       
  110 Format(a4,i4,4x,'(CENTRE',i3,')',1x,3f12.8,'  CHARGE =',f5.1)
C
      Write(Iunit, 120) (Icent(I), I = 1, Itot_prim)
      Write(Iunit, 130) (Itype(I), I = 1, Itot_prim)
      Write(Iunit, 140) (Alpha(I), I = 1, Itot_prim)

  120 Format('CENTRE ASSIGNMENTS',2x,20i3)
  130 Format('TYPE ASSIGNMENTS',4x,20i3)
  140 Format('EXPONENTS',1x,1p,5e14.7)
     
      Do Iorb = 1, Maxocc
         Write(Iunit, 150) Iorb, Occnm(Iorb), Eigvals(Iorb)
         Write(Iunit, 160) (Wfn(Iprim, Iorb), Iprim =1, Itot_prim)
      Enddo 

  150 Format('MO',i3,21x,'OCC NO = ',f12.8,' ORB. ENERGY =',f13.8 )
  160 Format(1p,5e16.8)
  
      Write(Iunit, "(a)") "END DATA"
      Close(Iunit)
      Return
      End
