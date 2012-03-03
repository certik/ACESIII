C---------------------------------------------------------------------
C
      SubroutinE a2read_jarc(FlgACESGeom,FlgACESElec,FlgACESGrad,
     &                       FlgACESHess,FlgACESForc,Natom,Nreal,
     &                       Imap, Coord, Etot, AtmMass, AtmLabel, 
     &                       Vgrad, Grad, Vhess, Hess, Vomega, 
     &                       Omega, IUhf,Nirrep, Nvibs, B2ang, 
     &                       Au2Invcm, Wrt_extrnl)
C
      Implicit double precision (a-h,o-z)
C
#include "machsp.com"
#include "mxatms.par"
C
      Logical FlgACESGeom,FlgACESElec,FlgACESGrad,FlgACESHess,
     &        FlgACESForc, Wrt_extrnl
      Character*5  AtmLabel(Mxatms)
      character*10 Label
      Character*4  PTGRP
C
      Dimension NOcc(16),AtmMass(Natom),Coord(3,Natom),
     &          VGrad(3,Natom),Grad(3,Nreal),Imap(Nreal),
     &          Hess(3*Nreal,3*Nreal),Freq(3*Nreal),
     &          Vhess(9*Natom*Natom),VOmega(3*Natom),
     &          Omega(Nvibs)

C
      IUnitO=30
C
C---------------------------------------------------------------------
C
C     Read number of atoms
C---------------------------------------------------------------------
C
#ifdef _DEBUG_LVLM1
      Write(6, *) 
      Write(6,9000) Nreal
      Write(6,9100) Natom
      Write(6,9200) (AtmLabel(i),i=1,Natom)
*     -----------------
 9000 Format(4X,'--- ','number of real atoms =',I7)
 9100 Format(4X,'--- ','number of zmat atoms =',I7)
 9200 Format(4X,'--- ','atoms:',8X,10A5,('\n',8X,10A5))
#endif
C
C---------------------------------------------------------------------
C
C     Read coordinate and translate from bohr to angstrom
C---------------------------------------------------------------------
C
      ierr =  1
      If (flgACESGeom) then
C
         Do j=1,3
            Do i=1,Natom
               Coord(j,i)=Coord(j,i)/b2ang
            End do
         End do
         If (Wrt_extrnl) Then
             Open(IUnitO,File='Geom.rate',Status='Unknown')
             Do i=1,Natom
                Do j=2,5
                   jchar=ichar(AtmLabel(i)(j:j))
                   If (jchar .GE. 65 .AND. jchar .LE. 90) then
                       AtmLabel(i)(j:j)=char(jchar+32)
                   End if
                End do
                If (AtmLabel(i)(1:5) .NE. 'X    ') then
                    ierr = 0
                    Do while (ierr .eq. 0) 
                       Read(IUnitO,7000, iostat=ierr) Temp_buff
                    Enddo
                    Write(IUnitO,7000) AtmLabel(i),(Coord(j,i),j=1,3),
     &                                 AtmMass(i)
                End if
             End do
             Write(IUnitO,"") 
         Endif
         Close(IUnitO)
C
      End if
 7000 Format(5X,A5,F15.10,2F18.10,F18.10)
C
C---------------------------------------------------------------------
C
C     Read total energy and degeneracy
C---------------------------------------------------------------------
C
      If (flgACESElec) then
C
C        Read the number of electrons and estimate degeneracy
C        -------------------------------------------------------------
         Call Getrec(20,'JOBARC','OCCUPYA0',NIrRep,NOcc(1)) 
         Nalph=0
         If (iUHF .EQ. 0) then
            Do i=1,NIrRep
               Nalph=Nalph+NOcc(i)
            End do
            Ntotal=Nalph*2
            Mult=1
            Nbeta = Nalph 
            call icopy(8, nocc(1), 1, nocc(9), 1)
         Else if (iUHF .EQ. 1) then
            Call Getrec(20,'JOBARC','OCCUPYB0',NIrRep,NOcc(9))
            Nbeta=0
            Do i=1,NIrRep
               Nalph=Nalph+NOcc(i)
               Nbeta=Nbeta+NOcc(i+8)
            End do
            Ntotal=Nalph+Nbeta
            Mult=(Nalph-Nbeta)+1
         End if
C
C        Calculate degeneracy in the atomic case
C        -------------------------------------------------------------
C
*        --- print out ---
#ifdef _DEBUG_LVLM1
         Write(6,6905) iUHF
         Write(6,6910) NIrRep
         Write(6,6920) (NOcc(i),i=1,NIrRep)
         Write(6,6930) (NOcc(i+8),i=1,NIrRep)
         Write(6,6940) Nalph
         Write(6,6950) Nbeta
         Write(6,6960) Ntotal
         Write(6,6970) Mult
*        -----------------
 6905    Format(4X,'--- ','iUHF  =',I7)
 6910    Format(4X,'--- ','number of irr.  rep. =',I7)
 6920    Format(4X,'--- ','number of alpha orbs :',8I4)
 6930    Format(4X,'--- ','number of beta  orbs :',8I4)
 6940    Format(4X,'--- ','number of alpha elecs=',I7)
 6950    Format(4X,'--- ','number of beta  elecs=',I7)
 6960    Format(4X,'--- ','total number of elecs=',I7)
 6970    Format(4X,'--- ','spin multiplicity    =',I7)
#endif 
C        Read total energy
C        -------------------------------------------------------------
         CALL Getrec(20,'JOBARC','TOTENERG',IINTFP,Etot)
         If (Wrt_extrnl) Then
            Open(IUnitO,File='Elec.rate',Status='Unknown')
            ierr = 0
            Do while (ierr .eq. 0)
                Read(IUnitO,7000, iostat=ierr) Temp_buff
            Enddo
            Write(IUnitO,6000) Mult,Etot
            Write(IUnitO, "")
            Close(IUnitO)
         Endif
C
      End if
 6000 Format(I3,F20.10)
C
C---------------------------------------------------------------------
C
C     Read gradient
C---------------------------------------------------------------------
C
      If (flgACESGrad) then
C
C   --- IMAP transforms from VMol order to ZMAT order
C
         Call Getrec(20,'JOBARC','GRADIENT',3*Nreal*IINTFP,
     &               VGrad(1,1))
         Do i=1,Nreal
            k=Imap(i)
            Do j=1,3
               Grad(j,k)=VGrad(j,i)
            End do
         End do
         If (Wrt_extrnl) Then
            Open(IUnitO,File='Grad.rate',Status='Unknown')
            ierr = 0
            Do while (ierr .eq. 0)
               Read(IUnitO,7000, iostat=ierr) Temp_buff
            Enddo
            Do i=1,nreal
               Write(IUnitO,5000) (Grad(j,i),j=1,3)
            End do
            Write(IUnitO,"")
            Close(IUnitO)
         Endif
C
      End if
 5000 Format(7X,3F18.10)
C
C---------------------------------------------------------------------
C
C     Read Hessian
C---------------------------------------------------------------------
C
      If (flgACESHess) then
C
         Call Getrec(20,'JOBARC','CART_HES',9*Natom*Natom*IINTFP,
     &               VHess)
         Call output(Vhess,1,3*natom,1,3*natom,3*natom,3*natom,1)

c
c-----Transform to the ZMAT order

         Ioff = 1
         Do Ivmol = 1, Natom
            Iatmzmat = Imap(Ivmol)
            Do Ixyz = 1, 3
               Icol = Ixyz + (Iatmzmat - 1)*3
               Do Jvmol = 1, Natom
                  Jatmzmat = Imap(Jvmol)
                  Irow = 1 + (Jatmzmat - 1)*3
                  Call Blkcpy(Vhess(ioff), 3, 1, Hess, 3*Natom,
     &                        3*Natom, Irow, Icol)
                  Ioff = Ioff + 3
               Enddo
            Enddo
         Enddo
C
C-----Eleminate the entry for dummy atoms
C
         Ix = 0
         Jx = 1
         Do Iatm = 1, 3*Natom
            Do Jatm = 1, 3*Natom
               Ink = 1 + (Iatm - 0.1)*3
               Jnk = 1 + (Jatm - 0.1)*3
               If ((AtmLabel(Ink)(1:5) .NE. 'X    ') .AND.
     &            ( AtmLabel(Jnk)(1:5) .NE. 'X    ')) Then
                   Ix = Ix + 1
                   Hess(Ix, Jx) = Hess(Iatm, Jatm)
                   If (Ix .EQ. 3*Natom) Ix = 0
                   If (Ix .EQ. 0) Jx = Jx + 1
               Endif
           Enddo
         Enddo 
C
         Call output(hess,1,3*natom,1,3*natom,3*natom,3*natom,1)
C
         If (Wrt_extrnl) Then
            Open(IUnitO,File='Hess.rate',Status='Unknown')
C
            ierr = 0
            Do while (ierr .eq. 0)
               Read(IUnitO,7000, iostat=ierr) Temp_buff
            Enddo
C
            Ioff = 0
            Joff = 0
            Do Iatm = 1, Nreal
                  Ioff = Ioff + (Iatm - 1)*3
                  Do Jatm = 1, Iatm
                     Joff = Joff + (Jatm - 1)*3
                     Write(IUnitO,4100) ((Hess(Ioff, Joff), 
     &                                    Ioff = 1, 3), Joff=1,3)
                  Enddo
            Enddo
            Write(IUnitO, "")
            Close(IUnitO)
         Endif
C
      Endif
C
 4100 Format(7X, 3D20.10)
C
      If (FlgACESForc) Call Getrec(20, "JOBARC", 'FORCECON', 
     &                             3*Natom*IINTFP, VOmega)
C
      Ivib = 0
      Do Imode = 1, 3*Natom
C
         If (Vomega(Imode) .lt. 0.0D0) Then
            Vomega(Imode) = - DSQRT(DABS(Vomega(Imode)))*Au2Invcm
         Else
            Vomega(Imode) = DSQRT(Vomega(Imode))*Au2Invcm
         Endif
C
         If (DABS(Vomega(Imode)) .GT. 1.0D0) Then
            Ivib = Ivib + 1
            Omega(Ivib) = Vomega(Imode)
         Endif
C
      Enddo
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The vibrational Frequencies"
      Write(*,"(4(2X,F10.4))") (Omega(I), I = 1, Nvibs)
#endif
C
      Open(IUnitO,File='Freq.rate',Status='Unknown')
      ierr = 0
      Do while (ierr .eq. 0)
         Read(IUnitO,7000, iostat=ierr) Temp_buff
      Enddo
C
      Write(IUnitO,"(4(2X,F18.8))") (Omega(I), I = 1, Nvibs)
      Write(IUnitO, "")
      Close(IUnitO)
C
      Return
      End
