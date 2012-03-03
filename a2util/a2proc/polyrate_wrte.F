C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
      Subroutine polyrate_wrte(FlgPOLYGeom,FlgPOLYElec,FlgPOLYGrad,
     &                       FlgPOLYHess,Natom,Nreal,NOcc,AtmMass,
     &                       Coord,Vgrad,Grad,Hess,Freq,Imap,Iord,
     &                       iuhf)
C
C     -----
C     interface ACESII
C
C     -----
C     called by Main
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
      Implicit double precision (a-h,o-z)
C
#include "machsp.com"
#include "mxatms.par"
C
      Logical FlgPOLYGeom,FlgPOLYElec,FlgPOLYGrad,FlgPOLYHess,
     &        FRQARC_EXIST
      Character*5  AtmLabel(Mxatms)
      character*10 Label
      Character*4  PTGRP
      Character*80 FNAME
C
      Dimension NOcc(16),AtmMass(Natom),Coord(3,Natom),VGrad(3,Nreal),
     &          Grad(3,Nreal),Imap(Nreal),Iord(Nreal),
     &          Hess(9*Natom*Natom),Freq(3*Nreal)

C     translation factor from bohr to angstrom
      Parameter (b2ang=1.889725989)
C
      Logical debug, polyrate_exist
      Parameter (debug=.True.)
      If (debug) write(*,*) '\n--- In the ACESII'
C
      IUnitO=30
C
C---------------------------------------------------------------------
C
C     Read number of atoms
C---------------------------------------------------------------------
C
      Call Getcrec(1,'JOBARC','ZSYM',5*Natom,AtmLabel)
*     --- print out ---
      If (debug) then
      Write(*,9000) Nreal
      Write(*,9100) Natom
      Write(*,9200) (AtmLabel(i),i=1,Natom)
      End if
*     -----------------
 9000 Format(4X,'--- ','number of real atoms =',I7)
 9100 Format(4X,'--- ','number of zmat atoms =',I7)
 9200 Format(4X,'--- ','atoms:','\n',8X,10A5,('\n',8X,10A5))
C
C Open the file POLYRATE, if it does exists, then close it and
C reopen it.
C
      CALL GFNAME('POLYRATE',FNAME,ILENGTH)
      Inquire(FILE=FNAME(1:ILENGTH), EXIST=Polyrate_exist)
      If (Polyrate_exist) Then
         Open(UNIT=IUnitO, FILE=FNAME(1:ILENGTH), FORM='FORMATTED')
         Close(IUnito, Status="DELETE")
         Open(UNIT=IUnitO, FILE=FNAME(1:ILENGTH), FORM='FORMATTED',
     &         STATUS='NEW')
      Endif 
C
C---------------------------------------------------------------------
C
C     Read coordinate and translate from bohr to angstrom
C---------------------------------------------------------------------
C
      If (flgPOLYGeom) then
C
         Call Getrec(20,'JOBARC','ATOMMASS',Natom*IINTFP,AtmMass)
         Call Getrec(20,'JOBARC','COORD   ',3*Natom*IINTFP,Coord(1,1))
         Do j=1,3
            Do i=1,Natom
               Coord(j,i)=Coord(j,i)/b2ang
            End do
         End do
         Open(IUnitO,File='Geom.rate',Status='Unknown')
         Do i=1,Natom
            Do j=2,5
               jchar=ichar(AtmLabel(i)(j:j))
               If (jchar .GE. 65 .AND. jchar .LE. 90) then
                  AtmLabel(i)(j:j)=char(jchar+32)
               End if
            End do
            If (AtmLabel(i)(1:5) .NE. 'X    ') then
               Write(IUnitO,7000) AtmLabel(i),(Coord(j,i),j=1,3),
     &                            AtmMass(i)
            End if
         End do
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
      If (flgPOLYElec) then
C
C        Read the number of electrons and estimate degeneracy
C        -------------------------------------------------------------
         Call Getrec(20,'JOBARC','UHFRHF  ',1,iUHF)
         Call Getrec(20,'JOBARC','NIRREP  ',1,NIrRep)
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
         If (debug) then
         Write(*,6905) iUHF
         Write(*,6910) NIrRep
         Write(*,6920) (NOcc(i),i=1,NIrRep)
         Write(*,6930) (NOcc(i+8),i=1,NIrRep)
         Write(*,6940) Nalph
         Write(*,6950) Nbeta
         Write(*,6960) Ntotal
         Write(*,6970) Mult
         End if
*        -----------------
 6905    Format(4X,'--- ','iUHF  =',I7)
 6910    Format(4X,'--- ','number of irr.  rep. =',I7)
 6920    Format(4X,'--- ','number of alpha orbs :',8I4)
 6930    Format(4X,'--- ','number of beta  orbs :',8I4)
 6940    Format(4X,'--- ','number of alpha elecs=',I7)
 6950    Format(4X,'--- ','number of beta  elecs=',I7)
 6960    Format(4X,'--- ','total number of elecs=',I7)
 6970    Format(4X,'--- ','spin multiplicity    =',I7)

C        Read total energy
C        -------------------------------------------------------------
         CALL Getrec(20,'JOBARC','TOTENERG',IINTFP,Etot)
         Open(IUnitO,File='Elec.rate',Status='Unknown')
         Write(IUnitO,6000) Mult,Etot
         Close(IUnitO)
C
      End if
 6000 Format(I3,F20.10)
C
C---------------------------------------------------------------------
C
C     Read gradient
C---------------------------------------------------------------------
C
      If (flgPOLYGrad) then
C
         Call Getrec(20,'JOBARC','MAP2ZMAT',Nreal,IMAP)
         Do i=1,Nreal
            IORD(i)=i
         End do
         Call iSort2(Nreal,IMAP,IORD,2)
C        --- IORD transforms from VMol order to ZMAT order
C
         Call Getrec(0,'JOBARC','GRADIENT',Length,VGrad(1,1))
   
          If (Length .gt. 0) Call Getrec(20,'JOBARC','GRADIENT',
     &                                   3*Nreal*IINTFP,VGrad(1,1))
         Do i=1,Nreal
            k=IORD(i)
            Do j=1,3
               Grad(j,i)=VGrad(j,k)
            End do
         End do
         Open(IUnitO,File='Grad.rate',Status='Unknown')
         Do i=1,nreal
            Write(IUnitO,5000) (Grad(j,i),j=1,3)
         End do
         Close(IUnitO)
C
      End if
 5000 Format(7X,3F18.10)
C
C---------------------------------------------------------------------
C
C     Read Hessian
C---------------------------------------------------------------------
C
      If (flgPOLYHess) then
C
         Call Getrec(20,'JOBARC','CART_HES',9*Natom*Natom*IINTFP,
     &               Hess)
         Open(IUnitO,File='Hess.rate',Status='Unknown')
         m=0
         Do i=1,Natom
            Do ix=1,3
               m=m+1
               If (AtmLabel(i)(1:5) .NE. 'X    ') then
                  n=0
                  Do j=1,i
                     Do jx=1,3
                        n=n+1
                        If (AtmLabel(j)(1:5) .NE. 'X    ') then
                           k=Natom*3*(n-1)+m
                           If (n .LE. m) Write(IUnitO,4100) Hess(k)
                        End if
                     End do
                  End do
               End if
            End do
         End do
         Close(IUnitO)
C
      End if
 4100 Format(D20.10)
C
C---------------------------------------------------------------------
C
      Call Getrec(20, 'JOBARC', 'LINEAR  ', 1, Iamlinear)
      Print*, Iamlinear, natom
      If (Iamlinear .EQ. 1) Then
         Nvibs = 3*Nreal - 5
      Else
         Nvibs = 3*Nreal - 6
      Endif
C
      INQUIRE(FILE='FRQARC',EXIST=FRQARC_EXIST)
      IF(FRQARC_EXIST)THEN
        OPEN(UNIT=5,FILE='FRQARC',FORM='UNFORMATTED',STATUS='OLD')
        rewind(5)
      else
        write(6,20)
 20     format(T3,'@RDVIB-F, The FRQARC file not present.')
        call errex
      ENDIF
C
      READ(5) LABEL, PTGRP, NRX, NIMAG, ZMAS, ZIX, ZIY, ZIZ,
     &          (Freq(IVIB), IVIB=1, Nvibs)
      
C
C Note that Frequencies written in ascending order. Since all the
C imaginary frequencies have negative sign they are at the begining.
C
      If (NIMAG .GT. 0) Then
         DO Iimag = 1, NIMAG
            Freq(Iimag) =  -Freq(Iimag)
         End Do
      End If
C
      If (debug)  then
         Write(6,*) "The vibrational Frequencies"
         Write(*,"(4(2X,F10.4))") (freq(I), I = 1, Nvibs)
      Endif
      Open(IUnitO,File='Freq.rate',Status='Unknown')
      Write(IUnitO,"(4(2X,F18.8))") (freq(I), I = 1, Nvibs) 
      Close(IUnitO)
C
      If (debug) write(*,'(A)') '--- Out of the a2rate_wrte'
      Return
C
      End
