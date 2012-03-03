      subroutine a2rd_vib(natoms,iatchrg,freq,freqco,dnormmd,Nimag,
     &                    Nvib,iunit,Write_molden_file)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nelement=103)
      parameter (angtob=1.88972598772D0)
      character *80 wrk
      character *10 Label
      Character  *4 PTGRP
      character*2 celemol(nelement)
      logical yesno,Write_molden_file
      inTEger iatchrg(natoms)
      DOuble precision freq(nvib),freqco(3,natoms),
     &  dnormmd(nvib,3,natoms)

#include "machsp.com"
#include "flags.com"
      data (celemol(i),i=1,nelement)
     & /' H','He',
     & 'Li','Be',' B',' C',' N',' O',' F','Ne',
     & 'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     & ' K','Ca',
     & 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     & 'Ga','Ge','As','Se','Br','Kr',
     & 'Rb','Sr',
     & ' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     & 'In','Sn','Sb','Te',' I','Xe',
     & 'Cs','Ba','La',
     & 'Ce','Pr','Nd','Pm','Sm','Eu','Gd',
     & 'Tb','Dy','Ho','Er','Tm','Yb','Lu',
     & 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     & 'Tl','Pb','Bi','Po','At','Rn',
     & 'Fr','Ra','Ac',
     & 'Th','Pa',' U','Np','Pu','Am','Cm',
     & 'Bk','Cf','Es','Fm','Md','No','Lr'/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C If Writing the MOLDEN files then read the NORMCO file otherwise
C read the FRQARC file.
C
      If (Write_molden_file) Then
C
C Reset the number of frequncies to include rotational and
C vibrational modes since the NORMCO file contains all of
C them.
C
        Nvib  = 3*Natoms
        Nimag = 0
        INQUIRE(FILE='NORMCO',EXIST=YESNO)
        IF(YESNO)THEN
          OPEN(UNIT=4,FILE='NORMCO',FORM='FORMATTED',STATUS='OLD')
          rewind(4)
        else
          write(6,10)
10        format(T3,'@RDVIB-F, NORMCO file not present.')
          call errex
        ENDIF
C
C This is an error, molden visualization does not need mass weigthed 
C coordnates. Ajth Perera, 05/2010
        read(4,*)wrk
        do iatom=1,natoms
CSSS           read(4,*)(freqco(i,iatom),i=1,3)
           read(4,*) tmp
        end do 
      
        do 30 ivib=1,nvib
           read(4,*)wrk      
           read(4,'(a)') wrk(1:21)
           read(unit=wrk(1:20),fmt=*) freqy
           if (wrk(21:21) .eq. 'i') Then
              freqy = -freqy 
              Nimag = Nimag + 1
           endif
           freq(ivib) = freqy      
           read(4,*)wrk
           do 40 iatom=1,natoms
              read(4,*)(dnormmd(ivib,i,iatom),i=1,3)
#ifdef _DEBUG_LVLM1
              Write(6, "(F10.5)") (dnormmd(ivib,i,iatom),i=1,3)
#endif
         
   40      continue
   30   continue
C
        close(unit=4,status='KEEP')    

        do 50 ivib=1,nvib
           if (freq(ivib) .lt. 0.0D0) freq(ivib)= -freq(ivib)
           write(iunit,60)' ',freq(ivib)
60         format(A,f11.4)
50      continue
        write(iunit,70)'[FR-COORD] Bohr'
70      format(A)
C
C Watson: change to bohr for molden
C
        do 80 iatom=1,natoms
           write(iunit,90)' ',celemol(iatchrg(iatom)),
     &                        angtob * freqco(1,iatom),
     &                        angtob * freqco(2,iatom),
     &                        angtob * freqco(3,iatom)
90         format(A,A2,3x,f9.6,3x,f9.6,3x,f9.6)
80      continue
        write(iunit,100)'[FR-NORM-COORD]'
100     format(A)
        do 110 ivib=1,nvib
           write(iunit,115)'vibration    ',ivib
115        format(A,I2)
           do 120 iatom=1,natoms
              write(iunit,130)'   ',dnormmd(ivib,1,iatom),
     &        dnormmd(ivib,2,iatom),dnormmd(ivib,3,iatom)
130           format(A,f9.6,4x,f9.6,4x,f9.6)
120        continue 
110     continue
C
      Else ! If Write_Molden_File
C
        INQUIRE(FILE='FRQARC',EXIST=YESNO)
        IF(YESNO)THEN
          OPEN(UNIT=15,FILE='FRQARC',FORM='UNFORMATTED',STATUS='OLD')
          rewind(15)
        else
          write(6,20)
20        format(T3,'@RDVIB-F, The FRQARC file not present.')
          call errex
        ENDIF
        READ(15) LABEL, PTGRP, NRX, NIMAG, ZMAS, ZIX, ZIY, ZIZ,
     &          (Freq(IVIB), IVIB=1, Nvib)
        Nimag = 0
#ifdef _DEBUG_LVL0
        Print*, PTGRP, NRX, NIMAG, ZMAS, ZIX, ZIY, ZIZ,
     &         (Freq(IVIB), IVIB=1, Nvib)
#endif

C
C Scale coordinates to bohr from angstrom
C
C        call dscal(3*natoms,angtob,coord,1)
C Watson, this seems to already have been done and put in JOBARC
C         COORD isn't even sent over!
C
        close(unit=15,status='KEEP')
C
      Endif
C
      Return
      End
