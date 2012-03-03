      subroutine molden_rdvib(natoms,iatchrg,freq,freqco,dnormmd,nvib,
     &                        iunit)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nelement=103)
      character *80 wrk
      character*2 celemol(nelement)
      logical yesno
      integer iatchrg(natoms)
      double precision freq(nvib),freqco(3,natoms),
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
      INQUIRE(FILE='NORMCO',EXIST=YESNO)
      IF(YESNO)THEN
        OPEN(UNIT=4,FILE='NORMCO',FORM='FORMATTED',STATUS='OLD')
        rewind(4)
      else
        write(6,10)
10      format(T3,'@RDVIB-F, NORMCO file not present.')
        call errex
      ENDIF
      read(4,*)wrk
      do 20 iatom=1,natoms
        read(4,*)(freqco(i,iatom),i=1,3)
20    continue

      do 30 ivib=1,nvib
        read(4,*)wrk
        read(4,*)freq(ivib)
        read(4,*)wrk
        do 40 iatom=1,natoms
          read(4,*)(dnormmd(ivib,i,iatom),i=1,3)
40      continue
30    continue
      close(unit=4,status='KEEP')
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      do 50 ivib=1,nvib
        write(iunit,60)' ',freq(ivib)
60      format(A,f9.4)
50    continue
      write(iunit,70)'[FR-COORD]'
70    format(A)
      do 80 iatom=1,natoms
        write(iunit,90)' ',celemol(iatchrg(iatom)),freqco(1,iatom),
     &    freqco(2,iatom),freqco(3,iatom)
90      format(A,A2,3x,f9.6,3x,f9.6,3x,f9.6)
80    continue
      write(iunit,100)'[FR-NORM-COORD]'
100   format(A)
      do 110 ivib=1,nvib
        write(iunit,115)'vibration    ',ivib
115     format(A,I2)
        do 120 iatom=1,natoms
          write(iunit,130)'   ',dnormmd(ivib,1,iatom),
     &      dnormmd(ivib,2,iatom),dnormmd(ivib,3,iatom)
130       format(A,f9.6,4x,f9.6,4x,f9.6)
120     continue 
110   continue
c/////////////////////////////////////////////////////////////////////////////
      return
      end
