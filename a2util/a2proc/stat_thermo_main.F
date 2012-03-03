      Subroutine Stat_thermo_main
c
c Statistical thermodynamic calculation for aces program system.
c
c The original program was written by J. F. Stanton 6/20/89, and
c was modified by S. Beck 10/91, remodified, audited and incorporated
c to the general interface program by Ajith Perera, 04/2006.
c
c The equations used to calculate stat mech data were taken from
c "Statistical Mechanics", Donald A. McQuarrie, HarperCollins
c Publishers, New York, NY.  Chapter 8.
c
      implicit double precision (a-h,o-z)
C
#include "mxatms.par"
#include "machsp.com"
C
      Parameter(Mxdata = 100)
C
      dimension Dmoments(3), Freq(3*Mxatms), Rot(3), Chr_rot(3), 
     &          Chr_freq(3*Mxatms), Dnormmod(9*Mxatms*Mxatms), 
     &          Iatmchrg(Mxatms), Atommass(Mxatms), 
     &          Coord(3*Mxatms), Temp(Mxdata), Press(Mxdata)
C
      logical error, Thermo_exist
      Character*80 infile,line, Fname
      Character*4 Fpgrp
C
      Data Ione, Ithree, I_inunit, I_outunit/1, 3, 15, 25/
C
C Read the temperatures and pressures at which theromodynamic parameters
C are  computed from the "THERMO" file. The thermo file is a simple 
C formatted file. The first line indicating the number of termperatue
C and presure points followed by corresponding temperatures and 
C pressures (one line for each point). The output is written to 
C standard out and the formatted "THERMO.log" file. 
C
      Call Gfname("THERMO  ", FNAME, ILENGTH)
      Inquire(FILE=FNAME(1:ILENGTH), EXIST=Thermo_exist)
      If (Thermo_exist) Then
         Open(UNIT=I_inunit, FILE=FNAME(1:ILENGTH), FORM='FORMATTED')
         Read (I_inunit, *, END=9) Ndata
         If (Ndata .GT. Mxdata) Then
             Print*, "The number of currently allowd data points", 
     &               "are:", Mxdata
             Call Errex
         Endif
         Do Idata = 1, Ndata
            Read (I_inunit, *, END=9) Temp(Idata), Press(Idata)
         Enddo
      Else
         Ndata    = 1
         Temp(1)  = 298.15
         Press(1) = 1
      Endif 
C
   9  Close (I_inunit)
      Open(UNIT=I_outunit, FILE="THERMO.log", FORM='FORMATTED') 
C
      Write(6, 910)
 910  format(8x,//,'*** STATISTICAL THERMODYNAMIC CALCULATION ***',//,
     +'All calculations based on the rigid-rotor, harmonic-oscillator,'
     +,' ideal gas',/,'approximation at given temperatures (K)',
     +'pressures (atm.).', ' The default temperature',/, 
     + 'and pressure are 298 K and 1 atm respectively.',/, 
     +'Effects of nuclear spin on the rotational ',
     +'symmetry factor are ignored.',/,'Corrections to the ',
     +'electronic partition function (assumed to be unity) must'
     +,/,'be made by the user.')
C
C Read in basic data about the molecule from the JOBARC file.
C
      Call Getrec(20, 'JOBARC', 'NREALATM', Ione, Natoms)
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Natoms*Iintfp, Coord)
      Call Getrec(20, 'JOBARC', 'ATOMMASS', Natoms*Iintfp, Atommass)
      Call Getrec(20, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
      Call Getrec(20, 'JOBARC', 'ROTCONST', Ithree*Iintfp, Rot)
      Call Getrec(20, 'JOBARC', 'LINEAR  ', Ione, Iamlinear)
      Call Getcrec(20, 'JOBARC', 'PTGP    ', 4, Fpgrp)

#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "The number of real atoms", Natoms 
      Print*, "The atomic masses:", (Atommass(I), I=1, Natoms)
      Print*, "The atomic charges:", (Iatmchrg(I), I=1, Natoms) 
      Print*, "The Cartesian Coordinates", (Coord(I), 
     &         I= 1, 3*Natoms)
      Print*, "The Rotational constants", (Rot(I), I=1, 3)
#endif
C
C Read the vibrational frequncies, mass weighted Cartesian coordiantes,
C normal modes, etc.
C
      If (Iamlinear .EQ. 1) Then
         Nvibs = 3*Natoms - 5
         Ntrro = 5
      Else
         Nvibs = 3*Natoms - 6
         Nttro = 6
      Endif
      WRITE(*,*) 'BEFORE A2RD_VIB' ! Watson
      Call A2rd_vib(Natoms, Iatmchrg, Freq, Coord, Dnormmod,
     &              Nimag_freqs, Nvibs, Ione, .FALSE.)
      WRITE(*,*) 'AFTER  A2RD_VIB'
#ifdef _DEBUG_LVLM2
      Write(6,*)
      Print*, "The vib. data from NORMCO:"
      Print*, "The Frequencies:", (Freq(I), I=1, Nvibs) 
      Print*, "The normal modes:", (Dnormmod(I), 
     &         I=1, 9*Natoms)
      Write(6,*)
#endif
C 
C Convert the rotational constants to moments of inertias 
C
      Call Rot2moments(Dmoments, Rot)
C 
      Nreal_freqs = Nvibs
      Nreal_freqs = Nreal_freqs - Nimag_freqs 
C
      Do Ivibs = 1, Nreal_freqs 
         Freq(Ivibs) = Freq(Nimag_freqs + Ivibs)
      Enddo
C    
      If (Nimag_freqs .gt. 0) write(6,915) Nimag_freqs
 915  format(/,' There are ',i3,' imaginary frequencies which will not'
     +     ,' be considered.')
C 
C Obtain the symmetry factors
C
      Call Symfac(Fpgrp, Isigma, Ntrro) 

C Calculate the characteristic vibrational temperatures
C
      Call Calccvt(Freq, Nreal_freqs, Nimag_freqs, Chr_freq)
      write(6,*)
      write(6,*),' Vibrational frequencies (cm-1):'
c
      write(6,920)(Freq(i),i=1,Nreal_freqs)
 920  format(7(2x,f8.1))
C
      write(6,*),' Vibrational temperatures (K):'
      write(6,920),(Chr_freq(i),i=1,Nreal_freqs)
C
C Calculate the characteristic rotational temperatures
C
      Call Calccrt(Dmoments, Ntrro, Chr_rot)
C
      write(6,*)
      if (ntrro.eq.5) then
        write(6,925)rot(1)
        write(6,930)rot(1)*29.97925
        write(6,935)dmoments(1)
        write(6,940)Chr_rot(1)
      else
        write(6,925)(rot(i),i=1,3)
        write(6,930)(rot(i)*29.97925,i=1,3)
        write(6,935)(dmoments(i),i=1,3)
        write(6,940)(Chr_rot(i),i=1,3)
      end if
 925  format(' Rotational constants (cm-1):',3(2x,f9.5))
 930  format('                      (GHz) :',3(2x,f9.5))
 935  format(' Moments of inertia (a.u.)  :',3(2x,f9.5))
 940  format(' Rotational temperatures (K):',3(2x,f9.5))
C
C Calculate thermodynamic properties
C
      Call Calczpe(Nreal_freqs, Zpe, Freq)
C
      write(6,945) zpe
      write(6,950) zpe/4.184
 945  format(' Zero point energy (kJ/mole)  :',2x,f9.5)
 950  format('                   (kcal/mole):',2x,f9.5)
C
C Please look at the preamble of each of these routines
C to find out what the do. In general each one of them
C do something described in McQuarrie's statistical mechanics
C book. 
C
      Weight = 0.0D0
      Do Iatm = 1, Natoms
         Weight = Weight + Atommass(Iatm)
      Enddo 
C
#ifdef _DEBUG_LVLM2
      Write(6,*) 
      Print*, "The total weight of the molecule:", weight
      Print*, "The number of Tem. and Pres. points:", Ndata
      Print*, (Temp(I), I=1, Ndata)
      Print*, (Press(I), I=1, Ndata)
#endif     
C
      Write(I_outunit, 960) 
      Write(I_outunit, 961)
 960  Format(2X, 'Temp.', 2X, 'Press.', 2X, 'Enthalpy', 2X, 'Energy',
     &       3X, 'Heat capacity', 4X, 'Entropy') 
 961  Format(4X, 'K', 5X, 'Atm', 4X, 'Kcal/mol', 2X, 'Kcal/mol', 
     &       1X, 'Kcal/mol-K', 7X, 'Kcal/mol-K')
C   
      Do Idata = 1, Ndata
         Call Calce(Nreal_freqs, Etran, Erot, Evib, Etot, 
     &              Chr_Freq, Temp(Idata), Ntrro)
C
         Call Calccv(Nreal_freqs, Cvtran, Cvrot, Cvvib, Cvtot,
     &               Chr_Freq, Temp(Idata), Ntrro)
C
         Call Calcv(Press(Idata), Temp(Idata), Vol)
C
         Call Calcs(Nreal_freqs, Stran, Srot, Svib, Stot, 
     &              Chr_Freq, Temp, Ntrro, Weight, Vol,
     &              Chr_rot, Isigma)
C
        h=etot+6.022045d23*1.380662d-23*Temp(Idata)*0.001
 
         Call Write2_stdout(Etran, Cvtran, Stran, Erot, Cvrot,
     &                      Srot, Evib, Cvvib, Svib, Etot, Cvtot,
     &                      Stot, H)
C
C The most important data that user might want to use are H, E, Cv and S
C (Enthalpy, internal energy, heat capacity and entropy). Write them to
C a file along withe temperature and pressure.
C 
      Write(I_outunit, 970) Temp(Idata), Press(Idata), H/4.184,
     &                      Etot/4.184, Cvtot/4.184, Stot/4.184
 970  Format(1X, F6.1,1X, F5.1, 3X, F7.3, 3X, F7.3, 1X, 
     &       F7.3, 11X, F7.3)
C
      Enddo
C
      Close(I_outunit)
C
      Return 
      End
C
      Subroutine Rot2moments(Dmoments, Rot)
C
c This will extract two or three numbers from string (rotational constants
c separated by spaces), translate them to dmoments of inertia, and then
c return them in the array rot.
c
c The equation for this is given by
c                h              h = 6.626176e-34 J s
c  RC(i) = -------------        c = 2.997925e8   m s-1
c          8 pi c I(i,i)
c
c The rotational constant RC(i) is read in with units of (cm-1) and the
c moment of inertia I(i,i) is in atomic units (amu ao**2 -- ao is the
c bohr radius).  The following conversion factors are necessary:
c   1 ao = 5.291787e-11 m
c   1 amu = 1.660565e-27 kg
c   1 cm-1 = 1.986478e-23 J
c to get: I(i,i) = 60.19936/RC(i).
c
c For linear molecules, only 2 rotational constants are read in (and they
c will both have the same value).  For nonlinear molecules, all 3 rotational
c constansts are read in.
C
      Implicit Double Precision (a-h,o-z)
      dimension dmoments(3),rot(3)
      logical error
      fact=60.19936
      dmoments(1)=fact/rot(1)
      dmoments(2)=fact/rot(2)
      dmoments(3)=0.0d0
      if (rot(3).ne.0.0) dmoments(3)=fact/rot(3)
      return
      end
C
      Subroutine Symfac(Ptgrp, Isigma, Ntrro)
C
c This returns the rotation factor isigma used in stat mech calculations
c in the calculation of Qrot.  In addition, it returns ntrro which is the
c number of translational and rotational degrees of freedom (5 for a
c linear molecule, 6 for a nonlinear one).
C
      implicit integer (a-z)
      character*4 ptgrp
      character*1 c1,c2,c3,c4
      c1=ptgrp(1:1)
      c2=ptgrp(2:2)
      c3=ptgrp(3:3)
      c4=ptgrp(4:4)
      zero=ichar('0')
      ntrro=6
      if (c2.eq.'X') ntrro=5
      if (c4.eq.' ') then
        isigma=ichar(c2)-zero
        if (isigma.lt.1.or.isigma.gt.9) isigma=1
      else
        isigma=(ichar(c2)-zero)*10+(ichar(c3)-zero)
      end if
      if (c1.eq.'C') then
        if (c2.eq.'X'.or.c3.eq.'s'.or.c3.eq.'i') isigma=1
      else if (c1.eq.'D') then
        isigma=2*isigma
        if (c2.eq.'X') isigma=2
      else if (c1.eq.'S')then
        isigma=isigma/2
      else if (c1.eq.'O')then
        isigma=24
      else if (c1.eq.'T')then
        isigma=12
      end if
      return
      end
C
      Subroutine Calccvt(freq, nfreq, nimag, cvt)
C
c This calculates the characteristic vibrational temperatures
c given by:
c  theta(j) = h v(j)/k        k = 1.380662e-23 J K-1
c                             h = 6.626176e-34 J s
c and v(j) is the freqency of the vibration.  The frequencies are
c in wavenumbers, so the conversion factor:  c = 2.997925e10  cm s-1
c is used to get:  theta(j)=XXX*v(j).
c
      Implicit Double Precision (a-h,o-z)
      dimension freq(nfreq),cvt(nfreq)
      factor=1.438787
C As far as I can tell, this must have been an error. When
C frequencies were read, all the imaginaries were neglected
C and Freq array starts from the real frequencies. 04/2006.
CSSS      do 10 i=1+nimag,nfreq
CSSS        cvt(i)=factor*freq(i)
CSSS 10   continue
C
      Do Ivibs = 1, Nfreq
         Cvt(Ivibs)=factor*freq(Ivibs)
      Enddo
C
      return
      end
C
      Subroutine Calccrt(dmoments, ntrro, crt)
C
c This calculates the characteristic rotational temperatures given by:
c                  h**2           k = 1.380662e-23 J K-1
c  theta(A) = --------------      h = 6.626176e-34 J s
c             8 pi**2 k I(A)
c where the dmoments of inertia I(A) are in units of (amu ao**2) which
c involves the conversion factors (1 amu=1.660565e-27 kg,
c 1 ao=5.291787e-11 m).  Using these, we get that
c   theta(A) = 86.61403/moment(A)
c For a linear molecule (ntrro=5), there is only one characteristic
c temperature calculated.  For nonlinear molecules (ntrro=6), all 3 are
c calculated.
c
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      real*8 dmoments(3),crt(3)
      factor=8.661403d1
      crt(1)=factor/dmoments(1)
      if (ntrro.eq.6) then
        crt(2)=factor/dmoments(2)
        crt(3)=factor/dmoments(3)
      end if
      return
      end

      Subroutine Calczpe(nfreq, zpe, freq)

c This calculates the zero point energy given by the equation:
c      ZPE = 0.5 h * [ v(1) + v(2) + ... + v(n) ]    h = 6.626176e-34 J s
c      c = v*lambda                                  c = 2.997925e10  cm s-1
c where ZPE is in kJ/mole and the frequencies are in wavenumbers.
c Using the conversion factors 4.359788e-18 J = 627.5057 kcal/mole and
c 1 cal = 4.184 J, you get ZPE = 5.981329e-3 * [ v(1) + ... + v(n) ]
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension freq(nfreq)
      factor=5.981329d-3
      zpe=0.0
      do 10 i=1,nfreq
        zpe=zpe+freq(i)
 10   continue
      zpe=zpe*factor
      return
      end
C
      Subroutine Calce(nfreq, etran, erot, evib, etot, cvt,
     &                 t, ntrro)
C
c This calculates the different contributions to the stat. mech.
c property E given by equations 8-25 and 8-31 in McQuarrie.  All
c are returned in kJ/mole.
c
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension cvt(nfreq)
      factor1=1.380662d-23*t
      factor2=6.022045d23*factor1/1000.0d0
      etran=1.5
      eelec=0.0D0
      if (ntrro.eq.5) then
        erot=1.0
      else
        erot=1.5
      end if
      evib=0.0
      do 10 i=1,nfreq
        t1=cvt(i)/t
        evib=evib +0.5*t1 +t1/(exp(t1)-1.0d0)
 10   continue
      etran=etran*factor2
      erot=erot*factor2
      evib=evib*factor2
      etot=etran+erot+evib+eelec
      return
      end
C
      Subroutine Calccv(nfreq, cvtran, cvrot, cvvib, cvtot, 
     &                  cvt,t,ntrro)
c
c This calculates the different contributions to the stat. mech.
c property Cv given by equations 8-26 and 8-32 in McQuarrie.  All
c are returned in J/mole-K.
c
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension cvt(nfreq)
      factor1=1.380662d-23*6.022045d23
      cvtran=1.5
      if (ntrro.eq.5) then
        cvrot=1.0
      else
        cvrot=1.5
      end if
      cvvib=0.0
      do 10 i=1,nfreq
        t1=cvt(i)/t
        t2=exp(t1)
        t3=(t2-1.0)*(t2-1.0)
        cvvib=cvvib +t1*t1*t2/t3
 10   continue
      cvtran=cvtran*factor1
      cvrot=cvrot*factor1
      cvvib=cvvib*factor1
      cvtot=cvtran+cvrot+cvvib
      return
      end
C
      Subroutine Calcv(p, t, v)
c
c This calculates the molar volume of a gas at temperature t (in K) and
c pressure p (in atm.) using:  pV=NkT.  The following conersion
c factor is used:  1 atm = 1.013247e5 N m-2.
C
      implicit double precision (a-z)
      k=1.380662d-23
      n=6.022045d23
      v=n*k*t/(p*1.013247d5)
      return
      end
C
      subroutine calcs(nfreq,stran,srot,svib,stot,cvt,t,ntrro,
     +     weight,v,crt,isigma)
c
c This calculates the different contributions to the stat. mech.
c property S given by equations 8-27 and 8-33 in McQuarrie.  All
c are returned in J/mole-K.
c
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension cvt(nfreq),crt(3)
      cn=6.022045d23
      ck=1.380662d-23
      ch=6.626176d-34
      ce=exp(1.0)
      cpi=3.141593
      dmass=weight/1000.0/cn
      sigma=real(isigma)
      factor1=ck*cn
      stran=(2.0*cpi*dmass*ck*t/ch/ch)**3
      stran=sqrt(stran)*v*sqrt(ce**5)/cn
      stran=log(stran)
      if (ntrro.eq.5) then
        srot=log(t*ce/sigma/crt(1))
      else
        srot=sqrt(cpi*ce**3)/sigma*sqrt(t**3/(crt(1)*crt(2)*crt(3)))
        srot=log(srot)
      end if
      svib=0.0
      do 10 i=1,nfreq
        t1=cvt(i)/t
        t2=exp(t1)-1.0
        t3=1.0-exp(-t1)
        svib=svib +t1/t2 -log(t3)
 10   continue
      stran=stran*factor1
      srot=srot*factor1
      svib=svib*factor1
      stot=stran+srot+svib
      return
      end
#ifdef _KEEP
      integer function islen(str)
c Compute the length of a string.
      implicit integer (a-z)
      character *(*) str
      character*1 sp,tab
      i=len(str)
      sp=' '
      tab='	'
 10   if (i.ne.0) then
        if (str(i:i).eq.sp .or. str(i:i).eq.tab) then
          i=i-1
          goto 10
        end if
      end if
      islen=i
      return
      end

      real*8 function atod(str,ERROR)
      implicit double precision (a-h,o-z)
      character*80 istr,dstr,estr
      character*1 sp,tab
      logical ERROR
      character*(*) str
      sp=' '
      tab='	'
      ilen=islen(str)
      ERROR=.false.
c Error occurs if str is a null string.
      if (ilen.eq.0) then
        ERROR=.true.
        atod=1.0d0
        return
      end if
c Strip all blanks from the string and get sign.
 10   if (str(1:1).eq.sp .or. str(1:1).eq.tab) then
        str=str(2:ilen)
        ilen=ilen-1
        goto 10
      end if
      isign=1
      if (str(1:1).eq.'+') then
        str=str(2:ilen)
        ilen=ilen-1
      else if (str(1:1).eq.'-') then
        isign=-1
        str=str(2:ilen)
        ilen=ilen-1
      end if
      if (ilen.eq.0) then
        ERROR=.true.
        atod=2.0d0
        return
      end if
 20   if (str(1:1).eq.sp .or. str(1:1).eq.tab) then
        str=str(2:ilen)
        ilen=ilen-1
        goto 20
      end if
c Split the string into integer, decimal, and exponent parts.
      iint=0
      istr=''
      idec=0
      dstr=''
      ndec=0
      iexp=0
      estr=''
      id=if_cidx(str,'.')
      ie=if_clidx(str,'dDeE')
      if ((id.gt.ie.and.ie.ne.0) .or. ie.eq.ilen) then
        ERROR=.true.
        atod=3.0d0
        return
      end if
      if (id.eq.0) then
        if (ie.eq.0) then
          istr=str
        else if (ie.eq.1) then
          istr='1'
          estr=str(2:ilen)
        else
          istr=str(1:ie-1)
          estr=str(ie+1:ilen)
        end if
      else if (id.eq.1) then
        if (ie.eq.0) then
          dstr=str(2:ilen)
          if (dstr(1:1).eq.sp .or. dstr(1:1).eq.tab) then
            ERROR=.true.
            atod=4.0d0
            return
          end if
        else if (ie.eq.2) then
          ERROR=.true.
          atod=5.0d0
          return
        else
          dstr=str(2:ie-1)
          if (dstr(1:1).eq.sp .or. dstr(1:1).eq.tab) then
            ERROR=.true.
            atod=6.0d0
            return
          end if
          estr=str(ie+1:ilen)
        end if
      else if (id.eq.ilen) then
        i=ilen-1
        istr=str(1:i)
        if (istr(i:i).eq.sp .or. istr(i:i).eq.tab) then
          ERROR=.true.
          atod=7.0d0
          return
        end if
      else
        i=id-1
        istr=str(1:i)
        if (istr(i:i).eq.sp .or. istr(i:i).eq.tab) then
          ERROR=.true.
          atod=8.0d0
          return
        end if
        if (ie.eq.0) then
          dstr=str(id+1:ilen)
          if (dstr(1:1).eq.sp .or. dstr(1:1).eq.tab) then
            ERROR=.true.
            atod=9.0d0
            return
          end if
        else if (ie.eq.id+1) then
          estr=str(ie+1:ilen)
        else
          dstr=str(id+1:ie-1)
          if (dstr(1:1).eq.sp .or. dstr(1:1).eq.tab) then
            ERROR=.true.
            atod=10.0d0
            return
          end if
          estr=str(ie+1:ilen)
        end if
      end if
      if (islen(istr).ne.0) iint=i_atoi(istr,ERROR)
      if (ERROR) then
        atod=11.0d0
        return
      end if
      ndec=islen(dstr)
      if (ndec.ne.0) idec=i_atoi(dstr,ERROR)
      if (ERROR) then
        atod=12.0d0
        return
      end if
      if (islen(estr).ne.0) iexp=i_atoi(estr,ERROR)
      if (ERROR) then
        atod=13.0d0
        return
      end if
      d=real(iint)
      if (ndec.ne.0) d=d+ real(idec)/10.0**ndec
      if (iexp.ne.0) d=d*10.0**iexp
      atod=d*real(isign)
      return
      end

      integer function i_atoi(str,ERROR)
      implicit integer (a-z)
      character*(*) str
      character*1 c,sp,tab
      logical ERROR
      i_atoi=0
      ERROR=.false.
      ilen=islen(str)
c Error occurs if str is a null string.
      if (ilen.eq.0) then
        ERROR=.true.
        i_atoi=1
        return
      end if
c isign is the sign
c i is the position of the character currently being examined
c val is the value to date
      sign=1
      i=1
      zero=ichar('0')
      sp=' '
      tab='	'
c The sign section may consist of spaces, tabs, or a single sign.
c A digit means go to the integer section.
      do while (str(i:i).eq.sp .or. str(i:i).eq.tab)
        i=i+1
      end do
      if (str(i:i).eq.'+' .or. str(i:i).eq.'-') then
        if (str(i:i).eq.'-') sign=-1
        i=i+1
      end if
      if (i.gt.ilen) then
        i_atoi=2
        ERROR=.true.
        return
      end if
      do while (str(i:i).eq.sp .or. str(i:i).eq.tab)
        i=i+1
      end do
c The integer section may consist only of digits.
 10   c=str(i:i)
      if (c.lt.'0' .or. c.gt.'9') then
        ERROR=.true.
        i_atoi=3
        return
      end if
      i_atoi=10*i_atoi + ichar(c)-zero
      i=i+1
      if (i.le.ilen) goto 10
      i_atoi=sign*i_atoi
      return
      end

      integer function if_cidx(str,char)
      implicit integer (a-z)
      character *(*) str
      character*1 char
      i=1
      strlen=islen(str)
 10   if (i.gt.strlen) then
        if_cidx=0
        return
      end if
      if (str(i:i).eq.char) then
        if_cidx=i
        return
      end if
      i=i+1
      goto 10
      end

      integer function if_clidx (str,charlist)
      implicit integer (a-z)
      character *(*) str, charlist
      strlen=islen(str)
      i=1
 10   if_clidx=if_cidx(charlist,str(i:i))
      if (if_clidx.ne.0) then
        if_clidx=i
        return
      end if
      i=i+1
      if (i.gt.strlen) return
      goto 10
      end
#endif
