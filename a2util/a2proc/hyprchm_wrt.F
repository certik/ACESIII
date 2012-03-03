      subroutine wrt_hyprchm(natoms,coords,freqs,feigs,intens,anum,
     .                       tbohr, hfcrap, linear, iord, itag)
c
c     natoms = the total number of atoms in the calculation (INTEGER)
c
c     coords = a (natoms,3) array of cartesian coordinates (REAL)
c
c     freqs = the array that contains the eigenvalues of the Hessian (REAL)
c
c     feigs = the array that contains the eigenvectors of the Hessian (REAL)
c             (it is presumed the null vectors have been removed)
c
c     intens = the array that contains the IR intensities (REAL)
c
c     anum = the array that contains the atomic numbers (INTEGER)
c
c     atype = the array that contains the atomic symbols (CHARACTER)
c
c     itag = the array that contains the connectivity table (INTEGER)
c            (itag(natoms,10) designates the number of atoms bonded to
c             atom "x")
c
c     hfcrap = .true. if a frequency run was performed
c
c     linear = .true. if the molecular system is linear
c
c
c------ Warnings
c     The character array "line" in "subroutine entscr", dimensioned
c     currently to 10000, limits that array to a little more than 300
c     atomic centers. This is the only known hard-wired size limitation.
c------ Warnings
c

      implicit double precision (a-h,o-z)
      double precision intens(*)
      integer anum(*)
      character*2 atype(111)
      dimension coords(3,natoms),freqs(*),feigs(9*natoms*natoms),
     .          iord(*), itag(natoms,10)
      logical  tbohr, hfcrap, linear

      data bohr / 0.52917724924d0 /

      data atype / ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F',
     .             'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar',
     .             ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co',
     .             'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
     .             'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
     .             'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe',
     .             'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
     .             'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
     .             'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
     .             'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     .             'Pa', ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
     .             'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
     .             'Mt', 'Ds', 'Rg' /

c
c     if the supplied spatial coordinates are in units of
c     Bohr, then logical tbohr = .true., else tbohr = .false.
c
      if (tbohr) then
         do 100 i=1, natoms
            coords(1,i) = coords(1,i) * bohr
            coords(2,i) = coords(2,i) * bohr
            coords(3,i) = coords(3,i) * bohr
 100     continue
      endif

c
c     Generate the HYPERCHEM.hin file
c

      call hinfl(coords,atype,anum,natoms,itag)

c
c     logical hfcrap is set .true. if a frequency run was performed
c     and one wants the HyperChem "name.ent" and "name.scr" files
c     to be generated
c
      if (hfcrap) then
        call entscr(coords,atype,anum,natoms,itag,freqs,feigs,
     .              intens,linear,iord)
      endif

      return
      end




      subroutine entscr(coords,atype,anum,natoms,itag,freqs,feigs,
     .                  intens,linear,iord)
c
c     Subroutine entscr generates the "HYPERCHEN.ent" and "HYPERCHEM.scr"
c     files which HyperChem uses to display the vibrational modes
c
      implicit double precision (a-h,o-z)
      double precision intens(*)
      character*10000 line
      character*2 atype(*)
      integer anum(*)
      logical linear
      dimension coords(3,natoms), feigs(9*natoms*natoms), freqs(*), 
     .          itag(natoms,10), iord(*)
c
c     Generate the "ent" file
c
      open(unit=10,status='unknown',form='formatted',
     .     file='HYPERCHEM.ent')

      write(10,'(a)') 'HEADER'
      write(10,'(a)') 'REMARK From an Aces2 run'

      do 100 i=1, natoms
         write(10,'(a,i4,1x,a,10x,a,4x,3(1x,f7.3),23x,a)') 'HETATM ',i,
     .        atype(anum(i)),'1', (coords(j,i), j=1,3), atype(anum(i))
 100  continue

      do 110 i=1, natoms
         write(10,'(a,i4,9i5)')'CONECT ',i,(itag(i,j),j=1,itag(i,10))
 110  continue

      write(10,'(a)') 'END'
      write(10,'(a)') ' '
      close(unit=10,status='keep')

c
c     Generate the "scr" file
c

      open(unit=10,status='unknown',form='formatted',
     .     file='HYPERCHEM.scr')
c
c     If the system is linear then "linear = .true."
c
      if (linear) then
         nvec = natoms*3 - 5
      else
         nvec = natoms*3 - 6
      endif

      write(10,'(a)') 'file-format pdb'
      write(10,'(a)') 'open-file HYPERCHEM.ent'
      write(10,'(a,i5)') 'ir-band-count =  ', nvec


      do 150 i=1, nvec
         istr = 1
         istp = 10

         do 140 j=1, 3*natoms
            write(unit=line(istr:istp),fmt='(f9.6,a)') 
     .                 feigs((iord(i)-1)*3*natoms + j), ','
            istr = istr + 10
            istp = istp + 10
 140     continue

         write(10,'(a,i3,a,f12.4)') 'ir-frequency(',i,') = ',
     .        freqs(i)
         write(10,'(a,i3,a,f12.4)') 'ir-intensity(',i,') = ',
     .        intens(iord(i))
         write(10,'(a,i3,a,a)') 'ir-normal-mode(',i,') = ', 
     .        line(1:istp-11)
 150  continue

      write(10,'(a)') 'menu-compute-vibrational-spectrum'
      write(10,'(a)') 'exit-script'
      write(10,'(a)') ' '
      close(unit=10,status='keep')

      return
      end



      subroutine hinfl(coords,atype,anum,natoms,itag)
c
c     Subroutine hinfl generates the HyperChem "hin" file
c
      implicit double precision (a-h,o-z)
      character*2 atype(*)
      integer anum(*)
      dimension coords(3,natoms), radii(111), itag(natoms,10)
c
c     It is presumed no atom will have more than 9 valence bonds
c

c
c     The atomic covalent radii in A. 
c
c     From tabulations and averages of C(sp3)-X distances in Allen,
c     F.H., Kennard, O., Watson, D.G., Brammer, L., Orpen, A.G., &
c     Taylor, R. (1987) J.Chem. Soc. Perkin II, p. S1, subtracting
c     0.767 A for the radius of carbon. 
c
c     A value of 0.0d0 implies no value was found
c

      data radii /

c          H       He      Li      Be      B       C       N       O   
     . 0.299d0,0.000d0,1.562d0,1.060d0,0.830d0,0.767d0,0.702d0,0.659d0,
c          F       Ne      Na      Mg      Al      Si      P       S   
     . 0.619d0,0.000d0,1.911d0,1.602d0,1.180d0,1.090d0,1.088d0,1.052d0,
c          Cl      Ar      K       Ca      Sc      Ti      V       Cr  
     . 1.023d0,0.000d0,2.376d0,1.974d0,1.641d0,1.462d0,1.346d0,1.282d0,
c          Mn      Fe      Co      Ni      Cu      Zn      Ga      Ge  
     . 1.264d0,1.274d0,1.252d0,1.246d0,1.278d0,1.394d0,1.250d0,1.220d0,
c          As      Se      Br      Kr      Rb      Sr      Y       Zr  
     . 1.196d0,1.203d0,1.199d0,0.000d0,2.546d0,2.151d0,1.801d0,1.602d0,
c          Nb      Mo      Tc      Ru      Rh      Pd      Ag      Cd  
     . 1.468d0,1.400d0,1.360d0,1.339d0,1.345d0,1.376d0,1.445d0,1.568d0,
c          In      Sn      Sb      Te      I       Xe      Cs      Ba  
     . 1.410d0,1.390d0,1.370d0,1.391d0,1.395d0,0.000d0,2.731d0,2.243d0,
c          La      Ce      Pr      Nd      Pm      Sm      Eu      Gd  
     . 1.877d0,1.825d0,1.828d0,1.821d0,1.810d0,1.802d0,2.042d0,1.802d0,
c          Tb      Dy      Ho      Er      Tm      Yb      Lu      Hf  
     . 1.782d0,1.773d0,1.766d0,1.757d0,1.746d0,1.740d0,1.734d0,1.580d0,
c          Ta      W       Re      Os      Ir      Pt      Au      Hg  
     . 1.467d0,1.408d0,1.375d0,1.353d0,1.357d0,1.387d0,1.442d0,1.573d0,
c          Tl      Pb      Bi      Po      At      Rn      Fr      Ra  
     . 1.716d0,1.750d0,1.700d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,
c          Ac      Th      Pa      U       Np      Pu      Am      Cm  
     . 0.000d0,1.798d0,0.000d0,1.560d0,0.000d0,0.000d0,0.000d0,0.000d0,
c          Bk      Cf      Es      Fm      Md      No      Lr      Rf  
     . 0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,
c          Db      Sg      Bh      Hs      Mt      Ds      Rg  
     . 0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0,0.000d0 /

c
c     Initialize itag (connection matrix)
c
      call izero(itag,natoms*10)
c
c     Generate the connection table
c

      do 100 i=1, natoms-1
         do 200 j=i+1, natoms

            rab = r12(i,j,coords,natoms)
            rcv = rchk(i,j,radii,anum)

            if (rab .lt. rcv) then
               itag(i,10) = itag(i,10) + 1
               itag(j,10) = itag(j,10) + 1
               itag(i,itag(i,10)) = j
               itag(j,itag(j,10)) = i
            endif

 200     continue
 100  continue

      do 300 i=1, natoms
         if(itag(i,1) .eq. 0) then
            itag(i,1) = icls(i,coords,natoms)
            itag(i,10) = itag(i,10) + 1
            itag(itag(i,1),10) = itag(itag(i,1),10) + 1
            itag(itag(i,1), itag(itag(i,1),10)) = i
         endif
 300  continue

      call hinout(natoms,itag,anum,atype,coords)

      return
      end


      subroutine hinout(natoms,itag,anum,atype,coords)
c
c     Subroutine hinout writes out the "hin" file
c
      implicit double precision (a-h,o-z)
      dimension itag(natoms,10), coords(3,natoms)
      integer anum(*)
      character*2 atype(*)
      character*80 line

      open(unit=10,form='formatted',status='unknown',
     .     file='HYPERCHEM.hin')
      rewind 10

      write(10,'(a)') 'forcefield mm+'
      write(10,'(a)') 'sys 0'
      write(10,'(a)') 'view 40 0.25458 55 15 0.966211 0.2538206 0.044848
     .8 0.2456284 -0.9594546 0.1382519 0.07812156 -0.1225643 -0.9893811 
     .0.073402 -0.41257 -55.191'
      write(10,'(a)') 'seed -1111'
      write(10,'(a)') 'mol 1'

      do 100 i=1, natoms
         istr=1
         istp=6
         do 110 j=1, itag(i,10)
           write(unit=line(istr:istp),fmt='(a,i3,a)') ' ',itag(i,j),' s'
           istr = istr + 6
           istp = istp + 6
 110     continue

         istp = istp - 6

         write(10,'(a,i3,1x,a,1x,a,a,a,a,3(1x,f8.4),i3,a)') 'atom ',i,
     .        atype(anum(i)),atype(anum(i)),' ** ','h ','0 ',
     .        (coords(j,i), j=1,3), itag(i,10), line(1:istp)
 100  continue

      write(10,'(a)') 'endmol 1'
      write(10,'(a)') ' '
      close (unit=10,status='keep')

      return
      end


      integer function icls(a1,coords,natoms)
c
c     Function icls returns the nearest neighbor of center a1
c
      implicit double precision (a-h,o-z)
      integer a1, a2
      dimension coords(3,natoms)

      rsave = 1.0d10

      do 100 a2=1, natoms
         if(a2 .eq. a1) goto 100
         rtest = r12(a1,a2,coords,natoms)
         if (rsave .gt. rtest ) then
            rsave = rtest
            isave = a2
         endif         
 100  continue

      icls = isave

      return
      end


      double precision function rchk(a1,a2,radii,anum)
c
c     Function rchk returns the scaled summed covalent
c     radii of atoms a1 and a2
c
      implicit double precision (a-h,o-z)
      integer anum(*), a1, a2
      dimension radii(*)

      scal = 1.15d0
      rchk = scal * (radii(anum(a1)) + radii(anum(a2)))

      return
      end



      double precision function r12(a1,a2,coords,natoms)
c
c     Function returns the inter-nuclear distance between
c     centers a1 and a2
c
      implicit double precision (a-h,o-z)
      integer a1, a2
      dimension coords(3,natoms)

      x1 = coords(1,a1)
      y1 = coords(2,a1)
      z1 = coords(3,a1)
      x2 = coords(1,a2)
      y2 = coords(2,a2)
      z2 = coords(3,a2)

      r12 = dsqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

      return 
      end

