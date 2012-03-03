      subroutine molden_rdgeom(natoms,iatchrg,coord,ifuchrg,fucoord,
     &                         iunit)
c------------------------------------------------------------------------------
c read coordinates from jobarc 
c------------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nelement=103)
      parameter (btoang=0.52917724924D0)
      character*2 celemol(nelement)
      dimension iatchrg(natoms),ifuchrg(natoms)
      double precision coord(3,natoms), fucoord(3,natoms)
#include "machsp.com"
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
      call izero(iatchrg,natoms)
      call zero(coord,3*natoms)
c
c read in the atomic charges and coordinates
c get rid of dummy atoms
c convert from bohr to angstrom
c
      call getrec(20,'JOBARC','ATOMCHRG',natoms,ifuchrg)
      call getrec(20,'JOBARC','COORD',3*natoms*iintfp,fucoord)
      i=1
      do j = 1, natoms
         if (ifuchrg(j).ne.0) then
c            call scopy(3,fucoord(1,j),1,coord(1,i),1)
            coord(1,i)=fucoord(1,j)
            coord(2,i)=fucoord(2,j)
            coord(3,i)=fucoord(3,j)
            iatchrg(i)=ifuchrg(j) 
            i=i+1
         end if
      end do 

      call getrec(20,'JOBARC','NREALATM',1,nrealatm)
      if(i.ne.nrealatm+1)then
        write(6,*)' Problem reading coordinates from JOBARC file'
        call errex
      endif

      natoms = nrealatm
      do i = 1, natoms
         coord(1,i)=coord(1,i)*btoang
         coord(2,i)=coord(2,i)*btoang
         coord(3,i)=coord(3,i)*btoang
      end do

      do i = 1, natoms
         write(iunit,'(A2,4X,i2,1X,i2,3X,f10.6,3X,f10.6,3X,f10.6)')
     &        celemol(iatchrg(i)),
     &        i,
     &        iatchrg(i),
     &        (coord(j,i),j=1,3)
      end do

      return
      end

