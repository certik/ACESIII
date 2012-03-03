      subroutine hyprchm_rd_geomvib(natoms, iatchrg, coord,ifuchrg,
     &                            fucoord,freqco,freq,dnormmd, 
     &                            intensites,iord,normco_exsist,linear)
c------------------------------------------------------------------------------
c read coordinates from jobarc 
c------------------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      dimension iatchrg(natoms),ifuchrg(natoms),iord(3*natoms)
      character*1 askip
      character*80 line
      double precision coord(3,natoms), fucoord(3,natoms),
     &                 freqco(natoms, 3), dnormmd(9*natoms*natoms),
     &                 intensites(3*natoms), freq(3*natoms)
      logical normco_exsist, linear, tbohr
      data IONE /1/
C-----------------------------------------------------------------------------
#include "machsp.com"
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call izero(iatchrg,natoms)
      call zero(coord,3*natoms)
c
c read in the atomic charges and coordinates and remove dummy atoms
c convert from bohr to angstrom
c
      call getrec(20,'JOBARC','ATOMCHRG',natoms,ifuchrg)
      call getrec(20,'JOBARC','COORD',3*natoms*iintfp,fucoord)
C
      i=1
      do j = 1, natoms
         if (ifuchrg(j).ne.0) then
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
C
      natoms = nrealatm

      inquire(file='NORMCO',exist=normco_exsist)

      if(normco_exsist)then
        open(unit=4,file='NORMCO',form='formatted',status='old')
        rewind(4)
      else
         goto 5
      endif
C
      call getrec(20, 'JOBARC', 'LINEAR  ',IONE, ILINEAR)
      linear = .false.
      IF (ILINEAR .EQ. 1) linear = .TRUE.
C
C Read the mass weigthed Cartesian coordinates
C
      read(4,*) askip
      do 20 iatom=1,natoms
        read(4,*)(freqco(i,iatom),i=1,3)
20    continue

c      
c Read vibrational frequncies and normal coordinates
c 

      do 30 i=1,3*natoms
         read(4,'(a)') askip
         read(4,'(a)') line(1:21)
         read(unit=line(1:20),fmt=*) freqy
         if (line(21:21) .eq. 'i') freqy = -freqy
         freq(i) = freqy
         read(4,*) askip

         do 40 j=1,natoms
            read(4,*)(dnormmd((i-1)*3*natoms+(j-1)*3+k), k=1,3)
 40      continue

 30   continue

c
c     Identify the null vectors
c

      do 44 i=1, natoms*3
         iord(i) = i
 44   continue

      do 45 i=1, 3*natoms - 1
         do 46 j=i+1,3*natoms
            if(dabs(freq(i)) .lt. dabs(freq(j))) then
               save = freq(i)
               freq(i) = freq(j)
               freq(j) = save
               isave = iord(i)
               iord(i) = iord(j)
               iord(j) = isave
            endif
 46      continue
 45   continue

      if (linear) then
         istp = 3*natoms - 5
      else
         istp = 3*natoms - 6
      endif


      do 47 i=1, istp - 1
         do 48 j=i+1, istp
            if(freq(j) .lt. freq(i)) then
               save = freq(i)
               freq(i) = freq(j)
               freq(j) = save
               isave = iord(i)
               iord(i) = iord(j)
               iord(j) = isave
            endif
 48      continue
 47   continue

c
c Read intensities
c
      do 50 i = 1, 3*natoms
         read(4, *)  askip
         read(4, *) intensites(i)
50    continue

      close(unit=4,status='KEEP')

 5    return
      end

