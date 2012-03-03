
      subroutine sb_com_symm2_ks
      implicit none

#include "machsp.com"
#include "iuhf.com"
#include "symm2_ks.com"

      double precision nucrep
      integer length,ijunk,irrep,sqroff,trioff,n,tri,sqr,i,orb,iunit
      character*80 iiiifil
      character*8 title(24)
      logical bExist
      call gfname('IIII',iiiifil,length)
      inquire(file=iiiifil,exist=bExist)
      if (.not.bExist) then
         write(*,*) '@SB_COM_SYMM2_KS: No IIII integral file.'
         call errex
      end if
      iunit = 10
      open(iunit,file=iiiifil,form='unformatted',access='sequential')
      read(iunit) title,nirrep
      rewind(iunit)
      read(iunit) title,nirrep,(numbasir(i),i=1,nirrep),nucrep
      close(iunit,status='keep')
      call getrec(0,'JOBARC','NUCREP  ',length,n)
      if (length.eq.0) call putrec(1,'JOBARC','NUCREP  ',iintfp,nucrep)
      call putrec(1,'JOBARC','NUMBASIR',nirrep,numbasir)
      sqroff=1
      trioff=1
      maxirrsqr=0
      maxirrtri=0
      orb=1
      do irrep=1,nirrep

         n=numbasir(irrep)
         irrorboff(irrep)=orb
         orb=orb+n

         sqr=n*n
         tri=n*(n+1)/2

         irrsqrlen(irrep)=sqr
         irrsqroff(irrep)=sqroff
         sqroff=sqroff+sqr
         maxirrsqr=max(maxirrsqr,sqr)

         irrtrilen(irrep)=tri
         irrtrioff(irrep)=trioff
         trioff=trioff+tri
         maxirrtri=max(maxirrtri,tri)

      end do
      irrorboff(nirrep+1)=orb
      irrsqrtot=sqroff-1
      irrtritot=trioff-1

      call getrec(-1,'JOBARC','OCCUPYA0',nirrep,occup(1,1))
      call getrec(-1,'JOBARC','OCCUPYB0',nirrep,occup(1,2))
      totocc(1)=0
      totocc(2)=0
      do i=1,nirrep
         totocc(1)=totocc(1)+occup(i,1)
         totocc(2)=totocc(2)+occup(i,2)
      end do
      if (iuhf.ne.1) totocc(2)=totocc(1)
      totocca=totocc(1)
      totoccb=totocc(2)
      return
c     end subroutine sb_com_symm2_ks
      end

