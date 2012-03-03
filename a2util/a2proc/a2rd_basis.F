      subroutine A2rd_basis(natoms,iatchrg,gexp,coef,nshells,
     &                      iangmom,iunit,Write_Molden_file)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (mxatoms=100)
      parameter (nelement=103)
      parameter(toler=0.0000000D0)
      parameter(mxangmom=7)
      parameter(mxcoef=30)
      logical yesno, Write_Molden_file
      character*2 celeaces(nelement)
      character*1 cangmom(mxangmom)
      character*80 wrk
      character*80 wrk2
      character*80 wrk3 
      character*30 cbasnam(mxatoms)
      character*30 cbasnam2(mxatoms)
      integer iatchrg(natoms),iangmom(mxangmom),
     &  icbf(mxangmom),iexp(mxangmom),ireorder(mxatoms),
     &  lsearch(mxatoms)

      double precision gexp(mxangmom,mxcoef)
      double precision coef(mxangmom,mxcoef,mxcoef)

#include "flags.com"
      data (celeaces(i),i=1,nelement)
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
      data (cangmom(i),i=1,mxangmom)
     &  /'s','p','d','f','g','h','i'/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if(iflags(61).eq.1) then
c
c>>>>>>>>>>>>>>>>>>>>>>> BASIS SET IS NOT SPECIAL<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
        wrk = ' '
        call getrec(1,'JOBARC','BASNAMLN',1,lbasnam)
        call getcrec(1,'JOBARC','BASISNAM',lbasnam,wrk)
        do iatom=1,natoms
           i=2
           if (celeaces(iatchrg(iatom))(2:2).ne.' ') i=3
           cbasnam(iatom)=celeaces(iatchrg(iatom))
           cbasnam(iatom)(i:i)=':'
           cbasnam(iatom)(i+1:)=wrk
         end do
      else
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BASIS IS SPECIAL<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
c When BASIS=SPECIAL is choosen, the basis set entries are made at the
c bottom of ZMAT (after *ACES2 namelist), and the basis entires are made
c in the same order as the atom entries. A problem can arise if ACES II
c switches the first two atoms to get 2-1-3 connectivity pattern. We
c need to take care of that situation. Ajith Perera, 05/2005.
c
c get basis set name from ZMAT
c
        INQUIRE(FILE='ZMAT',EXIST=YESNO)
        IF(YESNO)THEN
          OPEN(UNIT=4,FILE='ZMAT',FORM='FORMATTED',STATUS='OLD')
          rewind(4)
        else
          write(6,'(T3,a)')'@RDBASIS-f, ZMAT file not present.'
          call errex
        ENDIF
c
300     READ(4,'(A)',END=900)WRK
        if (wrk(1:6).ne.'*ACES2') goto 300
c now look for a blank line
310     READ(4,'(A)',END=904)WRK
        i=1
        do while ((wrk(i:i).eq.' '.or.wrk(i:i).eq.achar(9)).and.
     &            i.le.len(wrk))
           i=i+1
        end do
        if (i.le.len(wrk)) goto 310
        do iatom=1,natoms
           read(4,'(A)',end=910)cbasnam(iatom)
        end do
        call getrec(-1,'JOBARC','12SWITCH',1,iswitch)
        if (iswitch.gt.0) then
           wrk        = cbasnam(2)
           cbasnam(2) = cbasnam(1)
           cbasnam(1) = wrk
        end if
c
100     CLOSE(UNIT=4,STATUS='KEEP')
      endif 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GENBAS<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INQUIRE(FILE='GENBAS',EXIST=YESNO)
      IF(YESNO)THEN
        OPEN(UNIT=4,FILE='GENBAS',FORM='FORMATTED',STATUS='OLD')
      else
        write(6,340)
340     format(T3,'@RDBASIS-F, GENBAS file not preset.')
        call errex
      ENDIF
c
c loop over atoms, searching for basis set 
c
      do 400 iatom=1,natoms
        rewind(4)
        wrk=cbasnam(iatom)
410     READ(4,'(A)',END=906)WRK2
        if(wrk2.ne.wrk)goto 410

        READ(4,'(A)',END=908)WRK
        read(4,'(A)',END=908)wrk
 
        read(4,*)nshells
        read(4,*)(iangmom(i),i=1,nshells)
        read(4,*)(icbf(i),i=1,nshells)
        read(4,*)(iexp(i),i=1,nshells)
c
c loop over shells reading the exponents and coefficents
c note the coefficents are in columns
c
        do 420 ishell=1,nshells
          read(4,*)
     &      (gexp(ishell,j),j=1,iexp(ishell))
          do 430 icoef=1,iexp(ishell)
            read(4,*)
     &        (coef(ishell,ibas,icoef),ibas=1,icbf(ishell)) 

430       continue
420     continue
c ----------------------------------------------------------------------
c debug
c        write(6,*)' here is what has been read in for: '
c        write(6,2000)wrk2(1:lsearch)
c        write(6,3000)nshells
c3000    format(I3)
c        write(6,3001)(iangmom(i),i=1,nshells)
c3001    format(I5)
c        write(6,3002)(icbf(i),i=1,nshells)
c3002    format(I5)
c        write(6,3003)(iexp(i),i=1,nshells)
c3003    format(I5)
c ----------------------------------------------------------------------
        If (Write_Molden_file) Then
           write(iunit,440)iatom
440        format(I3,' 0')
           do 450 ishell=1,nshells
              do 460 ibas=1,icbf(ishell)
c
c figure out how many primitives are needed to descripe this basis function
c
               ineed=0
               do 470 iprim=1,iexp(ishell)
                  if(abs(coef(ishell,ibas,iprim)).gt.toler)then
                     ineed=ineed+1
                 endif
470            continue
c
c write out the primitives and coefficents
c
               write(iunit,480)' ',Cangmom(iangmom(ishell)+1),ineed
480            format(A1,A1,3x,I2,' 1.00')
               do 490 iprim=1,iexp(ishell)
                   if(abs(coef(ishell,ibas,iprim)).gt.toler)then
                        write(iunit,500)' ',gexp(ishell,iprim),
     &                  ' ',coef(ishell,ibas,iprim)
500                     format(A,d18.10,A,d18.10)
                   endif
490            continue 
460           continue
450        continue
           write(iunit,*)
       Endif
c ----------------------------------------------------------------------
400   continue
      close(unit=4,status='KEEP')
      goto 1000 
900   WRITE(6,901)
901   FORMAT(T3,'@RDBASIS-F, *ACES2 namelist not found on ZMAT.')
      goto 911

902   WRITE(6,903)
903   FORMAT(T3,'@RDBASIS-F, keyword BASIS= not found in namelist.')
      goto 911

904   WRITE(6,905)
905   FORMAT(T3,'@RDBASIS-F, blank line is missing.  This specifies ',
     & 'the end of the namelist.')
      goto 911

906   write(6,907)
907   format('@RDBASIS-F, the following basis set not ',
     & 'present in GENBAS.')
      write(6,2000)cbasnam(iatom)
      goto 911

908   write(6,909)
909   format('@RDBASIS-F, GENBAS file terminates after ',
     & 'basis set specification for:')
      write(6,2000)cbasnam(iatom)
      goto 911

910   write(6,*)' ZMAT file terminates unexpectly'
      write(6,*)' basis sets read in: ',iatom
      write(6,*)' total basis sets needed: ',natoms
      goto 911

911   continue
      CLOSE(UNIT=4,STATUS='KEEP')
      CALL ERREX

1000  continue
      return
2000  format(A)
      end

