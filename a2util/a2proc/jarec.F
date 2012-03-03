
c This routine grabs a record from JOBARC and prints its contents.
c Recognized datatypes include: i/nt/eger, d/ouble, f/loat, r/eal,
c ad/ouble (2D array), ai/nt/eger (2D array)

      subroutine jarec(args,dimargs,bVerbose)
      implicit none

c ARGUMENT LIST
      integer dimargs
      character*80 args(dimargs)
      logical bVerbose

c INTERNAL VARIABLES
      double precision dTmp
      integer iLst, iTmp, i, j
      integer dim1, dim2
      character*8 szRecName

c EXTERNAL FUNCTIONS
      logical leq
      integer linblnk, atoi
      external leq, linblnk, atoi

c COMMON BLOCKS
#include "machsp.com"
#include "icore.com"
#include "istart.com"

c ----------------------------------------------------------------------

      if (dimargs.ne.3) then
         print *, '@JAREC: The jarec module requires 3 arguments'
         print '()'
         print *, '  xa2proc jarec datatype RECNAME dimension'
         print *, '  xa2proc jareq datatype RECNAME dimension (quiet)'
         call aces_exit(1)
      end if

      iLst = linblnk(args(1))
      if (bVerbose) print *, 'datatype  = ',args(1)(1:iLst)

      iLst = linblnk(args(2))
      iTmp = min(iLst,8)
      szRecName = args(2)(1:iTmp)
      if (bVerbose) print *, 'RECNAME   = ',szRecName

      iLst = linblnk(args(3))
      dim1 = atoi(args(3)(1:iLst))
      if (bVerbose) print *, 'dimension = ',dim1

c   o INTEGER VECTORS
      if (
     &    leq(args(1)(1:1),'i').or.
     &    leq(args(1)(1:3),'int').or.
     &    leq(args(1)(1:7),'integer')
     &   ) then
         if (dim1.eq.0) then
c         o read the whole thing
            call getrec(0,'JOBARC',szRecName,dim1,icore(i0))
         end if
         if (dim1.eq.1) then
            call getrec(1,'JOBARC',szRecName,1,iTmp)
            if (bVerbose) then
               print *, szRecName,' = ',iTmp
            else
               print *, iTmp
            end if
         else
            if (icrsiz.ge.dim1) then
               call getrec(1,'JOBARC',szRecName,dim1,icore(i0))
               if (bVerbose) then
                  do i = 0, dim1-1
                     print *, szRecName,'(',1+i,') = ',icore(i0+i)
                  end do
               else
                  do i = 0, dim1-1
                     print *, icore(i0+i)
                  end do
               end if
            else
               print *, '@JAREC: There is not enough memory to hold ',
     &                  'the array.'
               print *, '        requested = ',dim1,' integers'
               print *, '        available = ',icrsiz
               call aces_exit(1)
            end if
         end if
      end if

c   o DOUBLE VECTORS
      if (
     &    leq(args(1)(1:1),'d').or.
     &    leq(args(1)(1:1),'f').or.
     &    leq(args(1)(1:1),'r').or.
     &    leq(args(1)(1:6),'double').or.
     &    leq(args(1)(1:5),'float').or.
     &    leq(args(1)(1:4),'real')
     &   ) then
         if (dim1.eq.0) then
c         o read the whole thing
            call getrec(0,'JOBARC',szRecName,dim1,icore(i0))
            dim1 = dim1/iintfp
         end if
         if (dim1.eq.1) then
            call getrec(1,'JOBARC',szRecName,iintfp,dTmp)
            if (bVerbose) then
               print '(2a,e20.12)', szRecName,' = ',dTmp
            else
               print '(e20.12)', dTmp
            end if
         else
            if (icrsiz.ge.dim1*iintfp) then
               call getrec(1,'JOBARC',szRecName,iintfp*dim1,icore(i0))
               if (bVerbose) print *, szRecName,'(1:',dim1,') ='
               call dmat_prt(icore(i0),dim1,1,dim1)
            else
               print *, '@JAREC: There is not enough memory to hold ',
     &                  'the array.'
               print *, '        requested = ',dim1,' doubles'
               print *, '        available = ',icrsiz/iintfp
               call aces_exit(1)
            end if
         end if
      end if

c   o INTEGER ARRAYS
      if (
     &    leq(args(1)(1:2),'ai').or.
     &    leq(args(1)(1:4),'aint').or.
     &    leq(args(1)(1:8),'ainteger')
     &   ) then
         i = 1
         do while ((i.lt.80).and.
     &             ((args(3)(i:i).ne.',').and.
     &              (args(3)(i:i).ne.'x').and.
     &              (args(3)(i:i).ne.'X')
     &             )
     &            )
            i = i + 1
         end do
         if (i.lt.80) i = i + 1
         iLst = linblnk(args(3)(i:))
         dim2 = atoi(args(3)(i:i-1+iLst))
         if (dim2.eq.0) then
c         o read the whole thing
            call getrec(0,'JOBARC',szRecName,dim2,icore(i0))
            dim2 = dim2/dim1
         end if
         if (bVerbose) print *, 'dimension = ',dim2
         if (icrsiz.ge.dim1*dim2) then
            call getrec(1,'JOBARC',szRecName,dim1*dim2,icore(i0))
            if (bVerbose) print *, szRecName,'(1:',dim1,',1:',dim2,') ='
            do iTmp = 0, dim2-1, 7
            do i = 0, dim1-1
               print '(7(x,i10))',
     &               (icore(i0+dim1*(iTmp+j)+i),j=0,min(6,dim2-1-iTmp))
            end do
            print '(/)'
            end do
         else
            print *, '@JAREC: There is not enough memory to hold ',
     &               'the array.'
            print *, '        requested = ',dim1*dim2,' integers'
            print *, '        available = ',icrsiz
            call aces_exit(1)
         end if
      end if

c   o DOUBLE ARRAYS
      if (
     &    leq(args(1)(1:2),'ad').or.
     &    leq(args(1)(1:2),'af').or.
     &    leq(args(1)(1:2),'ar').or.
     &    leq(args(1)(1:7),'adouble').or.
     &    leq(args(1)(1:6),'afloat').or.
     &    leq(args(1)(1:5),'areal')
     &   ) then
         i = 1
         do while ((i.lt.80).and.
     &             ((args(3)(i:i).ne.',').and.
     &              (args(3)(i:i).ne.'x').and.
     &              (args(3)(i:i).ne.'X')
     &             )
     &            )
            i = i + 1
         end do
         if (i.lt.80) i = i + 1
         iLst = linblnk(args(3)(i:))
         dim2 = atoi(args(3)(i:i-1+iLst))
         if (dim2.eq.0) then
c         o read the whole thing
            call getrec(0,'JOBARC',szRecName,dim2,icore(i0))
            dim2 = dim2/dim1/iintfp
         end if
         if (bVerbose) print *, 'dimension = ',dim2
         if (icrsiz.ge.dim1*dim2*iintfp) then
            call getrec(1,'JOBARC',szRecName,dim1*dim2*iintfp,icore(i0))
            if (bVerbose) print *, szRecName,'(1:',dim1,',1:',dim2,') ='
            call dmat_prt(icore(i0),dim1,dim2,dim1)
         else
            print *, '@JAREC: There is not enough memory to hold ',
     &               'the array.'
            print *, '        requested = ',dim1*dim2,' doubles'
            print *, '        available = ',icrsiz/iintfp
            call aces_exit(1)
         end if
      end if

      return
      end

