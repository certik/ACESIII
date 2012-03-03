
c This routine fills JOBARC records with zeroes.

      subroutine zerorec(args,dimargs)
      implicit none

c ARGUMENT LIST
      integer dimargs
      character*80 args(dimargs)

c INTERNAL VARIABLES
      integer iArg
      integer nInts

c ----------------------------------------------------------------------

      do iArg = 1, dimargs
         call getrec(0,'JOBARC',args(iArg)(1:8),nInts,0)
         call putrec(0,'JOBARC',args(iArg)(1:8),nInts,0)
         print *, '@ZEROREC: zeroed ',args(iArg)(1:8)
      end do

      return
      end

