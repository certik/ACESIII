
c This routine prints the unique factors of a number.

      subroutine factor(args,dimargs)
      implicit none

c ARGUMENT LIST
      integer dimargs
      character*80 args(dimargs)

c INTERNAL VARIABLES
      double precision d
      integer iArg
      integer iMult, iFact1, iFact2, iList(100), nFact, i

c EXTERNAL FUNCTIONS
      integer atoi
      external atoi

c ----------------------------------------------------------------------

      if (dimargs.eq.0) return

      do iArg = 1, dimargs

         iMult = atoi(args(iArg))
         i = sqrt(1.*(iMult+0.1))
         nFact = 0
         do iFact1 = i, 1, -1
            iFact2 = iMult/iFact1
            if (iFact1*iFact2.eq.iMult) then
               iList(nFact+1) = iFact1
               iList(nFact+2) = iFact2
               nFact = nFact + 2
               if (iFact1.eq.iFact2) nFact = nFact - 1
            end if
         end do
         call isort(iList,0,nFact,1)
         print *, iMult,' -> ',(iList(i),i=1,nFact)

      end do

      return
      end

