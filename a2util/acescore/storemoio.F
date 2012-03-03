
cjp for management of multiple moio information for multiple references

      subroutine storemoio
      IMPLICIT INTEGER (A-Z)

#include "lists.com"
#include "bwcc.com"

      iLenMOIO = _MAX_IO_GRPS * _MAX_IO_FAMS
      call icopy(iLenMOIO, moio  (1,1), 1, bwmoio  (1,1,iref), 1)
      call icopy(iLenMOIO, moiowd(1,1), 1, bwmoiowd(1,1,iref), 1)
      call icopy(iLenMOIO, moiosz(1,1), 1, bwmoiosz(1,1,iref), 1)
      call icopy(iLenMOIO, moiods(1,1), 1, bwmoiods(1,1,iref), 1)
      call icopy(iLenMOIO, moiofl(1,1), 1, bwmoiofl(1,1,iref), 1)
      call icopy(iLenMOIO, moiomxsz(1,1), 1, bwmoiomxsz(1,1,iref), 1)
      do i = 1, _MAX_IO_LUNS
         bwpRec(i,iref) = pRec(i) 
         bwiIntOff(i,iref) = iIntOff(i) 
      end do

      return
      end

