
cjp for management of multiple moio information for multiple references

      subroutine loadmoio
      IMPLICIT INTEGER (A-Z)

#include "lists.com"
#include "bwcc.com"

      iLenMOIO = _MAX_IO_GRPS * _MAX_IO_FAMS
      call icopy(iLenMOIO, bwmoio  (1,1,iref), 1, moio  (1,1), 1)
      call icopy(iLenMOIO, bwmoiowd(1,1,iref), 1, moiowd(1,1), 1)
      call icopy(iLenMOIO, bwmoiosz(1,1,iref), 1, moiosz(1,1), 1)
      call icopy(iLenMOIO, bwmoiods(1,1,iref), 1, moiods(1,1), 1)
      call icopy(iLenMOIO, bwmoiofl(1,1,iref), 1, moiofl(1,1), 1)
      call icopy(iLenMOIO, bwmoiomxsz(1,1,iref), 1, moiomxsz(1,1), 1)
      do i = 1, _MAX_IO_LUNS
         pRec(i) = bwpRec(i,iref)
         iIntOff(i) = bwiIntOff(i,iref)
      end do

      return
      end

