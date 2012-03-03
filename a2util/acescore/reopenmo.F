
cjp closes moints and moabcd and reopens them with new iref

      subroutine reopenmo
      implicit integer (a-z)
#include "bwcc.com"

      if (nref.gt.1) call loadmoio

      close(50)
      close(52)
      call aces_io_open(1)
      call aces_io_open(3)

      return
      end

