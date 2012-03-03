
c This routine shuts down the ACES I/O subsystem in the SB environment.

      subroutine sb_io_fin
      implicit none
      call callstack_push('SB_IO_FIN')
cTODO      call aces_io_fin
      call callstack_pop
      return
      end

