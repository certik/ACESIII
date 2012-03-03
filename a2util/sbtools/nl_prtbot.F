
c This routine prints the bottom of a table of namelist values.

      subroutine nl_prtbot
      implicit none
      write(*,900)
      write(*,990)
  900 format(76('='))
  990 format(' ')
      return
      end

