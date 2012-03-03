
c This routine prints the top of a table of namelist values.

      subroutine nl_prttop(nltitle)
      implicit none
      character*(*) nltitle
      write(*,990)
      write(*,900)
      write(*,910) nltitle
      write(*,920) 'KEYWORD','TYPE','VALUE','DEFAULT'
      write(*,930)
      write(*,920) 'print_nl','logical','true'
  990 format(' ')
  900 format(76('='))
  910 format('NAMELIST: ',a)
  920 format(a20,2x,a10,2x,a20,2x,a20)
  930 format(20('-'),2x,10('-'),2x,20('-'),2x,20('-'),2x)
      return
      end

