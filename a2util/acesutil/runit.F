
      subroutine runit(string)
      character*(*) string
      integer length, linblnk
      external linblnk

c#ifdef _DEBUG
      length=linblnk(string)
      write(*,'(3a)') '@ACES2: Executing "',string(1:length),'"'
c#endif

      istat=ishell(string)
cYAU - Somehow, exit(istat) doesn't tell the OS the program failed.
c      So, running `xaces2 || echo failure` never worked.
c      if (istat.ne.0) call exit(istat)
      if (istat.ne.0) call c_exit(1)

      return
      end

