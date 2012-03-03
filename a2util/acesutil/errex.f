
c This routine was the common exit handler for all AMEs, but that function has
c been transferred to aces_exit(). This routine is now kept for compatibility.

      subroutine errex
      call aces_exit(1)
      end

