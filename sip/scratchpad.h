c-----------------------------------------------------------------------
c   Declarations for the data "scratchpad".
c-----------------------------------------------------------------------
      integer mx_scratchpad
      integer mx_fscratchpad
      parameter (mx_scratchpad = 2000)
      parameter (mx_fscratchpad = 2000)
      integer scratchpad
      double precision fscratchpad
      common /scratchpad/scratchpad(mx_scratchpad),
     *                   fscratchpad(mx_fscratchpad)

      
