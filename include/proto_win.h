c---------------------------------------------------------------------------
c   Include common block for window table.
c---------------------------------------------------------------------------

      integer nwintab
      parameter (nwintab = 100)

      integer*8 wintab
      common /proto_win/wintab(6,nwintab)

