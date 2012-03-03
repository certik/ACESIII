      Subroutine a2_reset_jarc()

      Implicit Double Precision (A-H, O-Z)
C

      Call Putrec(1,'JOBARC','HAVEGEOM', 1, 1)
      Call Getrec(0,'JOBARC','PASS1   ', Length, 1)
C
      Call aces_ja_truncate('JODAOUT ', 1)
      If (Length .gt. 0) Call Putrec(0,'JOBARC','PASS1   ', 1, 0)
C
      Call Putrec(1,'JOBARC','PES_SCAN', 1, 1)
C
      Call aces_ja_fin
      Call aces_ja_init
c
      Return 
      End
