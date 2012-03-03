



C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
C      This source file should not be edited.  Make
C      any necessary changes to the individual source
C      or include files.  This file has been produced
C      with 'make vmol2ja.f' from the original sources.
C      It is NOT the original source code itself.
C WARNING WARNING WARNING WARNING WARNING WARNING WARNING
      Subroutine NewPD(N, Cen, Ang, Map, Info)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: newpd.f,v 1.1.1.1 2003/04/02 19:21:47 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     newpd -- Create rearrangement necessary to group p, d, etc. fns.
C
C SYNOPSIS
      Integer N, Info
      Integer Cen(N), Ang(N), Map(N)
C ARGUMENTS
C DESCRIPTION
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer I, J, K, Cen1, Ang1, Blk1, LenBlk, NBlk
C
C
C     Statement fn for the number of cartesian functions a shel of ang.
C     momentum l.
C
      Integer NAng, l
      NAng(l) = (l+1)*(l+2)/2
C
C     First block must begin on first function
C
      Cen1 = Cen(1)
      Ang1 = Ang(1)
      Blk1 = 1
      LenBlk = 1
C
      Do 100 I = 2, N
         If ( Cen(i) .eq. Cen1 .AND. Ang(i) .eq. Ang1 ) then
C
C           **************************
C           * Same block -- continue *
C           **************************
C
            LenBlk = LenBlk + 1
         Else
C
C           *************************************
C           * New block.  First process old one *
C           *************************************
C           Best explained by example.  Say you have a block of
C           three p fns.  VMol, when using generally contracted
C           basis sets, will enumerate them as: x, x, x, y, y, y,
C           z, z, z.  We want them in the order x, y, z, x, y, z,
C           x, y, z.
C
C           To do this, we pick up every NBlk'th function from
C           the beginning of the block and place them in sequence.
C           We do this NAng(Ang1) times to collect a all the components
C           of a single function.  This is the inner loop.
C
C           The outer loop is merely to do this for every function in the
C           block.
C
C           By starting the inner loop with the value of the outer one,
C           we pick up, in turn, each of the NBlk functions with the
C           same (l,m).
C
C           This even works for s fns, but is a little inefficient
C
            K = NAng( Ang1 )
            NBlk = LenBlk/NAng(Ang1)
C
C           If this block does not have the structure we expect, 
C           leave it as is and flag an error.
C
            If ( Mod( LenBlk, NAng(Ang1) ) .ne. 0) then
               Info = Info + 1
               NBlk = LenBlk
               Ang1 = 0
            EndIf
C
            Do 200 j = Blk1, Blk1 + NBlk - 1
               Do 210 K = J, J + NBlk * (NAng(Ang1) - 1), NBlk
                  Map( Blk1 ) = K
                  Blk1 = Blk1 + 1
 210           Continue
 200        Continue
C
C           **************************
C           * Then setup the new one *
C           **************************
C
            Cen1 = Cen(i)
            Ang1 = Ang(i)
            Blk1 = i
            LenBlk = 1
         EndIf
 100  Continue
C
C     Do the last block
C
      NBlk = LenBlk/NAng(Ang1)
C     
C     If this block does not have the structure we expect, 
C     leave it as is and flag an error.
C     
      If ( Mod( LenBlk, NAng(Ang1) ) .ne. 0) then
         Info = Info + 1
         NBlk = LenBlk
         Ang1 = 0
      EndIf
C
      Do 300 j = Blk1, Blk1 + NBlk - 1
         Do 310 K = J, J + NBlk * (NAng(Ang1) - 1), NBlk
            Map( Blk1 ) = K
            Blk1 = Blk1 + 1
 310     Continue
 300  Continue
      Return
      End
