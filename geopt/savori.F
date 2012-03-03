C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      Subroutine SAVORI ( PtGrp, NAtms, Q, AtNr, At1, At2, XYOri, Err )
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: savori.F,v 1.3 2010/06/30 15:58:24 ponton Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:
C     Save information about the xy-orientation of the molecule.
C
C Arguments:
C     PtGrp    Point group to be used (Schoenflies symbol) (input)
C     NAtms    Number of atoms in molecule (input)
C     Q        Cartesian coordinates of centers (input)
C     AtNr     Atomic numbers of centers (input)
C     At1      First atom of "marker pair" (output)
C     At2      Second atom of "marker pair" (output)
C     XYOri    Direction of vector joining At1 to At2 (output)
C     Err      Error return code (output)
C
C Internal Variables:
C     ScrQ     Coordinates while they're being worked on. Length
C              2*(3*NAtms) because we need the workspace
C     RefOrd   Hold reference permutation list
C     NewOrd   Hold permutation list from sort
C     *** The GrpOp data structure ***
C     NrOps    Number of operations in the group.
C     PrinAx   Specify the principal axes for the molecule
C     OpTyp    Holds the type of symmetry operations from GrpOp
C     OpOrd    Holds the order of the axis for OpTyp = rotations
C     OpAxis   Holds the axis/plane to use for the sym. op.
C
C Dependents:
C     linblnk   Returns the index of the last non-blank character
C     GrpOps    Returns the operations of the given point group
C     MakGen   Updates the GenBy list
C     SortC    Sorts & mass weights coordinate vector.  We use it
C              to get the permutation list from the operations.
C     DoSyOp   The "symmetry operation interpreter" - does whatever
C              transformation its told.
C
C Diagnostic Messages:
C     "Unknown point group"  The given point group is unknown to
C     GrpOps. The C1 point group is used instead & is returned to
C     the caller in PtGrp.
C
C     "Unknown symmetry operation"  The sym. op. interpreter cannot
C     handle the indicated operation.  It is ignored.
C
C     "Unidentified error during..." The sym. op. interpreter returns
C     an error which SymEqv doesn't understand.  The operation is
C     ignored.  If this occurs, it probably means that the sym. op.
C     interpreter has been changed & SymEqv hasn't been updated to
C     match any new error codes.
C
C Limitations:
C     Parameter MxOps limits  the number of symmetry operations we
C     can handle from GrpOps.  It is presently set to 150, which should
C     be sufficient for anything but the I groups.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Maximum number of atoms currently allowed
#include "mxatms.par"
      Character*(*) PtGrp
      Dimension Q(3*NAtms)
      Integer AtNr(NAtms), At1, At2, XYOri, Err
C
      integer luout
      parameter (luout = 6 )
c
      Dimension OriVec(3)
C
C     To handle the coordinate transformations
      Dimension ScrQ(2*3*MxAtms)
      Integer RefOrd(MxAtms), NewOrd(MxAtms), GenBy(MxAtms)
C
C     For the operations of the point group.  MxOps=50 should handle
C     just about everything except I groups!
      Integer GrpOps
      Integer NrOps, MxOps, PrinAx(3)
      Parameter (MxOps = 150)
      Character*1 OpTyp(MxOps)
      Integer OpOrd(MxOps), OpAxis(MxOps)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     First, clear out the GenBy list, At1, At2, and XYOri.
C
      Do 100 i = 1, NAtms
         GenBy(i) = 0
 100  Continue
      At1 = 0
      At2 = 0
      XYOri = 0
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Now, get the operations of the group
C
      PrinAx(1) = 3
      PrinAx(2) = 1
      PrinAx(3) = 2
      NrOps = GrpOps (PtGrp, PrinAx, MxOps, OpTyp, OpOrd, OpAxis)
C
C     If GrpOps doesn't know this group, write an error message,
C     set the group to be C1, and return with everything marked
C     as unique in GenBy.
C
      If (NrOps .eq. -1) then
C         PtGrp = 'C1'
         Return
      EndIf
 9600 Format (' @SAVORI-W, Unknown point group ',A,'. Using C1')
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     If there are no non-trivial operations in the group, the
C     "right thing" is for GenBy to have everything as unique,
C     which it is, so we can return, printing the "final" GenBy
C     if required.
C
      If (NrOps .eq. 0) Return
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Generate the reference permutation list to compare with later.
C     Put resulting coordinates in ScrQ because we don't need them -
C     only the permutation list.
C
      Call SortC (NAtms, Q, AtNr, ScrQ, RefOrd, ScrQ(3*NAtms+1))
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     Look for a xz- or yz-plane which maps one center into another.
C
      Do 1000 Iop = 1, NrOps
         If (OpTyp(Iop) .eq. 'P' .AND. OpAxis(Iop) .ne. 3) then
            Call DoSyOp ( OpTyp(Iop), OpOrd(Iop), J, OpAxis(Iop),
     $         NAtms, Q, ScrQ, Err, 0)
C
C           Handle any errors in DoSyOp
C
            If (Err .eq. 1) then
               Write (Luout, 9610) OpTyp(Iop), OpOrd(Iop), OpAxis(Iop)
               Goto 1000
            ElseIf (Err .ne. 0) then
               Write (Luout, 9620) OpTyp(Iop), OpOrd(Iop), OpAxis(Iop)
               Goto 1000
            EndIf
 9610       Format (' @SAVORI-W, Unknown symmetry operation ',
     $         A,'(',I3,') about ',I3,': Operation ignored.')
 9620       Format (' @SAVORI-W, Unidentified error during, ',
     $         'operation ',A,'(',I3,') about ',I3,': Op. ignored.')
C
C           Generate the new permutation list & update GenBy.
C
            Call SortC (NAtms, ScrQ, AtNr, ScrQ, NewOrd,
     $         ScrQ(3*NAtms+1))
C
            Call MakGen ( NAtms, RefOrd, NewOrd, GenBy)
C
C     See if this operation maps one center onto another. If so, get
C     the centers, and the normalized vector between them.  This vector
C     should be in either the x or y direction.
C
            Do 1100 i = 1, NAtms
               If (GenBy(i) .ne. 0) then
                  At1 = i
                  At2 = GenBy(i)
                  Call VAdd(OriVec, Q(3*(At2-1)+1), Q(3*(At1-1)+1),3,
     $                            -1.D0)
                  Call Normal(OriVec,3)
                  XYOri = NInt(OriVec(1)) + 2*NInt(OriVec(2))
     $               + 3*NInt(OriVec(3))
cYAU                  Write (6,*) 'OriVec is ',OriVec, XYOri
                  Return
               Endif
 1100       Continue
         EndIf
 1000 Continue
C
C     If we leave this way, there are no equivalences we can use.
C     This is possible only for cases in which we don't need to worry
C     about x <--> y switching anyway, so its okay.
C
      Return
      End
