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
      Integer Function GrpOps ( PtGrp, PrimAx, Length, Type, Order,
     $   Axis )
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:
C     Given a point group, generate a list of operations of the group.
C
C Limitations:
C     Knows the abelian groups: C1, Cs, Ci, C2, D2, C2v, C2h, D2h
C     Knows all Cn, Cnh, Sn groups
C     Doesn't know any group with C2's or sigma's *not* along the
C     cartesian axes due to problems defining such operations with
C     the present data structure.
C
C Returns:
C     GrpOps   Returns the number of operations in the group.
C              In case of error, returns error codes a follows:
C              = -1  Group unknown to GrpOps (presently only abelian
C                    groups and selected others are known.)
C              NOTE: The operation count given here is *not* the same
C                    as the usual group theory usage.  For example,
C                    a C6 axis is given as a single operation by GrpOps
C                    and the code *using* this list is expected to
C                    iterate through all 6 of the C6's.
C
C Variables:
C     PtGrp    String containing the Schoenflies symbol for the
C              group of interest. (input)
C              The first character must be in [CDS].  The order
C              must begin with the second character of PtGrp.  The
C              last non-blank character is taken to be the "modifier"
C              in the group symbol. It may be in upper- or lowercase.
C              A group of infinite order is denoted with an 'X' for the
C              order.
C     PrimAx   Array indicating the primary, secondary, (input)
C              and tertiary axes.
C     Length   Length of Type/Order/Axis arrays (input)
C     Type     Type of symmetry operation. (output)
C              = 'C'   Rotation
C              = 'S'   Improper rotation
C              = 'P'   Reflection plane
C              = 'I'   Inversion center
C     Order    Order of rotation axis. (output)
C     Axis     Axis of rotation or plane of reflection. (output)
C              = 1     x-axis or yz-plane
C              = 2     y-axis or xz-plane
C              = 3     z-axis or xy-plane
C
C Bugs:
C     Doesn't check point group symbols carefully.
C     Requires that perpendicular C2s and reflection planes
C     coincide with the cartesian axes.
C
C Dependents:
C     AtoI     Converts string to integer
C     linblnk   Returns the index of the last non-blank character
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C $Log: grpops.f,v $
C Revision 1.4  2010/10/07 14:33:58  ponton
C Add symcor changes, fix bugs in optimization code, change numerical lib calls to use proper routines
C
C Revision 1.2  2010/02/10 17:20:48  ponton
C Add GNU GPL info to each source file
C
C Revision 1.1.1.1  2009/07/01 18:54:34  ponton
C Initial import for ACESIII Release 3.0
C
C Revision 1.1.1.1  2003/04/02 19:21:35  aces
C INITIAL 2.4 IMPORT
C
C Revision 4.0  89/03/14  01:15:38  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C 
C Revision 3.0  89/01/29  23:10:12  bernhold
C First working release for VAX
C 
C Revision 1.1  89/01/28  12:08:51  bernhold
C Initial revision
C 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Implicit Integer (a-z)
      Character*(*) PtGrp
      Integer PrimAx(3)
      Character*1 Type(length)
      Dimension Order(length), Axis(length)
C
      Character*1 PGMod
      Integer AtoI, linblnk
C
C     GrpOps serves as a counter for Type, Order, and Axis
C     PGMod is the "modifier" lifted from the name
C     OrdGrp is the order of the group (-1 if infinite order)
C
      GrpOps = 0
C
      i = linblnk(PtGrp)
      PGMod = PtGrp(i:i)
C
      If (PtGrp(2:2) .eq. 'X') then
         OrdGrp = -1
      Else
         OrdGrp = AtoI( PtGrp(2:) )
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     First check for the special cases: C1, Cs, Ci
C
      If (PtGrp(1:1) .eq. 'C' .AND. OrdGrp .eq. 1) then
C        C1 has *no* operations other than E
         GrpOps = 0
         Return
      ElseIf (PtGrp(1:2) .eq. 'Cs' .OR. PtGrp(1:2) .eq. 'CS'
     &   .OR. PtGrp .eq. 'C s') then
C        Cs has only Sigma-h (xy-plane)
         GrpOps        =  1
         Type(GrpOps)  = 'P'
         Order(GrpOps) =  0
         Axis(GrpOps)  =  PrimAx(1)
         Return
      ElseIf (PtGrp(1:2) .eq. 'Ci' .OR. PtGrp(1:2) .eq. 'CI'
     &      .OR. PtGrp .eq. 'C i') then
C        Ci has only i (inversion center)
         GrpOps        =  1
         Type(GrpOps)  = 'I'
         Order(GrpOps) =  0
         Axis(GrpOps)  =  0
         Return
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Handle the linear groups: Don't worry about C-infinity!
C     If we know its linear already, then only the non-axial
C     operations will help!
C
      If (PtGrp(1:3) .eq. 'CXv' .OR. PtGrp(1:3) .eq. 'CXV') then
C        CXv has no non-axial operations
         GrpOps = 0
         Return
      EndIf
      If (PtGrp(1:3) .eq. 'DXh' .OR. PtGrp(1:3) .eq. 'DXH') then
C        DXh has a sigma-h plane
         GrpOps        =  1
         Type(GrpOps)  = 'P'
         Order(GrpOps) =  0
         Axis(GrpOps)  =  PrimAx(1)
         Return
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     If we got this far, and its order is not at least 2, then
C     it must be T, O, I (which we can't handle) or a screw-up
C     with the point group symbol, so get out!
C
      If (OrdGrp .le. 1) then
         GrpOps = -1
         Return
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Sn groups are special: n must be even
C
      If (PtGrp(1:1) .eq. 'S') then
         If (mod(OrdGrp, 2) .eq. 0) then
            GrpOps        =  1
            Type(GrpOps)  = 'S'
            Order(GrpOps) =  OrdGrp
            Axis(GrpOps)  =  PrimAx(1)
            Return
         Else
            GrpOps = -1
            Return
         EndIf
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     All other groups (C & D) have a C(OrdGrp) about the z-axis
C
      If (PtGrp(1:1) .eq. 'C' .OR. PtGrp(1:1) .eq. 'D') then
         GrpOps        =  GrpOps + 1
         Type(GrpOps)  = 'C'
         Order(GrpOps) =  OrdGrp
         Axis(GrpOps)  =  PrimAx(1)
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     If its a 'D' group, add the perpendicular C2's
C
      If (PtGrp(1:1) .eq. 'D') then
         If (OrdGrp .eq. 2) then
C           Can handle D2, D2h
            Do 100 i = 2, 3
               GrpOps        = GrpOps + 1
               Type(GrpOps)  = 'C'
               Order(GrpOps) =  OrdGrp
               Axis(GrpOps)  =  PrimAx(i)
 100        Continue
         Else
C           Can't handle orders > 2 - they have C2's which
C           aren't on the cartesian axes.
            GrpOps = -1
            Return
         EndIf
      EndIf
C
C IF IT IS D2H, ADD ALL OF THE REMAINING OPS.
C
      IF(PTGRP.EQ.'D2h'.OR.PTGRP.EQ.'D2H')THEN
       DO 133 I=2,3
        GRPOPS=GRPOPS+1
        TYPE(GRPOPS)='P'
        ORDER(GRPOPS)=0
        AXIS(GRPOPS)=PRIMAX(I)
133    CONTINUE
       GRPOPS=GRPOPS+1
       TYPE(GRPOPS)='I'
       ORDER(GRPOPS)=0
       AXIS(GRPOPS)=0
      ENDIF
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Add the horizontal plane for 'h' groups
C
      If (PGMod .eq. 'h' .OR. PGMod .eq. 'H') then
         GrpOps        =  GrpOps + 1
         Type(GrpOps)  = 'P'
         Order(GrpOps) =  0
         Axis(GrpOps)  =  PrimAx(1)
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Add the vertical planes for 'v' groups
C
      If (PGMod .eq. 'v' .OR. PGMod .eq. 'V') then
         If (OrdGrp .eq. 2) then
C           Can handle C2v
            Do 200 i = 2, 3
               GrpOps        = GrpOps + 1
               Type(GrpOps)  = 'P'
               Order(GrpOps) =  0
               Axis(GrpOps)  =  PrimAx(i)
 200        Continue
         Else
C           Can't handle orders > 2 - they have C2's which
C           aren't on the cartesian axes.
            GrpOps = -1
            Return
         EndIf
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     The only other thing left is a 'd' group, which we can't
C     handle either
C
      If (PGMod .eq. 'd' .OR. PGMod .eq. 'D') then
         GrpOps = -1
         Return
      EndIf
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     If we get here and haven't got any operations, something's wrong!
C
      If (GrpOps .eq. 0) GrpOps = -1
C
      Return
      End
