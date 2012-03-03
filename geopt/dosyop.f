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
      Subroutine DoSyOp (Type, Order, Count, Axis, NAtms, Q, NewQ, Err,
     &                   IMode)
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose: Perform the given symmetry transformation on Q
C
C Usage:
C     To do all of the symmetry operations of a given point group:
C     Call GrpOps (...)
C     Do Op = 1, NrOps
C        Do Count = 1, Max(1, Order(Op))
C           Call DoSyOp (Type(Op), Order(Op), Count,...)
C           Actions using NewQ
C        EndDo
C     EndDo
C
C Limitations:
C     Knows all proper & improper rotations, reflections using the
C     cartesian axes.  Knows inversion centers. CANNOT do rotations
C     and reflections which do not coincide with the cartesian axes.
C
C Arguments:
C     Type    Type of symmetry operation (input)
C     Order   Order of the axis (0 if not an axis) (input)
C     Count   Do the Count-th rotation defined by the op. (input)
C             Ignored if Type is not a rotation.
C     Axis    Cartesian axis/plane for rotation/reflection (input)
C     NAtms   Number of atoms - length of Q (input)
C     Q       Cartesian coordinates to be transformed (input)
C     NewQ    Resulting cartesian coordinates (output)
C     Err     Non-zero in case of any error (output)
C             = 1   Unknown operation specified.
C     IMode   = 0 for real Sn.
C              = 1 means that S_n^q is actually C_n followed by
C                a reflection.  Although this is really not what
C                a S_n^q is, it makes generation of symops easier.
C
C Internal Variables:
C     RMat    Rotational transformation matrix
C
C Dependents:
C     RotM    Make a rotation matrix for any rotation
C     MatMulV Matrix-vector multiply
C     Reflect Reflect coordinates through a plane
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C $Log: dosyop.f,v $
C Revision 1.4  2010/10/07 14:33:58  ponton
C Add symcor changes, fix bugs in optimization code, change numerical lib calls to use proper routines
C
C Revision 1.2  2010/02/10 17:20:47  ponton
C Add GNU GPL info to each source file
C
C Revision 1.1.1.1  2009/07/01 18:54:34  ponton
C Initial import for ACESIII Release 3.0
C
C Revision 1.1.1.1  2003/04/02 19:21:35  aces
C INITIAL 2.4 IMPORT
C
C Revision 4.3  89/09/21  14:59:48  jstanton
C   Fixed bug in Sn operations where reflection was performed
C regardless of the operation "count".
C 
C Revision 4.2  89/08/19  13:42:13  jstanton
C Add missing declaration of array NewQ
C 
C Revision 4.1  89/05/15  11:01:05  jstanton
C Replaced ICCOD (bug!) with ICOORD.
C 
C Revision 4.0  89/03/14  01:15:28  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C 
C Revision 3.0  89/01/29  23:09:43  bernhold
C First working release for VAX
C 
C Revision 1.1  89/01/28  12:07:35  bernhold
C Initial revision
C 
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Implicit double precision (A-H, O-Z)
      Character*1 Type
      Double Precision NewQ
      Double Precision ONE,ZILCH
      Integer Order, Count, Axis, Natms, Err
      Dimension Q(3*NAtms), NewQ(3*NAtms)
C
      Dimension RMat(3,3)
C
      DATA ONE,ZILCH /1.0D+00,0.0D+00/
C
C     Initialize things
C
      Err = 0
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     
C     Proper and improper rotations (C and S)
C     
      If (Type .eq. 'C' .OR. Type .eq. 'S') then
         Angle = Count*360.0/Order
         Call RotM (Axis, Angle, 1, RMat)
C
CJDW 1/6/98. MATMULV call replaced by XGEMM call. NewQ=(RMat)^T * Q.
C
C        Call MatMulV (NewQ, Q, RMat, NAtms, 3, 3)
C
         CALL XGEMM('T','N',3,NAtms,3,ONE,RMat,3,Q,3,ZILCH,NewQ,3)
C
         If (Type.eq.'S'.and.Mod(Count,2).eq.1.Or.Type.eq.'S'.
     &         and.IMode.eq.1)then
            Call Reflect (NewQ, NewQ, NAtms, Axis)
         EndIf
C     
C     Reflection planes
C     
      ElseIf (Type .eq. 'P') then
         Call Reflect (Q, NewQ, NAtms, Axis)
C     
C     Inversion
C     
      ElseIf (Type .eq. 'I') then
         Do 1200 icoord = 1, 3*NAtms
            NewQ(icoord) = - Q(icoord)
 1200    Continue
C     
C     Unknown operation - set Err & return
C     
      Else
         Err = 1
      EndIf
      Return
      End
