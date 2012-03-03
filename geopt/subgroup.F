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
      SUBROUTINE SUBGROUP(INGROUP,ABGRP0,FUGRP0,IROTATE)
C
C THIS ROUTINE ORIENTS THE MOLECULE IN TERMS OF AN EXTERNALLY
C SPECIFIED ABELIAN SUBGROUP AND PROCEEDS TO CHECK ITS OPERATIONS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*3 GRPS(8),GROUP
      CHARACTER*4 ABGRP,ABGRP0,FUGRP,FUGRP0
C
      DATA GRPS /'C1 ','C2 ','C s','C i','C2v','C2h','D2 ','D2h'/
C
CSSS      LASTF=linblnk(FUGRP)
C A bug fix. This subroutine knows nothing about FUGRP. I am suprised
C that did not show up eariler than this, but I guess it is rare that
C we do frequencies with subgroup is specified. Ajith Perera, 05/2006.
C
      LASTF=linblnk(FUGRP0)
      GROUP=GRPS(INGROUP)
      ABGRP=ABGRP0(1:3)
      FUGRP=FUGRP0(1:LASTF)
c
#ifdef _DEBUG_LVL0
      Print*, "Abelian group, full group"
      Print*, ABGRP, FUGRP, GROUP
#endif

C
      IROTATE=0
C
      IF(ABGRP.EQ.GROUP)THEN
C
C QUIETLY RETURN IF SPECIFIED GROUP IS IDENTICAL TO ABELIAN SUBGROUP.
C
       ABGRP0(1:3)=GROUP
       RETURN
C
      ELSEIF(ABGRP.EQ.'D2h')THEN
C
C CASE I.  TRUE ABELIAN SUBGROUP IS D2H.  THEN ALL SUBGROUPS ARE
C OK AND ALL ARE PROPERLY ORIENTED.
C
       ABGRP0(1:3)=GROUP
       RETURN
C
      ELSEIF(ABGRP.EQ.'C2v')THEN
C
C CASE II.  TRUE ABELIAN SUBGROUP IS C2V.  ALLOWABLE GROUPS ARE
C C2 AND CS, BOTH OF WHICH ARE PROPERLY ORIENTED.
C
       IF(GROUP.EQ.'C2 '.OR.GROUP.EQ.'C s')THEN
        ABGRP0(1:3)=GROUP
        RETURN
       ELSEIF(GROUP.EQ.'D2 '.AND.FUGRP.EQ.'T d'.OR.
     &        GROUP.EQ.'D2 '.AND.FUGRP(LASTF:LASTF).EQ.'d')THEN
C
C TRICKY STUFF.  HAVE TO ROTATE BY 45 DEGREES ABOUT Z.
C
        IROTATE=1
        ABGRP0(1:3)=GROUP
       ELSE
        WRITE(6,1000)GROUP,FUGRP
        CALL ERREX
       ENDIF
C
      ELSEIF(ABGRP.EQ.'D2 ')THEN
C
C CASE III.  TRUE ABELIAN SUBGROUP IS D2.  ALLOWABLE GROUP IS C2.
C
       IF(GROUP.EQ.'C2 ')THEN
        ABGRP0(1:3)=GROUP
        RETURN
       ELSEIF(GROUP.EQ.'C2v'.AND.FUGRP.EQ.'T d'.OR.
     &        GROUP.EQ.'C2v'.AND.FUGRP(LASTF:LASTF).EQ.'d')THEN
        IROTATE=1
        ABGRP0(1:3)=GROUP
       ELSE
        WRITE(6,1000)GROUP,FUGRP
        CALL ERREX
       ENDIF
      ELSEIF(ABGRP.EQ.'C2h')THEN
C
C CASE IV.  TRUE ABELIAN SUBGROUP IS C2H.  ALLOWABLE GROUPS ARE C2,
C CS AND CI.
C
       IF(GROUP.EQ.'C2 '.OR.GROUP.EQ.'Ci '.OR.GROUP.EQ.'C s')THEN
        ABGRP0(1:3)=GROUP
        RETURN
       ELSE
        WRITE(6,1000)GROUP,FUGRP
        CALL ERREX
       ENDIF
C
C FOR ALL OTHER ABELIAN SUBGROUPS, THE INPUT GROUP MUST MATCH THE
C TRUE ABELIAN SUBGROUP, WHICH IS TRAPPED AT THE TOP OF THIS ROUTINE.
C
      ELSE
       WRITE(6,1000)GROUP,FUGRP
C       CALL ERREX
C
C There is really no reason to stop; calculation can proceed in the
C abelian sub-group, but we need to warn the user that we disobyed
C him/her when one or few finite difference points happen to have
C point group for which the user requested subgroup is no longer a
C a subgroup of the full group. 05/2006, Ajith Perera

      CALL WRITE_WARN_LOG("INVALID SUBGROUP")

      ENDIF
C
1000  FORMAT(T3,'@SUBGROUP-F, Group ',A,' is not a subgroup of ',A,'.')
      RETURN
      END
