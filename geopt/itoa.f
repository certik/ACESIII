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
       CHARACTER*(*) FUNCTION ITOA (NR, FRCPLS)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Convert NR to a left justified string
C
C Arguments:
C     NR       number to be converted (input only)
C     FRCPLS   Force leading '+' if NR positive (input only)
C
C Limitations:
C     May return with incomplete conversion if length of ITOA is too
C     short.  Puts '*' in last position of ITOA to indicate overlow.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C $Log: itoa.f,v $
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
C Revision 4.0  89/03/14  01:15:45  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C 
C Revision 3.0  89/01/29  23:10:22  bernhold
C First working release for VAX
C 
C Revision 2.1  89/01/02  20:36:12  bernhold
C To keep consistent with .u file just checked in.
C 
C     Revision 1.1  88/12/07  13:38:51  bernhold
C     Initial revision
C     
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      INTEGER NR, FRCPLS, NRABS
      INTEGER I, J
C
C     Clear out the string
C
      DO 10 I = 1, LEN(ITOA)
         ITOA = ' '
 10   CONTINUE
C
C     Start counting position in string
C
      J = 1
      NRABS = ABS (NR)
C
C     Put in sign as appropriate
C
      IF (NR .LT. 0) THEN
         ITOA(J:J) = '-'
         J = J + 1
      ENDIF
      IF (FRCPLS .NE. 0 .AND. NR .GT. 0) THEN
         ITOA(J:J) = '+'
         J = J + 1
      ENDIF
C
C     Check if we are about to overflow the string
C
      IF (J .GT. LEN(ITOA)) THEN
         ITOA(J-1:J-1) = '*'
         RETURN
      ENDIF
C
C     Loop over nr of digits in number
C
      NDIG = INT( LOG10( FLOAT(NRABS)) ) + 1
      DO 100 I = NDIG, 1, -1
         N = MOD ( ( NRABS / (10**(I-1) ) ), 10)
         ITOA(J:J) = CHAR(N + 48)
         J = J + 1
C
C        Check for overflow of the string, but if this is last digit
C        then its okay.
C
         IF (J .GT. LEN(ITOA) .AND. I .GT. 1) THEN
            ITOA(J-1:J-1) = '*'
            RETURN
         ENDIF
 100  CONTINUE
C
      RETURN
      END
