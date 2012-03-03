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

c This routine returns the DACOS of VALUE. If VALUE is within THRESH of 1.d0 or
c -1.d0, then VALUE is set to (-)1.d0 and the DACOS is evaluated and returned.

      double precision function dacosx(dValue,dThresh)
      implicit none
      double precision dValue, dThresh
      if (dabs(dabs(dValue)-1.d0).lt.dThresh) dValue=sign(1.d0,dValue)
      dacosx=dacos(dValue)
      return
      end

