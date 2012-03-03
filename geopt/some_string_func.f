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
      CHARACTER*4 FUNCTION some_string_func(NAME,Iord)
      CHARACTER*4 NAME, ItoA
      If (Name(2:2) .eq. 'N') then
         some_string_func(1:1) = Name(1:1)
         some_string_func(2: ) = ItoA(Iord, 0)
         some_string_func(linblnk(some_string_func)+1:) = Name(3:3)
      Else
         some_string_func = Name
      EndIf
      RETURN
      END
