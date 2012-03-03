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
      subroutine decode_cond_code(cond_str, cond_code)
c---------------------------------------------------------------------------
c   Converts from the SIAL compiler's internal string used for conditionals 
c   to the SIP runtime code for the corresponding conditional operator.
c---------------------------------------------------------------------------

      implicit none
      integer i, ncond_table, cond_code
      parameter (ncond_table = 6)

      character*2 cond_str_table(ncond_table)
      character*2 my_cond_str
      character*(*) cond_str

      cond_str_table(1) = '=='
      cond_str_table(2) = '>='
      cond_str_table(3) = '<='
      cond_str_table(4) = '> '
      cond_str_table(5) = '< '
      cond_str_table(6) = '!='
      
      my_cond_str(1:1) = cond_str(1:1)
      if (cond_str(2:2) .eq. char(0)) then
         my_cond_str(2:2) = ' '
      else
         my_cond_str(2:2) = cond_str(2:2)
      endif 

      do i = 1, ncond_table
         if (my_cond_str .eq. cond_str_table(i)) then
            cond_code = i
            return
         endif
      enddo

      print *,'No decoded value'
      return
      end
      
      
      
