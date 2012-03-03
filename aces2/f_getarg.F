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

c These routines localize the recasting of integer types for the system's iargc
c and getarg routines, which handle the command line arguments for Fortran
c programs.

c#define NO_GETARG
#if defined(_OSF) || defined(_CRAY_PVP)
#   define INTTYPE integer
#else
#   define INTTYPE integer*4
#endif

#ifdef NO_GETARG

      integer function f_iargc()
      f_iargc = 0
      return
      end

      subroutine f_getarg(ndx,sz)
      integer ndx
      character*(*) sz
      if (len(sz).gt.0) sz = " "
      return
      end

#else /* HAVE_GETARG */

      integer function f_iargc()
      INTTYPE iargc
#ifdef _UNICOS
      f_iargc = ipxfargc()
#else
      f_iargc = iargc()
#endif
      return
      end

      subroutine f_getarg(ndx,sz)
      integer ndx, ierr
      character*(*) sz
      INTTYPE i, i2, i3
      i = ndx
#ifdef _UNICOS
      call pxfgetarg(i,sz,i2,i3)
      if (i3.ne.0) then
         write(*,*) '@F_GETARG: pxfgetarg failed'
         write(*,*) '           i = ',i
         ierr = i3
         call c_exit(ierr)
      end if
#else
      call getarg(i,sz)
#endif
      return
      end

#endif /* NO_GETARG */

