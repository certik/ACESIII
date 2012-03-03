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
      function is_connected(atom1, atom2, connectiontable, natoms,
     &                      connectlist)
C
C Determines if two atoms are connected by constructing a full
C connection list and then searching it.
C
      implicit none
C
C Input Variables
C
      integer atom1,atom2,natoms,connectiontable(natoms,natoms)
C
C Function Variables
C
      logical is_connected
C
C Pre-allocated Local Variables
C
      integer connectlist(natoms), i
C
C External Procedures
C
      external followconnect

C Local Variables

      integer listatoms
      logical shortcircuit
C - - -- - - - - - - - - - - - - -- - - - - - - -- - - -  -- - 

      if (connectiontable(atom1,atom2).eq.1) then
         is_connected=.true.
         return
      else
         call izero(connectlist,natoms)
C
C Follow atom1 connections
C
         listatoms=0
         shortcircuit=.false.
         call followconnect(atom2,atom1,connectlist,listatoms,
     &      connectiontable,natoms,
     &      shortcircuit)
         is_connected=shortcircuit
      endif
C
#ifdef _DEBUG_LVLM1
      Print*, "The connectionlist in is_connected"
      Print*, "atom1 & atom2", (connectlist(i), i=1, natoms),
     &         atom1, atom2
#endif
c
      return
      end
C
      subroutine followconnect(atom,taratom,list,listatoms,
     &   connectiontable,natoms,
     &   shortcircuit)
C
C Constructs list of connections non-recursively
C
      implicit none
C Input Variables
      integer atom,natoms,connectiontable(natoms,natoms),taratom

C Input/Output Variables
      integer list(natoms),listatoms
      logical shortcircuit

C Local Variables
      integer ii, i, j, k
      integer paths(natoms,natoms)
      integer path_count(natoms)
      integer extension(natoms,natoms)
      integer extension_count(natoms)
      integer ipath, npaths, nextend
      

C - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - - -- - -
      
      npaths = 0

c----------------------------------------------------------------------------
c   Find all paths of length 1.
c----------------------------------------------------------------------------

      do ii = 1, natoms
         if (connectiontable(ii,atom) .eq. 1) then

c----------------------------------------------------------------------------
c   Add this path.
c----------------------------------------------------------------------------

            npaths = npaths + 1
            paths(1,npaths) = ii    ! path is from iatom to ii.
            path_count(npaths) = 1  ! length of path is 1.
         endif
      enddo

      call remove_invalid_paths(paths, path_count, natoms, npaths, 
     *                          atom)
      call check_paths(paths, path_count, natoms, npaths, atom, 
     *                           taratom, list, listatoms, 
     *                           shortcircuit)
      if (npaths .eq. 0 .or. 
     *    shortcircuit) return

  100 continue
      nextend = 0
      do ipath = 1, npaths

c----------------------------------------------------------------------------
c   Extend paths if possible.
c----------------------------------------------------------------------------

         ii = paths(path_count(ipath),ipath)   ! last atom in path.
         do i = 1, natoms
            if (connectiontable(i,ii) .eq. 1) then
               nextend = nextend + 1
               if (nextend .gt. natoms) then
                  print *,'Error: followconnect: Attempting to create ',
     *               'more path extensions than natoms ',natoms
                  print *,'atom, taratom ',atom,taratom
                  print *,'connectiontable: ii ',ii
                  do j = 1, natoms
                     print *,'ATOM ',j
                     print 10101,(connectiontable(j,k),k=1,natoms)
10101 format(10(1x,i3))
                  enddo      
                  call abort_job()
               endif

               extension_count(nextend) = path_count(ipath)+1
               do j = 1, path_count(ipath)
                  extension(j,nextend) = paths(j,ipath)
               enddo
               extension(extension_count(nextend),nextend) = i
            endif
         enddo   ! i

c---------------------------------------------------------------------------
c   Cull out the invalid extensions.
c---------------------------------------------------------------------------

         call remove_invalid_paths(extension, extension_count, 
     *                          natoms, nextend,
     *                          atom) 
      enddo ! ipath

c---------------------------------------------------------------------------
c   Copy the extended paths into the paths array.
c---------------------------------------------------------------------------

      npaths = nextend
      do j = 1, npaths
         path_count(j) = extension_count(j)
         do i = 1, path_count(j)
            paths(i,j) = extension(i,j)
         enddo 
      enddo
      
      call remove_invalid_paths(paths, path_count, natoms, npaths,
     *                          atom)
      call check_paths(paths, path_count, natoms, npaths, atom,
     *                           taratom, list, listatoms,
     *                           shortcircuit)
      if (npaths .eq. 0 .or.
     *    shortcircuit) return
      go to 100

      return
      end

      subroutine remove_invalid_paths(paths, path_count, natoms, 
     *                          npaths, iatom)
c-----------------------------------------------------------------------------
c   Checks for degenerate paths in the "paths" array and removes them.
c   A degenerate path is one that has looped back to the iatom, or one 
c   which contains duplicate atoms.
c-----------------------------------------------------------------------------

      implicit none
      integer natoms, iatom, npaths
      integer paths(natoms,natoms) 
      integer path_count(natoms)
      integer i, j, k, nremove, remove(natoms)
      integer kpath, lenpath
      logical remove_this

      if (npaths .eq. 0) return   ! nothing to do.

      nremove = 0
      do j = 1, npaths
         lenpath = path_count(j)

         remove_this = .false.
         do i = 1, lenpath

c----------------------------------------------------------------------------
c   Check path for loopback to iatom.
c----------------------------------------------------------------------------

            if (paths(i,j) .eq. iatom) then
               remove_this = .true.
               nremove = nremove + 1
               remove(nremove) = j
            endif
         enddo 

         if (.not. remove_this) then

c----------------------------------------------------------------------------
c   Check for duplicate atoms in path.
c----------------------------------------------------------------------------

            do i = 1, lenpath
               do k = i+1, lenpath
                  if (paths(i,j) .eq. paths(k,j)) then
                     nremove = nremove + 1
                     remove(nremove) = j
                  endif
               enddo     ! k
            enddo        ! i
         endif           ! .not. remove_this
      enddo              ! j

c----------------------------------------------------------------------------
c   Remove all paths marked for removal.
c----------------------------------------------------------------------------

      if (nremove .gt. 0) then
         do j = nremove, 1, -1  ! must be removed in reverse order.
            k = remove(j)     ! index of path to be removed.

c---------------------------------------------------------------------------
c   Move paths k+1 through npaths to k, npaths-1.
c---------------------------------------------------------------------------

            do kpath = k+1, npaths
               lenpath = path_count(kpath)
               path_count(kpath-1) = lenpath

               do i = 1, lenpath
                  paths(i,kpath-1) = paths(i,kpath)
               enddo
            enddo

            npaths = npaths - 1
         enddo
      endif

      return
      end

      subroutine check_paths(paths, path_count, natoms, npaths, iatom,
     *                           taratom, list, listatoms,
     *                           found_connection)
c---------------------------------------------------------------------------
c   Checks all paths to determine if a connection between iatom and taratom
c   has been found.  Returns found_connection = .true. if found, .false.
c   otherwise.  If found_connection = .true., the actual path is returned
c   in list and the length of the path is returned in listatoms. 
c---------------------------------------------------------------------------
      implicit none
      integer natoms, npaths, iatom, taratom, listatoms
      integer paths(natoms, natoms)
      integer path_count(natoms)
      integer list(natoms)
      logical found_connection

      integer i, j, k

      found_connection = .false.

      do j = 1, npaths
         do i = 1, path_count(j)
            if (paths(i,j) .eq. taratom) then

c---------------------------------------------------------------------------
c   We found the connection.
c---------------------------------------------------------------------------

               found_connection = .true.
               listatoms = path_count(j)
               do k = 1, listatoms
                  list(k) = paths(k,j)
               enddo
               return
            endif
         enddo
      enddo

      return
      end
