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
      subroutine bluegene_assign_companies(hosts, nprocs, 
     *                  nproc_per_node,
     *                  company_id, platoon_id)
      implicit none
      include 'company.h'
   
      integer nprocs, nproc_per_node
      integer hosts(nprocs), company_id(nprocs), platoon_id(nprocs)
      integer nodes, ionodes, nworker, nserver, nserver_per_ionode
      integer node, proc, proc_skip

      integer i, j

      do i = 1, nprocs
         company_id(i) = 0
         platoon_id(i) = 0
      enddo
 
      nodes = nprocs / nproc_per_node
      if (nprocs .ne. nodes * nproc_per_node) nodes = nodes + 1
      ionodes = nodes / 64
      if (nodes .ne. ionodes * 64) ionodes = ionodes + 1

      nworker = c_table(1,c_nwrkr)
      nserver = c_table(2,c_nwrkr)
      nserver_per_ionode = nserver / ionodes
      proc_skip = (64*nproc_per_node) / nserver_per_ionode 
      if (nserver .eq. 0) go to 1000

c---------------------------------------------------------------------------
c   Distribute the IOCOMPANY processes evenly among the I/O nodes.
c   On BLUEGENE, each group of 64 compute nodes is mapped to 1 I/O node.
c---------------------------------------------------------------------------      
      do i = 1, ionodes
         proc = (i-1)*64*nproc_per_node + 2
         do j = 1, nserver_per_ionode
            if (proc .gt. nprocs) then
               print *,'PROC overflow'
               print *,'i, j, proc ',i,j,proc
               print *,'nserver, nprocs, ionodes, nserver_per_ionode ',
     *           nserver, nprocs, ionodes, nserver_per_ionode
               print *,'proc_skip ',proc_skip
               call abort_job()
            endif

            company_id(proc) = 2
            platoon_id(proc) = 1 
            proc = proc + proc_skip
         enddo
      enddo
 
1000  continue
      do i = 1, nprocs
         if (company_id(i) .eq. 0) company_id(i) = 1
      enddo

      return
      end
