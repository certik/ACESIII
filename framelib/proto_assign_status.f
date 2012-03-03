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
      subroutine proto_assign_status(master, ihosts, nprocs, 
     *                               task_status,
     *                               master_is_worker,
     *                               managers_are_workers)

c------------------------------------------------------------------------
c   Determines task assignments based on host IDs and user parameters.
c   
c   Arguments:
c       master  	Global index of the master task.
c	ihosts		Array of processor host IDs.
c	nprocs		Number of processors.
c	task_status	Return array for the status of each task.
c
c   The task status takes on one of the following values:
c	worker_status		the task is a worker.
c	manager_status		the task is a manager.
c	both_status		the task is both worker and manager.
c       master_status           the task is the master.
c       master_worker_status    the task is both master and worker.
c
c-------------------------------------------------------------------------

      implicit none
      include 'proto_defines.h'

      integer master, nprocs
      integer ihosts(nprocs), task_status(nprocs)
      logical master_is_worker, managers_are_workers

      integer i, k
      logical manager_was_assigned

c------------------------------------------------------------------------
c   Initialize status array.
c------------------------------------------------------------------------

      do i = 1, nprocs
         task_status(i) = unassigned_status
      enddo

c-------------------------------------------------------------------------
c   Determine a node index for each task.
c-------------------------------------------------------------------------

      do i = 1, nprocs
         if (task_status(i) .ne. unassigned_status) go to 100

         if (i .eq. master+1) then

c-------------------------------------------------------------------------
c   Assign the master task.
c-------------------------------------------------------------------------

            manager_was_assigned = .false.
            if (master_is_worker) then
               task_status(i) = master_worker_status
            else
               task_status(i) = master_status
            endif
         else

c-------------------------------------------------------------------------
c   Assign this task as manager for its node.
c-------------------------------------------------------------------------

            manager_was_assigned = .true.
            if (managers_are_workers) then
               task_status(i) = both_status
            else
               task_status(i) = manager_status
            endif
         endif

c-------------------------------------------------------------------------
c   Assign all remaining tasks for this host.
c-------------------------------------------------------------------------

         do k = i+1, nprocs
            if (ihosts(k) .eq. ihosts(i)) then

c-------------------------------------------------------------------------
c   If the node's manager was not assigned prior to the loop, it must be 
c   done now.  Otherwise, we just have a worker.
c-------------------------------------------------------------------------
 
               if (.not. manager_was_assigned) then
                  manager_was_assigned = .true.
                  if (managers_are_workers) then
                     task_status(k) = both_status
                  else
                     task_status(k) = manager_status
                  endif
               else 
                  task_status(k) = worker_status 
               endif
            endif
         enddo
   
c---------------------------------------------------------------------------
c   IMPORTANT: If we get here and still haven't assigned a manager, we 
c              must do so now, even if the user has forbidden the use
c              of manager tasks as workers.  The only exception is if
c              the current task is the master, which can never be a 
c              manager.
c---------------------------------------------------------------------------

         if (.not. manager_was_assigned .and.
     *       i .ne. master+1) then
            task_status(i) = both_status       ! must do double duty.
         endif

  100 continue
      enddo

      return
      end

