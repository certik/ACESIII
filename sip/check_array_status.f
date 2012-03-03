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
      subroutine check_array_status(array, msgtype, source, 
     *                              array_table)
c----------------------------------------------------------------------------
c   Validates the message type with the array's current status.
c   
c   array     Array to be checked.
c   msgtype   1 = GET, 2 = PUT or PUT +=.
c
c   If the message type is not correct for the currnet status of the array,
c   an error message is printed, and the job is aborted.  Otherwise, the
c   subroutine returns, and processing can continue.
c----------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'parallel_info.h'
      include 'trace.h'

      integer array, msgtype
      integer source
      integer array_table(larray_table_entry, *)
      integer status

      status = array_table(c_array_status,array)
      if (msgtype .eq. -3) then

c--------------------------------------------------------------------------
c   GET message: Array should be either 0 or readonly.
c--------------------------------------------------------------------------

         if (status .eq. 0) then
            array_table(c_array_status,array) = read_only_array_status
            print *,'Task ',me,' line ',current_line,
     *         ' RESET ARRAY ',array,' TO READONLY source ',source
         else if (status .ne. read_only_array_status) then
            print *,'Task ',me,
     *        ' Error in check_array_status for array ',array,
     *        ' line ',current_line,' source ',source
            print *,'Message type is ',msgtype,' but array status is ',
     *              'write_only' 
            call abort_job()
         endif
      else if (msgtype .eq. -5 .or. 
     *         msgtype .eq. -6) then
      
c---------------------------------------------------------------------------
c   PUT or PUT += message: Array must be 0 or write_only.
c---------------------------------------------------------------------------

         if (status .eq. 0) then
            array_table(c_array_status,array) = write_only_array_status
            print *,'Task ',me,' line ',current_line,
     *          ' RESET ARRAY ',array,' TO WRITEONLY source ',
     *          source
         else if (status .ne. write_only_array_status) then
            print *,'Task ',me,
     *        ' Error in check_array_status for array ',array,
     *        ' line ',current_line,' source ',source
            print *,'Message type is ',msgtype,' but array status is ',
     *              'read_only' 
            call abort_job()
         endif
      else

c---------------------------------------------------------------------------
c   Invalid message type.
c---------------------------------------------------------------------------

         print *,'Task ',me,' check_array_status was called with an ',
     *           'invalid message type'
         print *,'msgtype = ',msgtype
         call abort_job()
      endif
 
      return
      end
