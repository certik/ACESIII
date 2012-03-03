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
      logical function check_where_conditions(ind, seg, nind, 
     *                               index_table, nindex_table)
c---------------------------------------------------------------------------
c   Checks to see whether a tuple of index-segment tuples satisfies a 
c   set of WHERE conditions.  Returns .true. if the conditons are
c   satisfied, .false. otherwise.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'where_table.h'
      include 'parallel_info.h'

      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer i, j, k, ind1, ind2, nind, cond, seg1, seg2
      integer ind(nind), seg(nind)

      check_where_conditions = .true.
      if (nwhere .le. 0) return

      do 200 i = iwhere, nwhere
         ind1 = where_table(i,c_where_ind1)
         ind2 = where_table(i,c_where_ind2)
         seg1 = index_table(c_current_seg,ind1) 
         seg2 = index_table(c_current_seg,ind2) 
         cond = where_table(i,c_where_cond)

         do j = 1, nind
            if (ind1 .eq. ind(j)) then
               seg1 = seg(j)
         
               do k = 1, nind
                  if (k .ne. j .and. ind2 .eq. ind(k)) then
                     seg2 = seg(k)
                     go to 100
                  endif
               enddo  
            else if (ind2 .eq. ind(j)) then
               seg2 = seg(j)

               do k = 1, nind
                  if (k .ne. j .and. ind1 .eq. ind(k)) then
                     seg1 = seg(k)
                     go to 100
                  endif
               enddo 
            endif    
         enddo   ! j
  100    continue

         if (seg1 .eq. -90909 .or.
     *       seg2 .eq. -90909) then
            check_where_conditions = .false.
            go to 200   ! undefined index.
         endif

         if (cond .eq. 1) then        ! ==
            if (seg1 .ne. seg2) 
     *          check_where_conditions = .false.
         else if (cond .eq. 2) then   ! >=
            if (seg1 .lt. seg2) 
     *          check_where_conditions = .false.
         else if (cond .eq. 3) then   ! <=
            if (seg1 .gt. seg2) 
     *          check_where_conditions = .false.
         else if (cond .eq. 4) then   !  >
            if (seg1 .le. seg2) 
     *          check_where_conditions = .false.
         else if (cond .eq. 5) then   !  <
            if (seg1 .ge. seg2) 
     *          check_where_conditions = .false.
         else if (cond .eq. 6) then   ! !=
            if (seg1 .eq. seg2) 
     *          check_where_conditions = .false.
         endif

         if (.not. check_where_conditions) return
  200 continue

      return
      end


      logical function pcheck_where_conditions(ind, seg, nind, 
     *                               index_table, nindex_table)
c---------------------------------------------------------------------------
c   Checks to see whether a tuple of index-segment tuples satisfies a 
c   set of WHERE conditions.  Returns .true. if the conditons are
c   satisfied, .false. otherwise.
c---------------------------------------------------------------------------
      implicit none
      include 'interpreter.h'
      include 'where_table.h'
      include 'parallel_info.h'

      integer nindex_table
      integer index_table(lindex_table_entry, nindex_table)
      integer i, j, k, ind1, ind2, nind, cond, seg1, seg2
      integer ind(nind), seg(nind)

      pcheck_where_conditions = .true.
      if (nwhere .le. 0) return

      do 200 i = iwhere, nwhere
         ind1 = where_table(i,c_where_ind1)
         ind2 = where_table(i,c_where_ind2)
         seg1 = index_table(c_current_seg,ind1) 
         seg2 = index_table(c_current_seg,ind2) 
         cond = where_table(i,c_where_cond)

         do j = 1, nind
            if (ind1 .eq. ind(j)) then
               seg1 = seg(j)
         
               do k = 1, nind
                  if (k .ne. j .and. ind2 .eq. ind(k)) then
                     seg2 = seg(k)
                     go to 100
                  endif
               enddo  
            else if (ind2 .eq. ind(j)) then
               seg2 = seg(j)

               do k = 1, nind
                  if (k .ne. j .and. ind1 .eq. ind(k)) then
                     seg1 = seg(k)
                     go to 100
                  endif
               enddo 
            endif    
         enddo   ! j
  100    continue

         if (seg1 .eq. -90909 .or.
     *       seg2 .eq. -90909) then
            pcheck_where_conditions = .false.
c            print *,'REJECT UNDEF: seg1, seg2 ',seg1,seg2
c            go to 200   ! undefined index.
         endif

         if (cond .eq. 1) then        ! ==
            if (seg1 .ne. seg2) 
     *          pcheck_where_conditions = .false.
         else if (cond .eq. 2) then   ! >=
            if (seg1 .lt. seg2) 
     *          pcheck_where_conditions = .false.
         else if (cond .eq. 3) then   ! <=
            if (seg1 .gt. seg2) 
     *          pcheck_where_conditions = .false.
         else if (cond .eq. 4) then   !  >
            if (seg1 .le. seg2) 
     *          pcheck_where_conditions = .false.
         else if (cond .eq. 5) then   !  <
            if (seg1 .ge. seg2) 
     *          pcheck_where_conditions = .false.
         else if (cond .eq. 6) then   ! !=
            if (seg1 .eq. seg2) 
     *          pcheck_where_conditions = .false.
         endif

         if (.not. pcheck_where_conditions) then
c             print *,'REJECT COND ',cond,' ind1, ind2, seg1,seg2 ',
c     *           ind1, ind2, seg1,seg2
             return 
         endif
  200 continue

      return
      end


