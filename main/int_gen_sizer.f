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

c This routine calculates the number of integrals in a shell quadruplet.


      subroutine int_gen_sizer(blk_m,blk_n,blk_r,blk_s,nInts)
      implicit none

      integer blk_m, blk_n, blk_r, blk_s, nInts

      include 'int_gen_parms.h'

      integer mn, rs
      integer full_m, full_n, full_r, full_s

c ----------------------------------------------------------------------

      full_m = nFpS(blk_m)
      full_n = nFpS(blk_n)
      full_r = nFpS(blk_r)
      full_s = nFpS(blk_s)

      if (restrict_mn_rs) then
         if (restrict_mn .and. .not. restrict_rs) then
            stop 'Assertion failed: (SR|N<=M); SR<=NM'
         endif

         if (restrict_rs .and. .not. restrict_mn) then
            stop 'Assertion failed: (S<=R|NM); SR<=NM'
         endif

         if (restrict_mn) then
            if ((blk_m.ne.blk_r).or.(blk_n.ne.blk_s)) then
               if (blk_m.ne.blk_n) then
                  mn = full_n*full_m
               else
                  mn = full_m*(full_m+1)/2
               end if

               if (blk_r.ne.blk_s) then
                  nInts = full_s*full_r*mn
               else
                  nInts = full_r*(full_r+1)/2*mn
               end if
            else
               if (blk_m.ne.blk_n) then
                  mn = full_n*full_m
               else
                  mn = full_m*(full_m+1)/2
               end if
               nInts = mn*(mn+1)/2
            end if
         else
            if ((blk_m.ne.blk_r).or.(blk_n.ne.blk_s)) then
               nInts = full_s*full_r*full_n*full_m
            else
               mn = full_n*full_m
               nInts = mn*(mn+1)/2
            end if
         endif

         return
      endif        ! restrict_mn_rs

c   o (M<=N|
      if (restrict_mn) then
         if (blk_m.ne.blk_n) then
            mn = full_n*full_m
         else
            mn = full_m*(full_m+1)/2
         end if
      else
c   o (MN|
         mn = full_n*full_m
      endif

c   o |R<=S)
      if (restrict_rs) then
         if (blk_r.ne.blk_s) then
            rs = full_s*full_r
         else
            rs = full_r*(full_r+1)/2
         end if
      else
c   o |RS)
         rs = full_s*full_r
      endif

c      print *,'%%%INT_GEN_SIZER: m,n,r,s = ',blk_m,blk_n,
c     *   blk_r, blk_s,' rs, mn, nInts = ',rs, mn, nInts
      nInts = rs*mn

c ----------------------------------------------------------------------

      return
c     end subroutine int_gen_sizer
      end

