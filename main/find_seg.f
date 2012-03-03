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
      subroutine find_seg(nocci,nvirti,norbi,nproci,o_seg,v_seg,n_seg) 
c
      Implicit none 
      integer nocci, nvirti, norbi, nproci  
      real*8 nocc, nvirt, norb, nproc  

      integer k, nmin, nmax, nb, wb(1000)  
      integer nseg, neff, rem, best, rem_min   
      integer o_seg,v_seg,n_seg 

      nocc = nocci
      nvirt = nvirti
      norb = norbi
      nproc = nproci
c
c norb segmentation
c ----------------- 
c
      nmin = 20 
      nmax = 35 

c----------------------------------------------------------------------------
c   Handle tiny jobs appropriately.
c----------------------------------------------------------------------------

      if (norb .le. nmin) then
         o_seg = norb
         v_seg = norb
         n_seg = norb
         return
      endif
c
c make sure range is applicable 
c ----------------------------- 
c
1     continue 
      if (nmin .ge. norb) then 
         nmin = norb/2 
         go to 1  
         else 
      endif 
c
2     continue 
      if (nmax .ge. norb) then 
         nmax = norb-1  
         go to 2  
         else 
      endif 
c 
      rem_min = nmax 
      best = nmax 
      nb = 0 
c
      do k = nmin, nmax 
c
         nseg = norb/k 
         neff = nseg*k 
         rem = norb - neff 
         if (rem .le. rem_min) then 
            rem_min = rem 
         endif 
c
      enddo 
c
      do k = nmax, nmin, -1  
c
         nseg = norb/k 
         neff = nseg*k 
         rem = norb - neff 
         if (rem .le. rem_min) then 
            nb = nb + 1 
            wb(nb) = k 
            rem_min = rem 
            go to 20 
         endif 
c
      enddo 
20    continue 
c
c ckeck if you like the result 
c ---------------------------- 
c
      do k = wb(1), wb(1) + nmax 
         nseg = norb/wb(1)  
         neff = nseg*k 
         rem = norb - neff  
         if (rem .le. 0) then 
            wb(1) = k 
            n_seg = k 
            go to 200 
         endif 
      enddo 

200   continue 
c
c nvirt segmentation
c ------------------ 
c
      nmin    = 22 
      nmax    = 35 
c
c make sure range is applicable 
c ----------------------------- 
c
3     continue 
      if (nmin .ge. nvirt) then 
         nmin = nvirt/2 
         go to 3 
         else 
      endif 
c
4     continue 
      if (nmax .ge. nvirt) then 
         nmax = nvirt-1  
         go to 4  
         else 
      endif 
c 
      rem_min = nmax 
      best    = nmax 
      nb      = 0 
c
      do k = nmin, nmax 
c
         nseg = nvirt/k 
         neff = nseg*k 
         rem  = nvirt - neff 
         if (rem .le. rem_min) then 
            rem_min = rem 
         endif 
c
      enddo 
c
      do k = nmax, nmin, -1  
c
         nseg = nvirt/k 
         neff = nseg*k 
         rem  = nvirt - neff 
         if (rem .le. rem_min) then 
            nb = nb + 1 
            wb(nb) = k 
            rem_min = rem 
            go to 10 
         endif 
c
      enddo 
10    continue 
c
c ckeck if you like the result 
c ---------------------------- 
c
      do k = wb(1), wb(1) + nmax 
         nseg = nvirt/wb(1)  
         neff = nseg*k 
         rem = nvirt - neff  
         if (rem .le. 0) then 
            wb(1) = k 
            v_seg = k 
            go to 100 
         endif 
      enddo 

100   continue 
c
c nocc segmentation
c ------------------ 
c
      nmin    = 10 
      nmax    = 30  
c
c make sure range is applicable 
c ----------------------------- 
c
5     continue 
      if (nmin .ge. nocc) then 
         nmin = nocc/2 
         go to 5 
         else 
      endif 
c
6     continue 
      if (nmax .ge. nocc) then 
         nmax = nocc-1  
         go to 6  
         else 
      endif 
c 
      rem_min = nmax 
      best    = nmax 
      nb      = 0 
c
      do k = nmax, nmin, -1  
c
         nseg = nocc/k 
         neff = nseg*k 
         rem  = nocc - neff 
         if (rem .le. rem_min) then 
            rem_min = rem 
         endif 
c
      enddo 
c
      do k = nmax, nmin, -1  
c
         nseg = nocc/k 
         neff = nseg*k 
         rem  = nocc - neff 
         if (rem .le. rem_min) then 
            nb = nb + 1 
            wb(nb) = k 
            rem_min = rem 
            go to 30 
         endif 
c
      enddo 
30    continue 
c
c ckeck if you like the result 
c ---------------------------- 
c
      do k = wb(1), wb(1) + nmax 
         nseg = nocc/wb(1)  
         neff = nseg*k 
         rem = nocc - neff  
         if (rem .le. 0) then 
            wb(1) = k 
            o_seg = k 
            go to 300 
         endif 
      enddo 

300   continue 
c
c   Attempt to improve segmentation.
c
      call refine_segs(nocci, 20, 35, o_seg)
      call refine_segs(nvirti, 20, 35, v_seg)
      call refine_segs(norbi, 20, 35, n_seg)

      return
      end 

      subroutine refine_segs(n, min_size, max_size, segsize)
c---------------------------------------------------------------------------
c   Attempts to improve the segmentation after the initial estimate.
c   Checks to determine if there is a smaller number of segments with a
c   better fit to the problem size.
c
c   Arguments:
c      n        Number of orbitals in either AO, virtual, or occupied space.
c      min_size Min. allowable segment size.
c      max_size Max. allowable segment size.
c      segsize  Initial segment size.
c
c   If a better fit is determined, the "segsize" argument is overwritten 
c   with the new value.
c---------------------------------------------------------------------------

      implicit none  
      integer n, min_size, max_size, segsize
      integer iseg, iseg1, iseg2, nseg 
      integer misfit, my_misfit
      integer total
      integer isize, my_size

c      print *,'Problem size ',n,' Original segsize ',
c     *           segsize
      nseg = (n +segsize - 1) / segsize
      if (nseg .eq. 1) return   ! pathologically small case
      if (n .lt. max_size) return

      total = nseg*segsize
      misfit = total - n

      iseg1 = (n + max_size - 1) / max_size
      iseg2 = (n + min_size - 1) / min_size
 
      if (misfit .gt. 0) then
         my_size = segsize
         do iseg = iseg1, iseg2

c--------------------------------------------------------------------------
c   Find the smallest segment size for with the smallest misfit for the 
c   problem size.
c--------------------------------------------------------------------------

            do isize = min_size, max_size
               total = isize * iseg
               my_misfit = total - n
               
               if (my_misfit .ge. 0 .and.
     *             my_misfit .lt. misfit) then 
                  misfit = my_misfit
                  my_size = isize
               endif
            enddo

            if (my_size .ne. segsize) then
c               print *,'New segsize ',
c     *            my_size,' for iseg ',iseg,' misfit = ',misfit
               segsize = my_size
               return
            endif
         enddo
      endif

      return
      end
