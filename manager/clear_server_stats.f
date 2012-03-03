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
      subroutine clear_server_stats()
c--------------------------------------------------------------------------
c   Initialize server statistics to 0.
c--------------------------------------------------------------------------
      implicit none
      include 'server_stat.h'
      include 'mpif.h'

      integer i

      stat_key = 1
      next_stat_key = 0

      do i = 1, mx_stat_keys
         sstat_tprep(i) = 0.
         sstat_tprepsum(i) = 0.
         sstat_treq(i)    = 0.
         sstat_tpreq(i)    = 0.
         sstat_trestore(i) = 0.
         sstat_tbackup(i)  = 0.

         sstat_tprep2(i) = 0.
         sstat_tprepsum2(i) = 0.
         sstat_treq2(i)    = 0.
         sstat_tpreq2(i)    = 0.
         sstat_trestore2(i) = 0.
         sstat_tbackup2(i)  = 0.

         sstat_nprep(i) = 0
         sstat_nprepsum(i) = 0
         sstat_nreq(i)     = 0
         sstat_npreq(i)     = 0
         sstat_nrestore(i) = 0
         sstat_nbackup(i)  = 0
      enddo

      return
      end
