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
      integer function get_pardo_master()
c----------------------------------------------------------------------------
c   Returns the "pardo master" of the current rank.
c----------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'
      integer pmaster
      integer get_pardo_cluster_size 
      integer cluster_size
 
      cluster_size = get_pardo_cluster_size()
      pmaster = (my_company_rank / cluster_size)*cluster_size
      get_pardo_master = pmaster
      return
      end 

      integer function get_pardo_cluster_size()
c---------------------------------------------------------------------------
c   Returns the pardo_cluster_size
c---------------------------------------------------------------------------
      implicit none
      include 'parallel_info.h'

      integer cluster_size
      parameter (cluster_size = 16)
      integer my_cluster_size, nclusters 

c      my_cluster_size = cluster_size
      my_cluster_size = my_company_size
      get_pardo_cluster_size = min(my_cluster_size, my_company_size)
      return
      end

