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
      subroutine push_clean_block(iblock, server_table, nserver_table)
c---------------------------------------------------------------------------
c   Push a block number onto the clean block stack.
c---------------------------------------------------------------------------
      implicit none
      include 'server.h'
      include 'parallel_info.h'

      integer i,j, ptr, iblock
      integer nserver_table
      integer server_table(lserver_table_entry,nserver_table)

      do i = 1, clean_block_ptr
         if (iblock .eq. clean_blocks(i)) then
            return  ! no need to add block
         endif
      enddo

      clean_block_ptr = clean_block_ptr + 1
      if (clean_block_ptr .gt. nserver_memblocks) then
         print *,'Error: clean block stack overflow '
         print *,'nserver_memblocks ',nserver_memblocks
         print *,'mx_server_memblocks ',mx_server_memblocks
         do i = 1, nserver_memblocks
            print *,'entry ',i,' clean_blocks(i) ',clean_blocks(i)
         enddo
         call server_abort_job(server_table, nserver_table)
      endif

      clean_blocks(clean_block_ptr) = iblock 
      nclean_blocks = nclean_blocks + 1
      return
      end
