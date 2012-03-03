c----------------------------------------------------------------------------
c   Definition of entries on BLOCK_INDEX file.
c----------------------------------------------------------------------------

      integer c_blk_index_array
      parameter (c_blk_index_array = 1)
      integer  c_blk_index_blkno
      parameter (c_blk_index_blkno = 2)
      integer c_blk_index_size
      parameter (c_blk_index_size = 3)
      integer c_blk_index_nind
      parameter (c_blk_index_nind = 4)
      integer c_blk_index_bsegs
      parameter ( c_blk_index_bsegs = 5)
      integer c_blk_index_esegs
      parameter (c_blk_index_esegs = c_blk_index_bsegs+mx_array_index)
      integer lblk_index_entry
      parameter (lblk_index_entry = c_blk_index_esegs+mx_array_index-1)

