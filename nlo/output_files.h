c------------------------------------------------------------------------
c  include file   : output_files
c  for module     : output files handler
c  description    : This include file contains the common block for
c                   the host system interface module.
c
c                         MAX_OUTPUT_FILES  =  the number of available
c                                              i/o units.
c
c                                free,used  =  status variables
c
c  author         : Norbert Flocke
c------------------------------------------------------------------------
c
c
c             ...declare variables.
c
c
         logical       logfile_exists

         character*13  logfile_name

         integer       free
         integer       used
         integer       used_by_ofh
         integer       MAX_OUTPUT_FILES

         parameter     (free             = 0              )
         parameter     (used             = 1              )
         parameter     (used_by_ofh      = 2              ) 
         parameter     (logfile_name     = 'aces3.logfile')
         parameter     (MAX_OUTPUT_FILES = 99             )

         integer       file_handle (1:MAX_OUTPUT_FILES)
c
c
c             ...define common blocks.
c
c
         common  /output_files/ file_handle
c
c
c------------------------------------------------------------------------
