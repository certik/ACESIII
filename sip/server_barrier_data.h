c-------------------------------------------------------------------------
c   The following common block contains the MPI request info for any 
c   descriptor messages sent to the servers of the I/O company.  The 
c   server_barrier subroutine must wait until all the MPI sends have
c   completed before sending the barrier signal to the servers.  This 
c   guarantees correct sequencing of the descriptor messages and allows
c   each server to correctly complete a set of messages before processing
c   additional data.
c--------------------------------------------------------------------------

      integer mx_msg
      parameter (mx_msg = 10000)
      integer server_requests
      common /server_barrier_data/server_requests(mx_msg)
