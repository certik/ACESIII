
      integer mx_stat_message
      parameter (mx_stat_message = 100)
      integer mx_stat_keys
      parameter (mx_stat_keys = 2000)
c---------------------------------------------------------------------------
c   Server performance statistics
c---------------------------------------------------------------------------
      double precision 
     *                 sstat_tprep(mx_stat_keys), 
     *                 sstat_tprep2(mx_stat_keys), 
     *                 sstat_tprepsum(mx_stat_keys), 
     *                 sstat_tprepsum2(mx_stat_keys), 
     *                 sstat_treq(mx_stat_keys), 
     *                 sstat_treq2(mx_stat_keys), 
     *                 sstat_tpreq(mx_stat_keys), 
     *                 sstat_tpreq2(mx_stat_keys), 
     *                 sstat_trestore(mx_stat_keys),
     *                 sstat_trestore2(mx_stat_keys),
     *                 sstat_tbackup(mx_stat_keys), 
     *                 sstat_tbackup2(mx_stat_keys), 
     *                 sstat_msg_time1
      integer sstat_nprep(mx_stat_keys), sstat_nprepsum(mx_stat_keys), 
     *                 sstat_nreq(mx_stat_keys),
     *                 sstat_npreq(mx_stat_keys),
     *                 sstat_nrestore(mx_stat_keys), 
     *                 sstat_nbackup(mx_stat_keys)

      integer stat_key, next_stat_key, lineno(mx_stat_keys)
      logical do_stats

      common /server_stat/sstat_tprep, sstat_tprep2,
     *                 sstat_tprepsum, sstat_tprepsum2,
     *                 sstat_treq, sstat_treq2, 
     *                 sstat_tpreq, sstat_tpreq2, 
     *                 sstat_trestore, sstat_trestore2,
     *                 sstat_tbackup, sstat_tbackup2,
     *                 sstat_nprep, sstat_nprepsum,
     *                 sstat_nreq, sstat_npreq,
     *                 sstat_nrestore, sstat_nbackup,
     *                 sstat_msg_time1(mx_stat_message),
     *                 stat_key, next_stat_key, lineno, do_stats 
