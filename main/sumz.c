/*
*  Copyright (c) 2003-2010 University of Florida
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  The GNU General Public License is included in this distribution
*  in the file COPYRIGHT.
*/ 
#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <mpi.h>
#include <unistd.h>
#include <sys/errno.h>
#include <sys/types.h>
#include "f77_name.h"
#include "f_types.h"
#include "c_blkmgr.h"
#include "c_timer.h"


#define msgquit -1
#define msg_break -2
#define msggett -3
#define msgdata -5
#define msg_copy_data -6
#define msg_barrier -7
#define msg_master_barrier -8
#define msg_pardo 4450
#define readytag 3456
#define pardo_readytag 4449
#define MAX_THREADS 1
#define MAX_MPI_TESTS 100
#define MAX_PARDO_LOCKS 50000

/*fortran functions definition*/
void one_pass_of_server();
void F77_NAME(sip_fmain, SIP_FMAIN)();
void F77_NAME(get_sip_blocksize, GET_SIP_BLOCKSIZE) (f_int *blocksize,
                f_int *number_blocks);
f_int F77_NAME(pst_get_company_comm, PST_GET_COMPANY_COMM) (f_int *proc);
void F77_NAME(print_opcounter, PRINT_OPCOUNTER)();

void F77_NAME(pardo_loadb_get_my_batch, PARDO_LOADB_GET_MY_BATCH)
             (f_int *iop, f_int *next_batch, f_int *last_batch);
void F77_NAME(send_master_barrier_msg, SEND_MASTER_BARRIER_MSG)();
void F77_NAME(abort_job, ABORT_JOB)();
f_int F77_NAME(get_pardo_master, GET_PARDO_MASTER)();

/*global variables*/
int bsize;
int *csarray;
f_int *array_table;
f_int *block_map_table;

int myrank, nprocs, dbg;
MPI_Comm newcomm;
int c_pardo_lock[MAX_PARDO_LOCKS];
int pardo_list_head, pardo_list_tail;
int barrier_count;
int master_barrier_msg_recvd;

int nmessage_count;      /* Running count of messages (GETs/PUTs) processed */

f_int F77_NAME(query_message_counter, QUERY_MESSAGE_COUNTER)()
{
   f_int n = (f_int) nmessage_count;
   return(n);
}

void F77_NAME(clear_message_counter, CLEAR_MESSAGE_COUNTER)()
{
   nmessage_count = 0;
}

void F77_NAME(set_sumz_tables, SET_SUMZ_TABLES) (f_int *artab, f_int *blkmtab)
{
   array_table = artab;
   block_map_table = blkmtab;
}

void F77_NAME(acquire_pardo_lock, ACQUIRE_PARDO_LOCK)(f_int *my_lock)
{
   if (*my_lock < 0 || *my_lock > MAX_PARDO_LOCKS)
   {
      printf("Task %d: Error in acquire_pardo_lock, lock index %d out of range\n",myrank,*my_lock);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   while (1)
   {
      if (c_pardo_lock[*my_lock-1] == 0)
      {
         c_pardo_lock[*my_lock-1] = 1;
         break;
      }

      one_pass_of_server();
   }
}

void F77_NAME(release_pardo_lock, RELEASE_PARDO_LOCK)(f_int *my_lock)
{
   if (*my_lock < 0 || *my_lock > MAX_PARDO_LOCKS)
   {
      printf("Task %d: Error in release_pardo_lock, lock index %d out of range\n",myrank,*my_lock);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   c_pardo_lock[*my_lock-1] = 0;
}

void F77_NAME(mutex_block, MUTEX_BLOCK) (f_int *iblk)
{
     int blk;

     return;

     if (nprocs == 1) return;   /* No synchronization necessary */

     blk = *iblk - 1;   /* Compensate for Fortran indexing */
     while (1)
     {
        if (csarray[blk] == 0)
        {
           csarray[blk] = 1;
           break;
        }

        one_pass_of_server();
     }

     return;
}

void F77_NAME(release_mutex_block, RELEASE_MUTEX_BLOCK) (f_int *iblk)
{
   int blk;
   int ierr;

   return;

   if (nprocs == 1) return;   /* No mutex necessary */

   blk = (int) (*iblk);
   csarray[blk] = 0;
}

int nreserved_tags = 34;
static int taglist[] = {3456, 3457, 2000, 9990, 9991, 9992, 9993, 9994, 9995,
                        9996, 9997, 9998, 9999, 2237, 2238, 2239, 4444, 4445,
                        3336, 3337,
                        3338, 3335, 3341, 3342, 3346,
                        3347, 3348, 3349, 3350,
                        4446, 4447, 4449, 4450, 4451};
int next_tag = 0;
int tag_limit = 0;
int tag_increment;

int init_msg_tag()
{
   int tag_ub, flag;
   void *v;

#ifdef MPI2
   MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &v, &flag);
   tag_ub = *(int *)v;
#else
   /* MPI_Attr_get(newcomm, MPI_TAG_UB, &tag_ub, &flag); */
   tag_ub = 32767;
#ifdef BLUEGENE 
   tag_ub = 1000000;
#endif
#endif
   if (myrank == 0) printf("MPI TAG_UB value is %d\n",tag_ub);
   tag_limit = tag_ub / 2;

   tag_increment = tag_limit/nprocs;
   next_tag = myrank * tag_increment;
   return 0;
}

int form_msg_tag()
{
   /* Returns a tag for the next message.  The tag is guaranteed not
      to conflict with any in a pre-defined reserved tag list. */

   int valid, i;

   valid = 0;
   while (!valid)
   {
      next_tag++;
      if (next_tag == tag_limit-1) next_tag = myrank*tag_increment + 1;

      valid = 1;
      for (i = 0; i < nreserved_tags; i++)
      {
         if (next_tag == taglist[i] ||
             next_tag+tag_limit == taglist[i])
         {
            valid = 0;
            break;
         }
      }
   }

   /* if (next_tag > tag_limit)
   {
      printf("Task %d Invalid tag %d tag_limit %d\n",myrank,next_tag,tag_limit);
      MPI_Abort(MPI_COMM_WORLD, 1);
   } */
   return next_tag;
}

void F77_NAME(f_usleep, F_USLEEP) (f_int *sleep_interval)
{
   int sleep_time;

   sleep_time = *sleep_interval;
   usleep(sleep_time);
}

f_int F77_NAME(f_form_msg_tag, F_FORM_MSG_TAG)()
{
   static int tag;

   tag = form_msg_tag();

   return tag;
}

f_int F77_NAME(f_form_msg_tag2, F_FORM_MSG_TAG2)()
{
	   static int tag;

	      tag = form_msg_tag() + tag_limit;

	         return tag;
}


void AddingData(double *z, double *buffer, int bsize)
{
    f_int fbsize;
    f_int fone;
    double alpha;
    int i;

    fone = 1;
    fbsize = bsize;
    alpha = 1.0;

    for (i=0; i<bsize; i++)
        z[i]+=buffer[i];

   /* F77_NAME(daxpy, DAXPY) (&fbsize, &alpha, buffer, &fone, z, &fone); */
}

#define MAX_BUFFERS 10000
#define MAX_PARDO_MESSAGES 20000

int quit_count;
int nprocs;
int mbuffer;
int timer_recv;
int timer_send;
int timer_wait;
int timer_sumz;
double *recv_buffers;
int msgbuf[MAX_BUFFERS][5];
MPI_Request ready_request[MAX_BUFFERS];

struct msginfo
{
   int source;
   int tag;
   int nwblock;
   int arrayno;
   int blocknumber;
   int state;
   int msgtype;
   MPI_Request request;
   int block;
   int flag;
   int forced_wait;
   int test_count;
   double *buffer;
   double *block_addr;
   double time;
};

struct msginfo info[MAX_BUFFERS];
int ProcessMsg(struct msginfo *info);

struct pardo_msginfo
{
   int source;
   f_int batch_number;
   f_int last_batch;
   f_int iop;
   f_int lock_index;
   MPI_Request request;
   int state;
   int next_entry;
};

struct pardo_msginfo pardo_info[MAX_PARDO_MESSAGES];
void ProcessPardoMsg(struct pardo_msginfo *info);

#define BeginState 999
#define WaitForCompletion 998
#define NullState 997
#define AddDataState 996

int nactive;
int nactive_pardo;

void init_thread_server()
{
    int i;

    quit_count = 0;
    nmessage_count = 0;
    barrier_count = 0;
    master_barrier_msg_recvd = 0;

    MPI_Comm_size(newcomm, &nprocs);

    nactive = 0;
    nactive_pardo = 0;

    pardo_list_head = -1;
    pardo_list_tail = -1;

    init_msg_tag();   /* Initialize message tag generator */

    for (i = 0; i < mbuffer; i++)
    {
       info[i].state = NullState;
       info[i].buffer = &recv_buffers[i*bsize];
    }

    for (i = 0; i < MAX_PARDO_MESSAGES; i++)
       pardo_info[i].state = NullState;
}

void pardo_add_list_entry(int index)
{
   /* Add an element to the pardo active message list. */

   if (pardo_list_tail == -1)
   {
      /* Empty list.  */

      pardo_list_head = index;
      pardo_list_tail = index;
   }
   else
   {
      /* Add to tail */

      pardo_info[pardo_list_tail].next_entry = index;
      pardo_list_tail = index;
   }
}

void pardo_remove_list_entry(int index)
{
   /* Remove the entry from the list of active pardo messages. */

   int ptr;

   if (pardo_list_head == index)
   {
      pardo_list_head = pardo_info[pardo_list_head].next_entry;
      if (pardo_list_tail == index) pardo_list_tail = pardo_list_head;
      return;
   }
   else
   {
      /* Find the entry that points to "index". */

      ptr = pardo_list_head;
      while (ptr > -1)
      {
         if (pardo_info[ptr].next_entry == index)
         {
            pardo_info[ptr].next_entry = pardo_info[index].next_entry;
            if (pardo_list_tail == index) pardo_list_tail = ptr;
            return;
         }

         ptr = pardo_info[ptr].next_entry;
      }

      printf("Error: Cannot remove entry %d from pardo list. \n",index);
      printf("head %d, tail %d\n",pardo_list_head, pardo_list_tail);
      ptr = pardo_list_head;
      while (ptr > -1)
      {
         printf("Entry %d points to %d\n",ptr,pardo_info[ptr].next_entry);
         ptr = pardo_info[ptr].next_entry;
      }

      MPI_Abort(MPI_COMM_WORLD, 1);
   }
}

void one_pass_of_server()
{
    MPI_Status status;
    MPI_Status status2;
    MPI_Status pardo_status;
    MPI_Status pardo_status2;
    int ierr;
    static int flag;
    static int flag_pardo;
    int i, index, ptr;
    int istat;
    int pardo_message[4];
    int msg[5];
    int pardo_master;

        /* Probe for a new descriptor message */

        flag = 0;
        if (nactive < mbuffer)
           MPI_Iprobe(MPI_ANY_SOURCE, readytag, newcomm, &flag, &status);

        pardo_master = (int) F77_NAME(get_pardo_master, GET_PARDO_MASTER)();
        flag_pardo = 0;
        if (nactive_pardo < MAX_PARDO_MESSAGES && myrank == pardo_master)
        do {  
           MPI_Iprobe(MPI_ANY_SOURCE, pardo_readytag, newcomm, &flag_pardo,
                      &pardo_status);

        /* Receive a pardo message */
        if (flag_pardo)
        {
           for (i = 0; i < MAX_PARDO_MESSAGES; i++)
              if (pardo_info[i].state == NullState)
              {
                 index = i;
                 MPI_Recv(pardo_message, 4, MPI_INT, pardo_status.MPI_SOURCE,
                          pardo_readytag, newcomm, &pardo_status2);
                 break;
              }

           /* Decode the message and begin its processing */
           nactive_pardo++;
           pardo_info[index].source = pardo_status.MPI_SOURCE;
           pardo_info[index].state  = BeginState;
           pardo_info[index].iop    = pardo_message[2];
           pardo_info[index].lock_index = pardo_message[3];
           pardo_info[index].next_entry = -1;
           pardo_add_list_entry(index);

           /* Begin processing this message. */

           ProcessPardoMsg(&(pardo_info[index]));
           if (pardo_info[index].state == NullState)
           {
              pardo_remove_list_entry(index);
              nactive_pardo--;  /* msg is done */
           }
        }
        else
        {
            /* Search for a pardo message to process */

           if (nactive_pardo)
           {
              index = pardo_list_head;
              while (index > -1)
              {
                 ptr = pardo_info[index].next_entry;
                 if (pardo_info[index].state != NullState)
                 {

                    ProcessPardoMsg(&(pardo_info[index]));
                    if (pardo_info[index].state == NullState)
                    {
                       pardo_remove_list_entry(index);
                       nactive_pardo--;
                    }
                 }

                 index = ptr;   /*  Point to next entry */
              }
           }
        }   /* flag_pardo */
        } while (flag_pardo); 

        /* Receive a data message (GET/PUT) */
        if (flag)
        {
           /* Find a free slot and receive the message. */

           for (i = 0; i < mbuffer; i++)
              if (info[i].state == NullState)
              {
                 index = i;
                 MPI_Recv(msg, 5, MPI_INT, status.MPI_SOURCE,
                       readytag, newcomm, &status2);
                 break;
              }

           /* A descriptor message was received for buffer index.
              Decode the message and start its processing.         */

           nactive++;
           info[index].source      = status.MPI_SOURCE;
           info[index].msgtype     = msg[0];
           info[index].arrayno     = msg[1];
           info[index].blocknumber = msg[2];
           info[index].tag         = msg[3];
           info[index].nwblock     = msg[4];
           info[index].state       = BeginState;
           info[index].request     = MPI_REQUEST_NULL;
           info[index].forced_wait = 0;
           info[index].test_count  = 0;
           info[index].time        = MPI_Wtime();

/*           printf("Task %d recv desc msg: source %d type %d tag %d node %d, nactive %d\n",
                   myrank, info[index].source, info[index].msgtype,
                   info[index].tag, index, nactive); */


           if (info[index].msgtype == msggett ||
               info[index].msgtype == msgdata ||
               info[index].msgtype == msg_copy_data)
           {
              info[index].block = c_get_block_number(
                                             info[index].arrayno,
                                             info[index].blocknumber,
                                           array_table, block_map_table);
              if (info[index].block == -1) 
              {
                 printf("Task %d, array %d, block %d   Block not found\n",myrank,info[index].arrayno,info[index].blocknumber);
                 printf("Message type %d\n",info[index].msgtype);
                 F77_NAME(abort_job, ABORT_JOB)();
              }
              info[index].block_addr = c_get_block_data_addr(info[index].arrayno,
                                            info[index].blocknumber, info[index].block);

              if (!info[index].block_addr)
              {
                 printf("Error: Array %d, block %d was requested, but was not found on processor %d\n",
                 info[index].arrayno, info[index].blocknumber, myrank);
                 MPI_Abort(MPI_COMM_WORLD, 1);
              }
           }
        }

        /* Loop over the remaining active messages, checking each of them
           for progress in processing. */

        if (nactive)
        {
           for (index = 0; index < mbuffer; index++)
           {
              if (info[index].state != NullState)
              {
                 istat = ProcessMsg(&(info[index]));
                 if (info[index].state == NullState) nactive--;
              }
            }
         }
}

void server_at_barrier(int code)
{
    MPI_Status status;
    MPI_Status status2;
    MPI_Status pardo_status;
    MPI_Status pardo_status2;
    MPI_Request barrier_request;
    int ierr;
    static int flag;
    static int flag_pardo;
    int i, index;
    int istat;
    int pardo_message[4];
    int msg[5];
    int pardo_master;
    int barrier_messages;
    int barrier_master;
    

    /* Process messages until the following conditions are met:
       1. There are no more messages active in either the msg queue
          or the pardo msg queue.
       2. Probes for remaining messages indicate there are no more
          messages coming from other workers.
       3. A barrier message from each other worker is received by the master.

       When these conditions occur, the master sends an "OK to proceed"
       message (msg_master_barrier) to each other worker.  When these
       messages have completed, all processes may continue processing.
    */

     while (1)
     {
        one_pass_of_server();
        if (nactive_pardo > MAX_PARDO_MESSAGES) MPI_Abort(MPI_COMM_WORLD, 0);

        if (nactive == 0 && nactive_pardo == 0)
        {

		/* Probe for a new descriptor message */

           flag = 0;
           if (nactive < mbuffer)
              MPI_Iprobe(MPI_ANY_SOURCE, readytag, newcomm, &flag, &status);

           flag_pardo = 0;
           pardo_master = (int) F77_NAME(get_pardo_master, GET_PARDO_MASTER)(); 
           if (nactive_pardo < MAX_PARDO_MESSAGES && myrank == pardo_master)
              MPI_Iprobe(MPI_ANY_SOURCE, pardo_readytag, newcomm, &flag_pardo,
                      &pardo_status);

           if (flag == 0 && flag_pardo == 0)
           {
              /* No active messages or pardo messages, and the probe indicates
		 that no messages are outstanding.  If code == 2, we can return,
		 and the calling routine will determine globally whether the net
                 number of sent and recvd messages is 0.  */

              if (code == 2) return;

              /* If code == 1, we must synchronize the master and workers with
                 the barrier arrival messages. */

              if (myrank == 0)
              {
                 barrier_messages = 3;
                 if (nprocs <=3) barrier_messages = nprocs-1;
                 if (barrier_count == barrier_messages + (nprocs-1) / 4)
                 {
                    /* Each worker has "checked in" with the master, indicating
                       that all workers have arrived at the barrier.  The master
                       responds with an acknowledgement signal to each worker. 
                       When these messages have all been received, and each
                       worker is "quiet", we can exit the thread server and
                       check to see if all messages have been processed. */

                     /* printf("Master sending master_barrier_msg to workers\n"); */
                    barrier_count = 0;
                    F77_NAME(send_master_barrier_msg,SEND_MASTER_BARRIER_MSG)();
                    return;
                 }
              }
              else
              {
                 /* Is this a "barrier master" rank? */
                 barrier_master = (myrank/4)*4;
                 if (myrank == barrier_master)
                 {
                    barrier_messages = 3;
                    if (nprocs - barrier_master <= 3) 
                       barrier_messages = nprocs-barrier_master-1;
                    if (barrier_count == barrier_messages) 
                    {
                       msg[0] = msg_barrier;
                       /* printf("Rank %d sent bmaster msg to master\n",myrank); */
                       MPI_Send(msg, 5, MPI_INT, 0, readytag, newcomm);
                       /* printf("Rank %d Completed sending bmaster msg to master\n",myrank); */
                       barrier_count = 0;
                    }
                 }

                 /* We are not on the master process.  Block until receiving
                    the "OK to proceed" msg from the master. */

                 if (master_barrier_msg_recvd == 1)
                 {
                    master_barrier_msg_recvd = 0;
                    return;
                 }
              }
           }

           /* Receive a pardo message */
           if (flag_pardo)
           {
              for (i = 0; i < MAX_PARDO_MESSAGES; i++)
                 if (pardo_info[i].state == NullState)
                 {
                    index = i;
                    MPI_Recv(pardo_message, 4, MPI_INT, pardo_status.MPI_SOURCE,
                          pardo_readytag, newcomm, &pardo_status2);
                    break;
                 }

              /* Decode the message and begin its processing */
              nactive_pardo++;
              pardo_info[index].source = pardo_status.MPI_SOURCE;
              pardo_info[index].state  = BeginState;
              pardo_info[index].iop    = pardo_message[2];
              pardo_info[index].lock_index = pardo_message[3];
              pardo_info[index].next_entry = -1;
              pardo_add_list_entry(index);

              /* Begin processing this message. */
              ProcessPardoMsg(&(pardo_info[index]));
              if (pardo_info[index].state == NullState)
              {
                 pardo_remove_list_entry(index);
		 nactive_pardo--;  /* msg is done */
	      }
           }   /* flag_pardo */

           /* Receive a data message (GET/PUT) */
           if (flag)
           {
              /* Find a free slot and receive the message. */

              for (i = 0; i < mbuffer; i++)
                 if (info[i].state == NullState)
                 {
                    index = i;
                    MPI_Recv(msg, 5, MPI_INT, status.MPI_SOURCE,
                          readytag, newcomm, &status2);
                    break;
                 }

              /* A descriptor message was received for buffer index.
                 Decode the message and start its processing.         */

              nactive++;
              info[index].source      = status.MPI_SOURCE;
              info[index].msgtype     = msg[0];
              info[index].arrayno     = msg[1];
              info[index].blocknumber = msg[2];
              info[index].tag         = msg[3];
              info[index].nwblock     = msg[4];
              info[index].state       = BeginState;
              info[index].request     = MPI_REQUEST_NULL;
              info[index].forced_wait = 0;
              info[index].test_count  = 0;
/*
           printf("Task %d recv desc msg at barrier: source %d type %d tag %d node %d, nactive %d\n",
                   myrank, info[index].source, info[index].msgtype,
                   info[index].tag, index, nactive);
*/

              if (info[index].msgtype == msggett ||
                  info[index].msgtype == msgdata ||
                  info[index].msgtype == msg_copy_data)
              {
                 info[index].block = c_get_block_number(
                                             info[index].arrayno,
                                             info[index].blocknumber,
                                             array_table, block_map_table);
                 info[index].block_addr = c_get_block_data_addr(info[index].arrayno,
                                            info[index].blocknumber, info[index].block);
                 if (!info[index].block_addr)
                 {
                    printf("Error: Array %d, block %d was requested, but was not found on processor %d\n",
                    info[index].arrayno, info[index].blocknumber, myrank);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                 }
              }
           }    /* flag */
        }       /* nactive == 0 && nactive_pardo == 0 */
     }          /* while */
}

int ProcessMsg(struct msginfo *info)
{
   double *block_addr;
   int i, flag, wait_for_block;
   MPI_Status status;
   int ierr;

   /* printf("Task %d processing message type %d, tag = %d, state %d\n", myrank, info->msgtype, info->tag, info->state); */
   if (info->msgtype == msggett)               /* GET message */
   {
      if (info->state == BeginState)
      {
         if (csarray[info->block] == 1) return 0;

         csarray[info->block] = 1;   /* set the lock */
         /* c_timer_start(timer_send); */
         ierr = MPI_Isend(info->block_addr, info->nwblock, MPI_DOUBLE,
                info->source, info->tag, newcomm, &(info->request));
         /* c_update_timer(timer_send); */

         info->state = WaitForCompletion;
         return 0;
      }
      else if (info->state == WaitForCompletion)
      {
         wait_for_block = info->forced_wait;
         info->test_count++;
#ifdef HP
         if (info->test_count > MAX_MPI_TESTS) wait_for_block = 1;
#endif

         if (wait_for_block)
         {
            ierr = MPI_Wait(&(info->request), &status);
            info->flag = 1;
         }
         else
            ierr = MPI_Test(&(info->request), &info->flag, &status);

         if (ierr)
         {
            printf("Task %d MPI Error: %d\n",myrank,ierr);
            MPI_Abort(MPI_COMM_WORLD,0);
         }

         if (info->flag)
         {
            csarray[info->block] = 0;   /* release the lock */
            info->state = NullState;
            nmessage_count--;
         }

         return 0;
      }
   }
   else if (info->msgtype == msgdata)          /* PUT += */
   {
       /* PUT_increment: Recv the data into a buffer, then sum it with
          the data existing in the blkmgr's memory. */

      if (info->state == BeginState)
      {
         /* c_timer_start(timer_recv); */
         ierr = MPI_Irecv(info->buffer, info->nwblock, MPI_DOUBLE,
                info->source, info->tag, newcomm, &(info->request));
         /* c_update_timer(timer_recv); */
         info->state = WaitForCompletion;
         return 0;
      }
      else if (info->state == WaitForCompletion)
      {
         wait_for_block = info->forced_wait;
         info->test_count++;
#ifdef HP
         if (info->test_count > MAX_MPI_TESTS) wait_for_block = 1;
#endif

         if (wait_for_block)
         {
            ierr = MPI_Wait(&(info->request), &status);
            info->flag = 1;
         }
         else
            ierr = MPI_Test(&(info->request), &info->flag, &status);

         if (ierr)
         {
            printf("Task %d MPI Error #2: %d\n",myrank,ierr);
            MPI_Abort(MPI_COMM_WORLD,0);
         }

         if (info->flag)
         {
           info->state = AddDataState;
         }
         return 0;
      }
      else if (info->state == AddDataState)
      {
         if (csarray[info->block] == 1) return 0;

         /* c_timer_start(timer_sumz); */
         AddingData(info->block_addr, info->buffer, info->nwblock);
         /* c_update_timer(timer_sumz); */

         info->state = NullState;   /* All done with this request. */
         nmessage_count--;
         return 0;
      }
   }
   else if (info->msgtype == msg_copy_data)    /* PUT */
   {
       /* PUT: Recv the data into a buffer, then store it in
          the blkmgr's memory. */

      if (info->state == BeginState)
      {
         if (csarray[info->block] == 1) return 0;

         csarray[info->block] = 1;   /* set the lock */
         /* c_timer_start(timer_recv); */
         MPI_Irecv(info->block_addr, info->nwblock, MPI_DOUBLE,
                info->source, info->tag, newcomm, &(info->request));
         /* c_update_timer(timer_recv); */

         info->state = WaitForCompletion;
         return 0;
      }
      else if (info->state == WaitForCompletion)
      {
         wait_for_block = info->forced_wait;
         info->test_count++;
#ifdef HP
         if (info->test_count > MAX_MPI_TESTS) wait_for_block = 1;
#endif

         if (wait_for_block)
         {
            ierr = MPI_Wait(&(info->request), &status);
            info->flag = 1;
         }
         else
            ierr = MPI_Test(&(info->request), &info->flag, &status);

         if (ierr)
         {
            printf("Task %d MPI error #3: %d\n",myrank,ierr);
            MPI_Abort(MPI_COMM_WORLD,0);
         }

         if (info->flag)
         {
            csarray[info->block] = 0;   /* release the lock */
            info->state = NullState;
            nmessage_count--;
         }
         return 0;
      }
   }
   else if (info->msgtype == msg_barrier)      /* barrier wait */
   {
      barrier_count++;
      info->state = NullState;
      return 0;
   }
   else if (info->msgtype == msg_master_barrier)  /* master OK to proceed */
   {
      master_barrier_msg_recvd = 1;
      info->state = NullState;
   }
   else if (info->msgtype == msgquit)         /* Quit job */
   {
      if (info->state == BeginState)
      {
         quit_count++;
         info->state = NullState;

         if (quit_count == nprocs - 1)
            return 1;   /* Shut down the thread. */
         else
            return 0;   /* Looking for more "quit" messages. */
      }
   }
   else
   {
      printf("Task %d Invalid message (%d) received by thread server.\n",
             myrank, info->msgtype);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   return 0;
}

void ProcessPardoMsg(struct pardo_msginfo *info)
{
   MPI_Status status;
   int ierr;
   f_int my_lock;
   int flag;

   /* This is a "load-balancing pardo" request for the next batch of work
         to process.  */

   if (info->state == BeginState)
   {
      my_lock = info->lock_index;
      if (c_pardo_lock[my_lock-1] == 0)
      {
         c_pardo_lock[my_lock-1] = 1;
         F77_NAME(pardo_loadb_get_my_batch, PARDO_LOADB_GET_MY_BATCH)(&info->iop, &info->batch_number, &info->last_batch);
         c_pardo_lock[my_lock-1] = 0;

         /* Send the batch number back */

         MPI_Isend(&info->batch_number, 2, MPI_INT,
             info->source, msg_pardo, newcomm, &(info->request));

         /* Next state is "WaitForCompletion" */
         info->state = WaitForCompletion;
      }
   }
   else if (info->state == WaitForCompletion)
   {
      MPI_Test(&(info->request), &flag, &status);
      if (flag) info->state = NullState;
   }

   return;
}

void sendquit(int target, MPI_Comm *comm)
{
    int temp[5];
    temp[0]=msgquit;
    temp[1]=target;
    if (myrank!=target)
    {
        MPI_Send(temp, 5, MPI_INT, target, readytag, *comm);
    }
}

void F77_NAME(fsendquit, FSENDQUIT)(f_int *target)
{
   int mytarget;

   mytarget = *target;
   sendquit(mytarget, &newcomm);
}

void init_available_msg_buf()
{
   int i;

   for (i = 0; i < MAX_BUFFERS; i++)
      ready_request[i] = MPI_REQUEST_NULL;
}

int find_available_msg_buffer()
{
   int i, flag;
   MPI_Status status;

   while (1)
   {
      for (i = 0; i < MAX_BUFFERS; i++)
      {
         if (ready_request[i] == MPI_REQUEST_NULL)
         {
            return i;
         }
         else
         {
            MPI_Test(&ready_request[i], &flag, &status);
            if (flag) return i;
         }
      }
   }
}

void sendrealdata(int arrayno, int blocknumber, int blkndx,
                  int target, int tarray,
                  int tblock,int nwblock, MPI_Comm *comm, int msgtype,
                  MPI_Request *request)
{

    /* Sends the data in arrayno, blocknumber to the target processor's tarray, tblock.  The data will be summed into tarray, tblock upon arrival. */

    int i, tag;
    static double *buffer;
    MPI_Status status;
    double d;

    buffer = c_get_block_data_addr(arrayno, blocknumber, blkndx);
    if (buffer)
    {

       i = find_available_msg_buffer();
       msgbuf[i][0]=msgtype;
       msgbuf[i][1]=tarray;
       msgbuf[i][2]=tblock;
       tag = form_msg_tag() + tag_limit;
       msgbuf[i][3] = tag;
       msgbuf[i][4] = nwblock;

       MPI_Isend(buffer, nwblock, MPI_DOUBLE, target, tag, *comm, request);
       MPI_Isend(msgbuf[i], 5, MPI_INT, target, readytag, *comm, &ready_request[i]);
       nmessage_count++;
    }
    else
    {
       printf("Error: sendrealdata cannot find array, block %d %d in blkmgr.\n", arrayno, blocknumber);
       MPI_Abort(MPI_COMM_WORLD, 0);
    }
}


void F77_NAME(fsumdata, FSUMDATA) (f_int *arrayno, f_int *block, f_int *blkndx,
              f_int *target,
              f_int *tarray, f_int *tarray_block, f_int *nwblock, f_int *fcomm,
              f_int *frequest)
{
    /*Fortran-callable routine to send one block to a particular target
      processor. */
   int array, blocknum, iblk, mytarget, my_nwblock;
   int mytarray, mytarray_block;
   MPI_Comm comm;
   MPI_Request request;

   array = *arrayno;
   blocknum = *block;
   iblk     = *blkndx - 1;   /* Adjust for FORTRAN indexing */
   mytarget = *target;
   mytarray = *tarray;
   mytarray_block = *tarray_block;
   my_nwblock = *nwblock;
#ifdef MPIF2C
   comm = MPI_Comm_f2c( (MPI_Fint)(*fcomm) );
#else
   comm     = (MPI_Comm) (*fcomm);
#endif

   sendrealdata(array, blocknum, iblk, mytarget, mytarray, mytarray_block,
                my_nwblock, &comm, msgdata, &request);
#ifdef MPIF2C
   *frequest = (f_int) MPI_Request_c2f(request);
#else
   *frequest = (f_int) request;
#endif
}

void F77_NAME(fcopydata, FCOPYDATA) (f_int *arrayno, f_int *block, f_int *blkndx,
              f_int *target, f_int *tarray, f_int *tarray_block,
              f_int *nwblock, f_int *fcomm, f_int *frequest)
{
    /*Fortran-callable routine to send one block to a particular target
      processor. */
   int array, blocknum, iblk, mytarget, my_nwblock;
   int mytarray, mytarray_block;
   MPI_Comm comm;
   MPI_Request request;

   array = *arrayno;
   blocknum = *block;
   iblk     = *blkndx - 1;  /* Compensate for FORTRAN indexing of blocks */
   mytarget = *target;
   mytarray = *tarray;
   mytarray_block = *tarray_block;
   my_nwblock = *nwblock;
#ifdef MPIF2C
   comm = MPI_Comm_f2c( (MPI_Fint)(*fcomm) );
#else
   comm     = (MPI_Comm) (*fcomm);
#endif

   sendrealdata(array, blocknum, iblk, mytarget, mytarray, mytarray_block,
                my_nwblock, &comm, msg_copy_data, &request);
#ifdef MPIF2C
   *frequest = (f_int) MPI_Request_c2f(request);
#else
   *frequest = (f_int) request;
#endif
}

void requestblock(int arrayno, int blocknumber, int blkndx, int size,
           int source,
           MPI_Request *request, MPI_Comm *comm)
{
    static int tag;
    int i, ierr;
    double *t;
    MPI_Status status;

    if (source==myrank)
    {
        *request=MPI_REQUEST_NULL;
        return;
    }

    /* Post a recv for the block. */

    t = c_get_block_data_addr(arrayno, blocknumber, blkndx);
    if (!t)
    {
       printf("Error: requestblock is unable to find arrayno, blocknumber %d %d/n", arrayno, blocknumber);
       MPI_Abort(MPI_COMM_WORLD, 0);
    }

    tag = form_msg_tag();

    MPI_Irecv(t, size, MPI_DOUBLE, source, tag, *comm, request);

    /* Now send a msg to the thread server on the other processor,
       requesting that he send us the block. */

    i = find_available_msg_buffer();
    msgbuf[i][0]=msggett;
    msgbuf[i][1]=arrayno;
    msgbuf[i][2]=blocknumber;
    msgbuf[i][3] = tag;
    msgbuf[i][4] = size;
    MPI_Isend(msgbuf[i], 5, MPI_INT, source, readytag, *comm,
              &ready_request[i]);
    nmessage_count++;
}

void F77_NAME(frequestblk, FREQUESTBLK) (f_int *array, f_int *source,
        f_int *block, f_int *blkndx, f_int *size, f_int *request, f_int *fcomm)
{
    /*Fortran-callable routine to request one block from a particular source
      processor. */

    int arrayno;
    int blocknumber;
    int iblk;
    int mysize;
    int mysource;
    MPI_Request myrequest;
    MPI_Comm comm;

    arrayno     = *array;
    blocknumber = *block;
    iblk        = *blkndx - 1; /* Account for FORTRAN indexing of blocks */
    mysize      = *size;
    mysource    = *source;
#ifdef MPIF2C
   comm = MPI_Comm_f2c( (MPI_Fint)(*fcomm) );
#else
   comm     = (MPI_Comm) (*fcomm);
#endif

    requestblock(arrayno, blocknumber, iblk, mysize, mysource, &myrequest, &comm);
#ifdef MPIF2C
   *request = (f_int) MPI_Request_c2f(myrequest);
#else
   *request = (f_int) myrequest;
#endif
}

void F77_NAME(exec_thread_server, EXEC_THREAD_SERVER)(f_int *bflag)
{
   int barrier_flag = (int) (*bflag);

   if (barrier_flag)
   {
      server_at_barrier(barrier_flag);
   }
   else
   {
      one_pass_of_server();
   }
}

void F77_NAME(sumz_work, SUMZ_WORK) (f_int *dryrun_flag, f_int *fmbuffer,
                                     f_int *dbg_flag, double *totalrecvbuffer)
{

    int mworkingbuffer;
    int i, j;
    int provide_level;
    int me, nblocks;
    int *arg;
    int fcomm;

    f_int blocksize, number_blocks;
    char tstr[41];

    if (*dryrun_flag ) return;   /* Nothing to do for dryrun */

    mbuffer = *fmbuffer;
    mworkingbuffer = mbuffer;
    recv_buffers   = totalrecvbuffer;
    dbg            = *dbg_flag;

    init_available_msg_buf();

    /* Get the company communicator for the current processor's company. */

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    {
       f_int me2=me;
#ifdef MPIF2C
       fcomm = F77_NAME(pst_get_company_comm, PST_GET_COMPANY_COMM)(&me2);
       newcomm = MPI_Comm_f2c( (MPI_Fint)(fcomm) );
#else
       newcomm = (MPI_Comm) F77_NAME(pst_get_company_comm, PST_GET_COMPANY_COMM) (&me2);
#endif
    }

    MPI_Comm_rank(newcomm, &myrank);
    MPI_Comm_size(newcomm, &nprocs);
    if (nprocs == 1)
    {
       /* Only one processor, no need for thread servers */

       F77_NAME(sip_fmain, SIP_FMAIN) ();
       return;
    }

    F77_NAME(get_sip_blocksize, GET_SIP_BLOCKSIZE) (&blocksize,
                 &number_blocks);
    bsize = blocksize;         /* Standard size of each managed block. */
    nblocks=number_blocks;

    /* Initialize a critical section for each block managed by the blkmgr. */
    csarray=(int *)malloc(sizeof(int)*nblocks);
    if (!csarray)
    {
       printf("Task %d: Unable to malloc %d blocks for csarray.\n",myrank,nblocks);
       MPI_Abort(MPI_COMM_WORLD, 0);
    }

    for (i=0; i<nblocks; i++)
       csarray[i] = 0;

    /* Initialize mutexes for load-balancing pardo loops. */
    for (i = 0; i < MAX_PARDO_LOCKS; i++)
       c_pardo_lock[i] = 0;


   /* Thread-based timers */

/*
   sprintf(tstr, "Thread recv (E)");
   timer_recv = c_register_timer(tstr, ELAPSED_TIME_TIMER);
   sprintf(tstr, "Thread send (E)");
   timer_send = c_register_timer(tstr, ELAPSED_TIME_TIMER);
   sprintf(tstr, "Thread wait (E)");
   timer_wait = c_register_timer(tstr, ELAPSED_TIME_TIMER);
   sprintf(tstr, "Thread sum  (C)");
   timer_sumz = c_register_timer(tstr, ELAPSED_TIME_TIMER);
*/

   init_thread_server();

   if (dbg) printf("Task %d: Thread server was initialized with %d buffers, size %d.\n",
          myrank,mbuffer,bsize);

   /* Enter computational execution of SIAL code */
   F77_NAME(sip_fmain, SIP_FMAIN) ();

   /* Clean up memory */

   free(csarray);

    if (dbg) printf("%d End of sumz\n", myrank);
    return;
}
