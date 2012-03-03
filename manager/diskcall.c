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
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include "f77_name.h"
#include "f_types.h"

#ifdef BLUEGENE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#endif

void F77_NAME(f_creatfile, F_CREATFILE)(char *fn, f_int *fh);
void F77_NAME(f_openfile, F_OPENFILE)(char *fn, f_int *fh);
void F77_NAME(f_deletefile,F_DELETEFILE)(char *fn, f_int *fh);
f_int F77_NAME(f_backupram, F_BACKUPRAM)(int *fhandle, f_int *diskloc,
                     f_int *blksize, double *buffer, f_int *fsize);
void F77_NAME(f_restoreram, F_RESTORERAM) (f_int *fhandle, f_int *diskloc, 
                  f_int *blksize, double *buffer, f_int *fsize);
void F77_NAME(f_close_file, F_CLOSE_FILE) (f_int *fh);
void F77_NAME(f_write_disk, F_WRITE_DISK) (f_int *fh, long long *startpos,
                                   double *buffer, f_int *fsize);
void F77_NAME(f_read_disk, F_READ_DISK) (f_int *fh, long long *startpos,
                                   double *buffer, f_int *fsize);

void F77_NAME(f_close_file, F_CLOSE_FILE) (f_int *fh)
{
   int handle = *fh;
   if (!handle) close(handle);
}

int backupram(int fh, long long startpos, double* buffer, 
               long size);

int restoreram(int fh, long long startpos, double* buffer, 
               long size);


int creatfile (char* filename, int *fh)
{
    int current=0;
    int filenumber=1;

    if( (*fh = open( (const char*)filename, O_CREAT|O_RDWR, S_IRUSR | 
                            S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
        printf( "failed on creating input file %s, errno=%d, EEXIST=%d\n", 
            filename, errno, EEXIST );
        return -1;
    }
    return 1;
}

int openfile (char* filename, int *fh)
{
    int current=0;
    int filenumber=1;

    if( (*fh = open( (const char*)filename, O_RDWR, S_IRUSR | 
                            S_IWUSR | S_IRGRP | S_IROTH )) == -1)
    {
        printf( "failed on opening input file %s, errno=%d, EEXIST=%d\n", 
            filename, errno, EEXIST );
        return -1;
    }
    return 1;
}

void F77_NAME(f_creatfile, F_CREATFILE)(char *fn, f_int *fh)
{
   int cfh, ierr;
   ierr = creatfile(fn, &cfh);
   *fh = (f_int) cfh;
}

void F77_NAME(f_openfile, F_OPENFILE)(char* filename, f_int *fh)
{
    int handle;

    if( (handle = open( (const char*)filename, O_RDWR, S_IRUSR | 
                            S_IWUSR)) == -1)
    {
        printf( "failed on opening input file %s, errno=%d, EEXIST=%d\n", 
            filename, errno, EEXIST );
    }

    *fh = handle;
}

/*==========================================================================
copy double data from buffer with the size of size.
return -1 if can not open file
return 1 if successful.
==========================================================================*/
int backupram(int fh, long long startpos, double* buffer, 
              long size)
{
    int me;

#ifdef AIX
    if((lseek64(fh, startpos*sizeof(double), SEEK_SET)) < 0)
#else
#ifdef BLUEGENE
    if((lseek64(fh, startpos*sizeof(double), SEEK_SET)) < 0)
#else
    if((lseek(fh, startpos*sizeof(double), SEEK_SET)) < 0)
#endif
#endif
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Backup: Error moving file pointer, errorno=%d, fh %d startpos %ld, size %ld\n", me,errno, fh, startpos, size);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if((write( fh, (char*)buffer, size*sizeof(double))) < 0 )
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Backup: Error writing to file, errorno=%d\n", me, errno); 
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
}

f_int F77_NAME(f_backupram, F_BACKUPRAM)(f_int *fhandle, f_int *diskloc,
                     f_int *blksize, double *buffer, f_int *fsize)
{
   int fh = *fhandle;
   long size = *fsize;
   long long startpos;
   long long loc, bsize;
   f_int rc;

   rc = 0;
   loc = *diskloc-1;
   bsize = *blksize;
   startpos = loc * bsize;
   if (startpos < 0) 
   {
      printf("Error: Integer overflow in f_backupram: startpos %lld\n",startpos);
      printf("       diskloc, fsize = %d %d\n",*diskloc, *fsize);
      rc = -1;
      return rc;
   }

   backupram(fh, startpos, buffer, size);

   return rc;
}

void F77_NAME(f_write_disk, F_WRITE_DISK) (f_int *fh, long long *startpos, 
                                   double *buffer, f_int *fsize)
{
   int size = *fsize;
   int handle = *fh;
   long long offset = *startpos;

   backupram(handle, offset, buffer, size);
}

void F77_NAME(f_write_disk_ints, F_WRITE_DISK_INTS) (f_int *fh, long long *startpos, 
                                   f_int *buffer, f_int *fsize)
{
   int size = *fsize;
   int handle = *fh;
   long long offset = *startpos*sizeof(f_int);

    int me;

#ifdef AIX
    if((lseek64(handle, offset, SEEK_SET)) < 0)
#else
#ifdef BLUEGENE
    if((lseek64(handle, offset, SEEK_SET)) < 0)
#else
    if((lseek(handle, offset, SEEK_SET)) < 0)
#endif
#endif
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Backup: Error moving file pointer, errorno=%d, fh %d offset %ld, size %ld\n", me,errno, fh, offset, size);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if((write( handle, (char*)buffer, size*sizeof(double))) < 0 )
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Backup: Error writing to file, errorno=%d\n", me, errno); 
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

/*==========================================================================
restore file from disk to ram, 
make up by the file name:
/(tmpdirectory)/Asmgr<current>.tmp
<current> will be filled with 0 to the length of filenumber.
copy data of double type from buffer with the size of (int) size.
return -1 if can not open file
return 1 if successful.
==========================================================================*/
int restoreram(int fh, long long startpos, double* buffer, 
               long size)
{
    int me;

#ifdef AIX
    if((lseek64(fh, startpos*sizeof(double), SEEK_SET)) < 0)
#else
#ifdef BLUEGENE
    if((lseek64(fh, startpos*sizeof(double), SEEK_SET)) < 0)
#else
    if((lseek(fh, startpos*sizeof(double), SEEK_SET)) < 0)
#endif
#endif
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Restore: Error move file pointer, errorno=%d, startpos %ld\n", me,errno,startpos);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if((read( fh, (char*)buffer, size*sizeof(double))) <0 )
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Restore: Error reading from file, errorno=%d\n", me, errno); 
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
}

void F77_NAME(f_restoreram, F_RESTORERAM) (f_int *fhandle, f_int *diskloc, 
                  f_int *blksize, double *buffer, f_int *fsize)
{
   int fh = *fhandle;
   long size = *fsize;
   long long startpos;
   long long loc, bsize;

   loc = *diskloc-1;
   bsize = *blksize;
   startpos = loc * bsize;

   restoreram(fh, startpos, buffer, size); 
}

void F77_NAME(f_read_disk, F_READ_DISK) (f_int *fh, long long *startpos, 
                                   double *buffer, f_int *fsize)
{
   long size = *fsize;
   int handle = *fh;
   long long offset = *startpos;

   restoreram(handle, offset, buffer, size);
}

void F77_NAME(f_read_disk_ints, F_READ_DISK_INTS) (f_int *fh, long long *startpos, 
                                   f_int *buffer, f_int *fsize)
{
   /* Size is in words, not bytes. */

   long size = *fsize;
   int handle = *fh;
   long long offset = *startpos*sizeof(f_int);
   int me; 

#ifdef AIX
    if((lseek64(handle, offset, SEEK_SET)) < 0)
#else
#ifdef BLUEGENE
    if((lseek64(handle, offset, SEEK_SET)) < 0)
#else
    if((lseek(handle, offset, SEEK_SET)) < 0)
#endif
#endif
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d f_read_disk_ints: Error move file pointer, errorno=%d, offset %ld\n", me,errno,offset);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if((read( handle, (char*)buffer, size*sizeof(f_int))) <0 )
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
        printf("Server %d Restore: Error reading from file, errorno=%d\n", me, errno);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

}

void deletefile (char *filename, int fh)
{
    if (fh >= 0) close(fh);
    if( remove( (const char* )filename) == -1 )
    {
        printf("Error in deleting %s, errno value is %d\n", filename, errno);
    }
}

void F77_NAME(f_deletefile,F_DELETEFILE)(char *fn, f_int *fh)
{
    int handle = *fh;

     deletefile(fn, handle);
}

void F77_NAME(f_renamefile, F_RENAMEFILE)(char *f1, char *f2)
{
   rename(f1, f2);
}

void F77_NAME(get_filelen, GET_FILELEN) (f_int *fh, long long *filelength)
{
     
   int handle = *fh;
   off_t offset = 0;
   off_t len;
#ifdef AIX
    len = lseek64(handle, offset, SEEK_END);
#else
#ifdef BLUEGENE
    len = lseek64(handle, offset, SEEK_END);
#else
    len = lseek(handle, offset, SEEK_END);
#endif
#endif

   *filelength = len;
}

