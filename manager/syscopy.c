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
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>

#ifdef XT3
/* The Cray XT3 has only a micro-kernel, with no ability to execute a
   "system" call from within an executable program.  However, joda uses
   the system command to execute copy commands when doing its coares backup
   and coarse restore.  Therefore, this dummy system call intercepts the
   normal one and calls routines that do a manual file copy. 
*/ 

#include <dirent.h>

#define BUFSIZE 2048

void xt3_copy_file(char *input, char *output);
int creatfile(char *file, int *fh);
int openfile(char *file, int *fh);

int system(char *cmd)
{
   char *s1, *s2, *s3, *s4;
   char *line;
   char *original_dir; 
   int i, n, ierr;
   struct dirent **namelist;

   line = (char *) malloc(strlen(cmd)+1);
   strcpy(line, cmd);

   /* Parse the command line. */

   s1 = strtok(line, " ");
   s2 = strtok(NULL, " ");
   s3 = strtok(NULL, " ");
   s4 = strtok(NULL, " ");

   /* Do we have a valid command?  The only ones we can handle are "cp" and 
      "rm". */

   if (!strcmp(s1, "/bin/cp"))
   {
      /* Copy command. */

      if (!(s1 && s2 && s3 && s4) )   /* Must have 4 fields. */
      {
         printf("Unable to parse system command line: s1 = %s, s2 = %s, s3 = %s, s4 = %s4\n",s1,s2,s3,s4);
         MPI_Abort(MPI_COMM_WORLD, 1);
      }

      xt3_copy_file(s3, s4);
   }
   else if (!strcmp(s1, "/bin/rm"))
   {
      /* Remove command. */
 
      if (!(s1 && s2 && s3) )   /* Must have 3 fields. */
      {
         printf("Unable to parse system command line: s1 = %s, s2 = %s, s3 = %s\n",s1,s2,s3);
         MPI_Abort(MPI_COMM_WORLD, 1);
      }

      /* We must first remove each file in the directory before 
         removing the directory itself.  The scandir routine gives 
         us a structure with the directory entries.
      */

      original_dir = getcwd(NULL, 0);   /* Save the original working dir */

      n = scandir(s3, &namelist, 0,alphasort);
      if (n < 0)
      {
         perror("scandir");
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else
      {
         /* The scandir worked, now remove each file. */
 
         ierr = chdir(s3);  /* chdir to the new directory. */
         if (ierr < 0) 
         {
            printf("Error: Cannot chdir to %s\n",s3);
            perror("chdir");
            MPI_Abort(MPI_COMM_WORLD, 1);
         }

         for (i = 0; i < n; i++)
         {
            if (strcmp(namelist[i]->d_name, ".") && 
                strcmp(namelist[i]->d_name, "..")  )
            {
               ierr = remove(namelist[i]->d_name);
               if (ierr < 0)
               {
                  printf("Error: Unable to remove file %s from dir %s\n",namelist[i]->d_name, s3);
                  perror("remove");
               }
            }

            free(namelist[i]);
         }
      }
      free(namelist);   /* free up the directory structure */

      /* Now remove the directory. */

      rmdir(s3);

      ierr = chdir(original_dir);   /* Return to the original working dir */
      free(original_dir);           /* clean up the memory */ 
   }
   else
   {
      printf("XT3 SYSTEM CALL: Invalid command:\n");
      printf("command line is %s\n",cmd);
      MPI_Abort(MPI_COMM_WORLD, 1);
   } 

   free(line);
   return 0;
}

void xt3_copy_file(char *input, char *output)
{
   /* Copies file "input" to "output". */

   int ihandle, ohandle, ierr;
   char buf[BUFSIZE];
   size_t size, datasize;

   ierr = creatfile(output, &ohandle);
   if (ierr < 0) 
   {
       printf("Error creating file %s\n",output);
       MPI_Abort(MPI_COMM_WORLD, 1);
   }

   ierr = openfile(input, &ihandle);
   if (ierr < 0)
   {
       printf("Error opening file %s\n",input);
       MPI_Abort(MPI_COMM_WORLD, 1);
   }

   datasize = BUFSIZE;

   while (1)
   {
      size = read( ihandle, (char*)buf, datasize);
      if (size == 0) break;

      if (size < 0) 
      {
         printf("Error: read error while reading %s\n",input);
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
      write(ohandle, buf, size); 
   }
 
   close(ihandle);
   close(ohandle);
}

#endif


