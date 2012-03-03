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

/*
 * This routine saves the precious files from geometry optimizations and
 * finite difference calculations to szSaveDir/CURRENT.
 */

#include <stdio.h>	/* for *printf */
#include <stdlib.h>	/* for system */
#include <unistd.h>	/* for access */
#include <sys/types.h>	/* for mkdir */
#include <sys/stat.h>	/* for mkdir */
#include <string.h>	/* for strcmp */
#include <strings.h>	/* for strncpy */

#include "f77_name.h"	/* for F77_NAME */
#include "f_types.h"	/* for f_int */

void
F77_NAME(coarse_backup,COARSE_BACKUP)
(const char * szSaveDirIn, f_int * iErr)
{

#define BUFLEN 80
#define FILLEN 80
#define CMDLEN 256
#define CURLEN 128
#define OLDLEN 128
    char szSaveDir[80], szCmd[CMDLEN], szCurrent[CURLEN], szOld[OLDLEN];
    if (!strncpy(szSaveDir,szSaveDirIn,80))
    {
        printf("\nERROR: unable to copy directory string\n");
        *iErr=1;
        return;
    }
    snprintf(&szCurrent[0],CURLEN,"%s/CURRENT",szSaveDir);
    snprintf(&szOld[0],    OLDLEN,"%s/OLD",    szSaveDir);

 /* CONDITION szSaveDir */
    if (access(szSaveDir,F_OK))
    {
        printf("Attempting to create %s for coarse grain restart . . . ",
               szSaveDir);
        if (mkdir(&szSaveDir[0],S_IRWXU))
        {
            printf("\nERROR: unable to create backup directory\n");
            *iErr=1;
            return;
        }
        printf("done\n");
    }
    else /* szSaveDir exists */
    {
        if (access(szSaveDir,R_OK|W_OK|X_OK))
        {
            printf("ERROR: %s does not have the proper permissions\n",
                   szSaveDir);
            *iErr=1;
            return;
        }
    }

#include "coarse_precious.h" /* for coarse_precious[] */
 /* MOVE CURRENT FILES TO OLD DIRECTORY */
    if (!access(&szCurrent[0],F_OK))
    {
        printf("Renaming %s to %s . . . ",szCurrent,szOld);
    	if (!access(&szOld[0],F_OK))
    	{
    		/* OLD directory exists, remove files in it */
            int i=0;
        	char szFileOld[FILLEN];
            while (*coarse_precious[i])
            {
            	snprintf(&szFileOld[0],FILLEN,"%s/%s",szOld,coarse_precious[i]);
            	FILE * oldfile = fopen(szFileOld,"rb");
            	if (oldfile != NULL)
            	{
            		fclose(oldfile);
            		if (remove(szFileOld) != 0)
            		{
                        printf("\nERROR: unable to remove %s/%s\n",
                               szOld,coarse_precious[i]);
                        *iErr=1;
                        return;
            		}
            	}
                i++;
            }
    	}
    	else if (mkdir(&szOld[0],S_IRWXU))
    	{
    		printf("ERROR: unable to create %s\n",szOld);
    		*iErr=1;
    		return;
    	}
        int i=0;
    	char szFileOld[FILLEN],szFileCur[FILLEN];
        while (*coarse_precious[i])
        {
        	snprintf(&szFileOld[0],FILLEN,"%s/%s",szOld,coarse_precious[i]);
        	snprintf(&szFileCur[0],FILLEN,"%s/%s",szCurrent,coarse_precious[i]);
        	FILE * oldfile = fopen(szFileCur,"rb");
        	if (oldfile != NULL)
        	{
        		fclose(oldfile);
        		if (rename(szFileCur,szFileOld) != 0)
        		{
        			printf("\nERROR: unable to move %s/%s to %s/%s\n",
        					szCurrent,coarse_precious[i],szOld,coarse_precious[i]);
        			*iErr=1;
        			return;
        		}
        	}
            i++;
        }
    }
    else if (mkdir(&szCurrent[0],S_IRWXU))
    {
        printf("ERROR: unable to create %s\n",szCurrent);
        *iErr=1;
        return;
    }

 /* SAVE PRECIOUS FILES */
   int i=0;
   char szBuf[BUFLEN], szFile[FILLEN], szFileOut[FILLEN];

	printf("Files copied to %s:",szCurrent);
    while (*coarse_precious[i])
    {
    	strncpy(&szFile[0],coarse_precious[i],FILLEN);
    	if (access(&szFile[0],F_OK))
    		printf(" !%s",coarse_precious[i]);
    	else
    	{
    		snprintf(&szFileOut[0],FILLEN,"%s/%s",szCurrent,coarse_precious[i]);
    		FILE * infile = fopen(szFile,"rb");
    		FILE * outfile = fopen(szFileOut,"wb");
    		if (infile == NULL || outfile == NULL)
    		{
    			printf("\nERROR: unable to copy %s to %s/%s\n",
    					szFile,szCurrent,coarse_precious[i]);
    			*iErr=1;
    			return;
    		}
    		int c;
    		while((c=getc(infile))!=EOF)
    			putc(c,outfile);
    		fclose(infile);
    		fclose(outfile);
    		printf(" %s",coarse_precious[i]);
    	}
    	i++;
    }
    printf("\n");

 /* REMOVE szSaveDir/OLD */
    if (!access(&szOld[0],F_OK))
    {
        printf("Removing old back up directory . . . ");
		/* OLD directory exists, remove files in it */
        int i=0;
    	char szFileOld[FILLEN];
        while (*coarse_precious[i])
        {
        	snprintf(&szFileOld[0],FILLEN,"%s/%s",szOld,coarse_precious[i]);
        	FILE * oldfile = fopen(szFileOld,"rb");
        	if (oldfile != NULL)
        	{
        		fclose(oldfile);
        		if (remove(szFileOld) != 0)
        		{
                    printf("\nERROR: unable to remove %s/%s\n",
                           szOld,coarse_precious[i]);
                    *iErr=1;
                    return;
        		}
        	}
            i++;
        }
        printf("done\n");
    }

    *iErr=0;
    return;
}


