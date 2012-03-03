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
 * This routine restores the precious files for geometry optimizations and
 * finite difference calculations from szSaveDir/(OLD|CURRENT).
 *
 * It is called EVERY TIME joda runs since the first joda does not know
 * whether the program is a restart and the files are simply missing.
 */

#include <stdio.h>	/* for *printf */
#include <unistd.h>	/* for access */
#include <stdlib.h>	/* for system */
#include <string.h>	/* for strcmp */
#include <strings.h>	/* for strncpy */

#include "f77_name.h"   /* for F77_NAME */
#include "f_types.h"    /* for f_int */

void
F77_NAME(coarse_restore,COARSE_RESTORE)
(const char * szSaveDirIn, f_int * iErr)
{

#define CMDLEN 256
#define BAKLEN 128
    char szSaveDir[80], szCmd[CMDLEN], szBak[BAKLEN];
    if (!strncpy(szSaveDir,szSaveDirIn,80))
    {
        printf("\nERROR: unable to copy directory string\n");
        *iErr=1;
        return;
    }

 /* TEST IF szSaveDir/OLD EXISTS */
    snprintf(&szBak[0],BAKLEN,"%s/OLD",szSaveDir);
    if (access(szBak,F_OK))
    {
     /* SWITCH TO szSaveDir/CURRENT */
        snprintf(&szBak[0],BAKLEN,"%s/CURRENT",szSaveDir);
        if (access(szBak,F_OK))
        {
         /* no directories in szSaveDir -> firstrun */
            *iErr=0;
            return;
        }
    }

 /* TEST IF szBak IS USEABLE */
    if (access(szBak,R_OK|X_OK))
    {
        printf("ERROR: %s does not have the proper permissions\n",szBak);
        *iErr=1;
        return;
    }

 /* RESTORE PRECIOUS FILES */
    {
#include "coarse_precious.h" /* for coarse_precious[] */
#define BUFLEN 80
#define FILLEN 80
        int i=0; char szBuf[BUFLEN], szFile[FILLEN], szFileBak[FILLEN];

        printf("Files restored from %s:",szBak);
        while (*coarse_precious[i])
        {
            strncpy(&szFile[0],coarse_precious[i],FILLEN);
            snprintf(&szCmd[0],CMDLEN,"%s/%s",szBak,coarse_precious[i]);
            if (access(&szCmd[0],F_OK))
                printf(" !%s",coarse_precious[i]);
            else
            {
            	snprintf(&szFileBak[0],FILLEN,"%s/%s",szBak,coarse_precious[i]);
            	FILE * infile = fopen(szFileBak,"rb");
            	FILE * outfile = fopen(szFile,"wb");
            	if (infile == NULL || outfile == NULL)
            	{
                    printf("\nERROR: unable to restore %s/%s to %s\n",
                           szBak,coarse_precious[i],szFile);
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
    }

    *iErr=0;
    return;
}

