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
#include <errno.h>

#define linesize 150

/**************************************************************************
 *
 *   top_mem_free: Gets the free memory by executing a "top" command and
 *                 decoding the "Memory free" field from the resulting
 *                 output.  If the top command fails, the function 
 *                 returns 0. 
 *
 *                 If the top command executes correctly, and the Memory 
 *                 free field is found, this function returns the free 
 *                 memory in megabytes.
 *
 **************************************************************************/ 
int top_mem_free()
{
    FILE* memfile;
    char temp[linesize];
    char *temp1;

    int maxlength;
    int freefound=0;
    int start=0;
    int end=0;
    int i, ii;
    int value=0;

    if (system("top >mem.txt") != 0) 
    {
       perror("Error executing top command");
       return 0;
    }

    memfile=fopen ("mem.txt", "r");
    
    do
    {
        temp1=fgets(temp, linesize, memfile);
        
        if (temp1)
        {
            i=0;
            while (temp[i]!='\0')
                i++;
                
            maxlength=i;
            if ((strncmp((const char*)temp, "Memory", 6))==0)
            {
                /*find keyword "free"*/
                i=6;
                while (i<maxlength-5&&!freefound)
                {
                    temp1=&(temp[i]);
                    if (strncmp((const char*)temp1, "free", 4)==0)
                        freefound=1;
                    else
                        i++;
                }
                
                ii=0;
                /*find first keword 0-9*/
                while (i-ii>0&&temp[i-ii]>'9'||temp[i-ii]<'0')
                    ii++;
                end=i-ii;
                
                while (i-ii>0&&temp[i-ii]<='9'&&temp[i-ii]>='0')
                    ii++;
                start=i-ii+1;
                
                value=0;
                for (i=start; i<=end; i++)
                {
                    value=temp[i]-48+value*10;
                }
            }
        }
    }
    while(temp1!=NULL&&!freefound);
    
    printf("value=%d\n", value);
    remove("mem.txt");
    return value;
}
