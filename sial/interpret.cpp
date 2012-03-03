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
/****************************************************
  Lei Wang
  wang@qtp.ufl.edu
  May 2004

  Program which does lexical analysis and syntax analysis of a 
  super instruction assembly language(SIAL).
****************************************************/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "lexical.h"
#include "keywdcnt.h"
#include "variable.h"
#include "syntax.h"
#include "qcarray.h"
#include "variable.h"
#include "instruct.h"
#include "interpret.h"
#include "parser_interface.h"
using namespace std;

/*#define debuglex*/
extern int outArray[];
extern int outarray_size;
extern int errors;
QCArrayClass QCArray;
VariableClass Variable;
InstructCls Instruct;
VariableClass InsVar;
extern int intconst[];
extern double floatconst[];
int idno=0;
string idname [max_id_size];


void allocidname()
{
    return ;
    //idname=new string [max_id_size];
    // if (!idname)
    // {
    //     printf("error in alloc idname\n");
    //     exit(1);
   //  }
}


/*====================================
readsourcefile reads the infile as characters, and stores them
into sourceArray, and returns the total number of characters
as integer. return -1 if there is an error.
====================================*/
int readsourcefile(char* infilename, char sourceArray[])
{
    int i=0;
    char temp;
    ifstream infile;

    infile.open(infilename);
    if(!infile)
    {
        cerr <<"Error. can not open source file:" <<infilename<<endl<<endl;
        return -1;
    }

    while(infile.get(temp)&&i<max_input_size)
    {
        /*convert to lower case*/
        if (temp<='Z'&&temp>='A')
            temp=temp+'a'-'A';

        sourceArray[i++]=temp;
    }
    infile.close();
    infile.clear();     /*have to use this in NT, but don't have to in Unix*/
    return i;
}



/*====================================
main function: will call some other functions to realize the requirement of 
this project. Because for this project the output is to the screen, some 
of the lines are set for comments for further usage of output to file.
====================================*/
int interpret (char* infilename, int lext, int indext, int arrayt, int inst)
{
/*
    char infilename[256];
    string outfilename ="out.txt";
    ofstream outfile (outfilename.c_str());
*/
    char *sourceArray;

    sourceArray=new char[max_input_size];
    if (!sourceArray)
    {
        printf("error in allocate sourceA. Please contact code writer with your input.");
        exit(1);
    }
    int sizeofsource = 0;

    int rnumber=1;

    // Set up global objects  

    QCArray.initializer();
    Variable.initializer();
    Instruct.initializer();
    InsVar.initializer();

    sizeofsource = readsourcefile(infilename, sourceArray);

    if (sizeofsource <0)
        return sizeofsource;

    lexical(sourceArray, sizeofsource, idno, idname, rnumber);


    int pos=0;
    Variable.insertconstant();
    functionsyntax(pos, rnumber);
/*
    if (errors==0)
    {
        printf("=           No error is printed by now,                   =\n");
        printf("=           it is safe to go to the next step.            =\n");
    }
*/
    if (errors)
    {
    printf("===========================================================\n");
        printf("=           errors found.                                 =\n");
        printf("=           please correct it before go to the next step. =\n");
    printf("===========================================================\n");
    }
    else
    {
    printf("   No compilation errors\n");
    }

    if (lext)
    {
        printf("\nThe following is printed for testing.\n\n");
        printf("\nlexical analysis\n");
        printf("================\n");
        printf("idno=%d\n", idno);

        int iconstprint=0;
        int fconstprint=0;
        for (int i=0, ii=0; i<outarray_size; i++)
        {
            if (outArray[i]==id)
            {
                printf("id{%s} ", idname[ii].c_str());
                ii++;
            }
            else if (outArray[i]==iconst)
            {
                printf("cosnt (%d)", intconst[iconstprint++]);
            }
            else if (outArray[i]==fconst)
            {
                printf("cosnt (%f)", floatconst[fconstprint++]);
            }

            else
            {
                printf("%s ", keywords[outArray[i]]);
            }
        }
    }

    if (indext)
    {
        printf("\nindex table\n");
        printf("================\n");
        Variable.print();
    }
    if (arrayt)
    {
        printf("\narray table\n");
        printf("================\n");
        QCArray.print();
    }
    if (inst)
    {
        printf("\ninstruction table\n");
        printf("================\n");
        Instruct.output();
    }
    
    if (errors==0)
    {
        Variable.outputll();
        QCArray.outputll();
        Instruct.outputll();
        return 0;
    }

    delete []sourceArray;
    return 1;
}
