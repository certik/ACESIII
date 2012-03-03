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
/******************************
lexical.cpp
Lei Wang
wang@qtp.ufl.edu
May 2004
******************************/
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iomanip>
//#include "lexical.h"
#include "keywdcnt.h"
#include "variable.h"

using namespace std;

int outArray[max_output_size]={0};
//the array that contains the index of output.

int outarray_size=0;
//real size of outArray

int errorArray[max_output_size]={0};
//the array that contains the index of output when meet errors.

int errorarray_size=0;
//real size of errorArray

#define max(a,b) ((a)>(b))?(a):(b)
/*
struct symbol
{
    string name;
    int type;
    int location;
};
*/

double floatconst[max_output_size];
int floatno=0;// for input
//int floatno1=0;for getout
int intconst[max_output_size];
int intno=0;


/*====================================
The function insertsymbol insert variable id into symboltable.

lexaryp:  the position of lexarray
  idno:     total number of ids
  symbolno: total number of symbols
  idname:    the name of the id need to be inserted into symbol table
  stackno:  the total number of stackitems

int insertsymbol(string id,int lexaryp, 
                  int &symbolno,  int &stackno, int type)
{
    int i;
    int identical =0;
    for (i=0; i<=symbolno; i++)
    {
        if(symboltable[i].name==id)
        {
            identical =1;
            break;
        }
    }
    if (!identical)
    {
        stackno++;
        symbolno++;

        
        symboltable[symbolno].name=id;
        //cout <<"new idname="<<id<<endl;
        symboltable[symbolno].type = type;
        symboltable[symbolno].location = stackno;
        return stackno;
    }
    else //the array exists in symboltable
    {
        cout<<"redeclaration====="<< i<<endl;
        return i;
    }
}
====================================*/





/************************************
The above is the syntax part.
Below is the lexical part
************************************/



/*====================================
insertp function inserts i as index of keywords
into the next available position of outArray 
====================================*/
void insertp(int i)
{
    outArray[outarray_size]=i;
    if (outarray_size<max_output_size)
        outarray_size ++;
    else
    {
        printf("outSourceArray length is out of bounds. Please contact code writer.\n");
        printf("   outarray_size = %d, max_output_size = %d\n",outarray_size,max_output_size);
        exit(1);
    }
}



/*====================================
inserterror function inserts i as the position of
the unknown character into the next 
available position of errorArray 
====================================*/
void inserterror(int i)
{
    errorArray[errorarray_size]=i;
    errorarray_size++;
}



/*====================================
insertid function inserts an id (character string) from position start 
to position end,and converts to string, then stores into the next 
available position string name[]. idno idicates the position
of the inserted string in string idname[].  
====================================*/
void insertid(int &idno, int start, int end, 
              string idname[], char sourceArray[])
{
    int i=0;
    int ii=0;
    
    char temp[128]="";
    for (i=start; i<=end;i++)
    {
        temp[ii++]=sourceArray[i];
    }
    //temp[ii]='\n';
    idname[idno]=temp;
    idno++;
}



/*====================================
insertfloat function inserts a float from char from position start to 
position end,and converts to float, then stores into the next available 
position in floatconst[]. floatno idicates the position of the inserted 
float in floatconst[].  
====================================*/
void insertfloat(char sourceArray[], int start, int end)
{
    char *temp = new char[end-start+2];
    int i;
    int ii=0;
    for (i=start; i<=end; i++)
        temp[ii++]=sourceArray[i];
    temp[ii]='\0';

    floatconst[floatno]=atof((const char *)temp);
    floatno++;
    delete temp;
}



/*====================================
insertint function inserts an int from char from position start to position
end, and converts to int, then stores into the next available position in 
intconst[]. intno idicates the position of the inserted int in intconst[].  
====================================*/
void insertint(char sourceArray[], int start, int end)
{
    char *temp = new char[end-start+2];
    int i;
    int ii=0;
    for (i=start; i<=end; i++)
        temp[ii++]=sourceArray[i];
    temp[ii]='\0';

    intconst[intno]=atoi((const char *)temp);
    intno++;
    delete temp;
}



/*====================================
lexical function lexicalizes the characters stored in sourceArray,
variable sAsize tells the function the size of sourceArray is. 
idno is the number of total Ids (including duplications).
====================================*/
void lexical(char sourceArray[], int sAsize, 
             int &idno, string idname[], int &rnumber)
{
    int psource =0;               //position of sourceArray
    int temp=0;
    int i;                        //for iteration
    int foundeq;
    int strsize=0;

    /*find start of sial*/
    while (psource+3<sAsize&&
        strncmp((const char*)&(sourceArray[psource]),"sial", 4))
    {
        if (sourceArray[psource]=='\n')
            rnumber++;
        psource++;
    }

    while (psource <sAsize)
    {
        if (sourceArray[psource]>='a'&& sourceArray[psource]<='z')
        //possible identifier
        {
            temp=psource+1;
            while (temp<sAsize&&
                    ((sourceArray[temp]>='a'&& sourceArray[temp]<='z')
                    ||(sourceArray[temp]<='9'&& sourceArray[temp]>='0')
                    ||sourceArray[temp]=='_')
                  )
            {
                temp++;//check the rest are 
            }

            //see if the possible identifier is a key word
            foundeq=0;
            for (i=0; i<noperat_keywords; i++)
            {

                strsize=max(temp-psource, signed(strlen(keywords[i])));
                if(!strncmp((const char*)&(sourceArray[psource]), 
                        keywords[i], 
                        strsize))
                {
                    insertp(i);
                    foundeq=1;
                    break;
                }
            }
            if (foundeq==0)
            {
                insertp(id);
                if (psource+maxidlen>=temp)
                    insertid(idno, psource, temp-1, idname, sourceArray);
                else
                    insertid(idno, psource, psource+maxidlen, idname, sourceArray);
            }
            psource=temp; /*move pointer to the next position*/
        }
        else if (sourceArray[psource] == '=')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')   //==
            {
                insertp(eq);
                psource+=2;
            }
            else
            {
                insertp(assign);
                psource++;
            }
        }
        else if (sourceArray[psource] == '!')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')   //!=
            {
                insertp(ne);
                psource+=2;
            }
            else
            {
                insertp(error);
                psource++;
            }
        }
        else if (sourceArray[psource] == '#')                   /*comment*/
        {
            temp=psource+1;
            while (temp<sAsize&&sourceArray[temp]!='\n')
                temp++;
            psource=temp;
        }
        else if (sourceArray[psource] == '+')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')    //+=
            {
                insertp(addassign);
                psource+=2;
            }
            else                                  //+
            {
                insertp(add);
                psource++;
            }
        }
        else if (sourceArray[psource] == '-')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')    //-=
            {
                insertp(subassign);
                psource+=2;
            }
            else                                  //-
            {
                insertp(sub);
                psource++;
            }
        }
        else if (sourceArray[psource] == '*')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')    //*=
            {
                insertp(multassign);
                psource+=2;
            }
            else                                  //*
            {
                insertp(mult);
                psource++;
            }
        }
        else if (sourceArray[psource] == '/')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')    // /=
            {
                insertp(divassign);
                psource+=2;
            }
            else                                  // /
            {
                insertp(divide);
                psource++;
            }
        }
        else if (sourceArray[psource] == '^')
        {
            insertp(tensor);
            psource++;
        }
        else if (sourceArray[psource] == '<')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')    // <=
            {
                insertp(le);
                psource+=2;
            }
            else                                  // <
            {
                insertp(lt);
                psource++;
            }
        }
        else if (sourceArray[psource] == '>')
        {
            if (psource+1<sAsize && sourceArray[psource+1] == '=')    // >=
            {
                insertp(ge);
                psource+=2;
            }
            else                                  // >
            {
                insertp(gt);
                psource++;
            }
        }
        else if (sourceArray[psource] <= char('9')
            &&sourceArray[psource]>= char ('0'))    //number
        {
            temp=psource+1;
            while (temp<sAsize && sourceArray[temp] <='9'
                &&sourceArray[temp] >='0')
            {
                temp++;
            }
            if (sourceArray[temp] == '.')               // float
            {
                temp=temp+1;
                while (temp<sAsize && sourceArray[temp] <='9'
                    &&sourceArray[temp] >='0')    
                {
                    temp++;
                }
                insertp(fconst);
                insertfloat(sourceArray, psource, temp-1);
                psource = temp;
            }
            else                                           //integer
            {
                insertp(iconst);
                insertint(sourceArray, psource, temp-1);
                psource = temp;//something
            }
        }//end of else if, end of number
        else if (sourceArray[psource] == ',')    // ,
        {
            insertp(comma);
            psource ++;
        }//end of else if, end of ','
        else if (sourceArray[psource] == '\n')   //\n
        {
            insertp(returnn);
            psource ++;
        }
        else if (sourceArray[psource] == ' '
            ||sourceArray[psource] == '\t')
        {
            psource ++;
        }
        else if (sourceArray[psource] == '(')
        {
            insertp(lparen);
            psource ++;
        }
        else if (sourceArray[psource] == ')') 
        {
            insertp(rparen);
            psource ++;
        }
        else if (sourceArray[psource] == '[')
        {
            insertp(lbracket);
            psource ++;
        }
        else if (sourceArray[psource] == ']')
        {
            insertp(rbracket);
            psource ++;
        }
        else                                    //any other case
        {
            insertp(error);
            inserterror(psource);
            psource ++;
        }

    }//end of while, end of sourceArray
}



/*====================================
output function printout lexicalized information to the output device,
the information are index stored in outArray.
index are changed meanings by string temp.
Once index of error meets, original character will be print out
using sourceArray.
====================================*/
void output(ostream &output, string temp[], char sourceArray[])
{
    int i;
    int errorcount=0;         //keep tracking how many errors have printout.

    output<<endl;
    for (i=0; i<outarray_size; i++)
        if (outArray[i]!=error)
            output<<temp[outArray[i]]<<endl;
        else                 //output is unrecongnized symbol
        {
            output<<"**unrecognized symbol: \""
                <<sourceArray[errorArray[errorcount]]<<"\""<<endl;
            errorcount++;
        }
    output<<endl;
}



/*====================================
output function printout one token to the output device,
the information are index stored in outArray. index are changed 
meanings by string temp. Once index of error meets, original 
character will be print out using sourceArray.
====================================*/
void outputone(ostream &output, string temp[], char sourceArray[], 
               int i, int &errorcount)
{
   
     //keep tracking how many errors have printout.

        if (outArray[i]!=error)
            output<<temp[outArray[i]]<<endl;
        else                 //output is unrecongnized symbol
        {
            output<<"**unrecongnized symbol: \""
                <<sourceArray[errorArray[errorcount]]<<"\""<<endl;
            errorcount++;
        }
        
}
