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
instruct.cpp
Lei Wang
wang@qtp.ufl.edu
May 2004
******************************/
#include <iostream>
#include <string>
#include <cstdlib>
#include "qcarray.h"
#include "variable.h"
#include "keywdcnt.h"
#include "instruct.h"
//#include "testf.h"
//#include "parsetest.h"
#include "interfacec.h"
#include "errorhl.h"
#include "parser_interface.h"

using namespace std;

extern QCArrayClass QCArray;
extern VariableClass Variable;
extern string idname[];
extern int glno;
extern int mx_array_dim;

extern VariableClass InsVar;
/*
extern "C"{
extern void* fn1();
extern void* fn2();
}
*/


InstructCls::InstructCls()
{
    numberofins=1;
    stackno=1000;
	stacknof=1000;
    instruction = NULL;
}

void InstructCls::initializer()
{
    numberofins=1;
    stackno=1000;
	stacknof=1000;

    // Allocate memory for instructions

    instruction = (int **) malloc( maxinstructionnumber*sizeof(int *));
    for (int i = 0; i < maxinstructionnumber; i++)
       instruction[i] = (int *) malloc((6+mx_array_dim)*sizeof(int));
}



InstructCls::~InstructCls()
{
    numberofins=0;
    if (instruction) 
    {
       for (int i = 0; i < maxinstructionnumber; i++)
          if (instruction[i]) free(instruction[i]);
       free(instruction);
    }
}



int InstructCls::HaveEnoughSize()
{
    if (numberofins>=maxinstructionnumber)
    {
        outerror(glno, 
            "More istr id than instr can hold, please contact writer.");
        return 0;
    }
    return 1;

}



/*insert allocate, deallocate*/
void InstructCls::insertallocate(int type, int idno, int nindex, int*indexarray)
{
    int i;

    if (!HaveEnoughSize())  return;

    if (type!=allocatee&&type!=deallocatee && 
        type != create && type != deletee)
    {
        outerror(glno, 
            "Internal error in instruction insertallocate.");
        return;
    }

    if (CheckConvertAllocate(idno, nindex, indexarray)!=0)
        return; 
    if ((type == allocatee || type == deallocatee) && 
       QCArray.gettype(idno)!=local)
    {
        outerror(glno, "allocate has to use local type array.");
        return;
    }
    
    if ((type == create || type == deletee) &&
        QCArray.gettype(idno)!=distributed)
    {
       outerror(glno, "create has to use a distributed array.");
       return;
    }

    instruction[numberofins][0]=type;
    instruction[numberofins][1]=idno;
    for (i=0; i<nindex; i++)
        instruction[numberofins][i+2]=indexarray[i];
    /*add 0 to meet the ll needs*/
    for (i=nindex; i<mx_array_dim; i++)
        instruction[numberofins][i+2]=0;
    instruction[numberofins][2+mx_array_dim]=nindex;
    instruction[numberofins][4+mx_array_dim]=glno;    
    numberofins++;
}

/*setindex, get, computeint*/
void InstructCls::insert(int type, int idno, int nindex, int*indexarray)
{
    int i;

    if (!HaveEnoughSize())  return;

    if (type!=setindex&&type!=get&&type!=computeint)
    {
        outerror(glno, 
            "Internal error in instruction insert, please contact writer");
        return;
    }

    if (get==type)
    {
        if (CheckConvert(idno, nindex, indexarray)!=0)
            return;
        if (QCArray.gettype(idno)!=distributed)
        {
            outerror(glno, "get has to use distributed type array.");
            return;
        }
    }
    else if (computeint==type)
    {
        if (CheckConvert(idno, nindex, indexarray)!=0)
            return;
        /* if (QCArray.gettype(idno)!=served)
        {
            outerror(glno, "compute_integrals has to use served type array.");
            return;
        } */
    }

    instruction[numberofins][0]=type;
    instruction[numberofins][1]=idno;
    for (i=0; i<nindex; i++)
        instruction[numberofins][i+2]=indexarray[i];
    /*add 0 to meet the ll needs*/
    for (i=nindex; i<mx_array_dim; i++)
        instruction[numberofins][i+2]=-1;
    instruction[numberofins][4+mx_array_dim]=glno;    
    numberofins++;
}



/*put, prepare, assign, collective, prequest*/
void InstructCls::insert(int type, int idno, int nindex, int*indexarray,
                         int idno1, int nindex1, int* indexarray1)
{

    if (!HaveEnoughSize())  return;

    if (type!=put&&type!=putadd&&type!=prepare&&type!=prepareadd
        &&type!=assign&&type!=collective && type != prequest)
    {
        outerror(glno, 
            "Internal error in instruction insert, please contact writer");
        return;
    }

    if(CheckConvert(idno, nindex, indexarray))
        return;
    if (QCArray.istempscalar(idno))
    {
        outerror(glno, "constant can not be modified.");
        return;
    }

    if(CheckConvert(idno1, nindex1, indexarray1))
        return;

    if(collective==type)
    {
        if (QCArray.Checkcollective(nindex, indexarray, nindex1, indexarray1))
            return;
    }

    if(assign==type)
    {
        if (QCArray.CheckAssign(nindex, indexarray, nindex1, indexarray1))
            return;
    }

    instruction[numberofins][0]=type;
    instruction[numberofins][1]=idno;
    instruction[numberofins][2]=idno1;
/*
    if (put==type && QCArray.CheckAccIndex(arrayno, array1no,
        nindex, indexarray,
        nindex1, indexarray1))
        return;
    if (assign==type)
    {
        if (QCArray.CheckAssignIndex(arrayno, array1no,
        nindex, indexarray,
        nindex1, indexarray1))
            return;
    }
*/
    instruction[numberofins][4+mx_array_dim]=glno;    
    numberofins++;
}



int InstructCls::CheckConvertAllocate(int &idno, int nindex, int indexarray[])
{
    int result;
    result=QCArray.CheckArrayAllocate(idno, nindex, indexarray);
    if (result) return result;
    return 0;
}



int InstructCls::CheckConvert(int &idno, int nindex, int indexarray[])
{
    int result;
    result=QCArray.CheckArray(idno, nindex, indexarray);
    if (result) return result;

    if (!QCArray.isscalar(idno))
        insert(setindex, idno, nindex, indexarray);
    return 0;
}



/*contraction, summation, substraction, tensor*/
void InstructCls::insert(int type, int idno, int nindex1, int indexarray1[],
                         int id2no, int nindex2, int indexarray2[],
                         int id3no, int nindex3, int indexarray3[])
{
    int nindexr=0, indexarrayr[mx_array_dim];
    int i, j;

    if (!HaveEnoughSize())  return;

    if(type!=mult&&type!=add&&type!=sub&&type!=divide_array&&type!=tensor)
    {
        outerror(glno, "unknown error in instruct.insert #1.");
        return;
    }

    if (CheckConvert(idno, nindex1, indexarray1))
        return;

    if (QCArray.istempscalar(idno))
    {
        outerror(glno, "constant can not be modified.");
        return;
    }

    if (CheckConvert(id2no, nindex2, indexarray2))
        return;
    if (CheckConvert(id3no, nindex3, indexarray3))
        return;
    if (type==add||type==sub)
    {
        if (QCArray.CheckSumDiff(nindex1, indexarray1,
            nindex2, indexarray2,
            nindex3, indexarray3))
            return;
    }
    if (type==mult)
        if (QCArray.CheckMult(nindex1, indexarray1,
            nindex2, indexarray2, 
            nindex3, indexarray3,
            nindexr, indexarrayr))
            return;
    if (type==tensor)
        if (QCArray.CheckTensor(nindex1, indexarray1,
            nindex2, indexarray2, 
            nindex3, indexarray3,
            nindexr, indexarrayr))
            return;
    if (type==divide_array)
        if (QCArray.CheckDivide(nindex1, indexarray1,
            nindex2, indexarray2, 
            nindex3, indexarray3))
            return;

    if (mult==type)
        instruction[numberofins][0]=contract;
    else if (add==type)
        instruction[numberofins][0]=sum;
    else if (sub==type)
        instruction[numberofins][0]=diff;
    else if (divide_array==type)
        instruction[numberofins][0]=divide_array;
    else if (tensor==type)
        instruction[numberofins][0]=tensor;
    else
        outerror(glno, "unknown error in instruct.insert #2.");

    instruction[numberofins][1]=idno;
    instruction[numberofins][2]=id2no;
    instruction[numberofins][3]=id3no;

    for (i=0; i<nindexr; i++)
        instruction[numberofins][i+4]=indexarrayr[i];

    for (j=i; j<mx_array_dim; j++)
        instruction[numberofins][j+4]=-1;

    instruction[numberofins][4+mx_array_dim]=glno;    
    numberofins++;
}


void InstructCls::insert(int inscode, int nindex, int*indexarray, 
                         int flagarray[])
{
    /* Note: flagarray[j] = 0 if the jth index is referenced normally,
                          = 1 if the index is a "ii in i" reference. */

    int i, index;

    if (!HaveEnoughSize())  return;

    if(inscode==pardo||inscode==endpardo)
    {
        instruction[numberofins][0]=inscode;
        instruction[numberofins][1]=nindex;
        for (i=0; i<nindex; i++)
        {
            if (Variable.getidnumber(indexarray[i])<0)
                outerror(glno, "index not found.");
            else
            {
                index = Variable.getidnumber(indexarray[i]);
                // If this is a "in" reference, or on the in_index_mask

                if (flagarray[i]) index |= in_index_mask;
                instruction[numberofins][i+2] = index;
            }
        }
        /*add 0 to meet the ll needs*/
        for (i=nindex; i<mx_array_dim; i++)
            instruction[numberofins][i+2]=-1;

    instruction[numberofins][4+mx_array_dim]=glno;    
    numberofins++;
    }
    else
        outerror(glno, "unknown error in instruct.insert(int, int, int*).");
}


/* create, delete, destroy, proc, endproc, call, start, cycle, exit*/
void InstructCls::insert(int inscode, int idno)
{
    int variableno;
    int i;

    if (!HaveEnoughSize())  return;

    if (inscode==create||inscode==deletee || inscode == destroye)
    {
        variableno=QCArray.getidnumber(idno);
        if (variableno==-1)
        {
            outerror(glno, "create/delete/destroy: index can not be found.");
            printf("idno = %d, variableno %d\n",idno,variableno);
            printf("Name for variable is %s\n",QCArray.getidname(idno).c_str());
            return;
        }
        else
        {
            instruction[numberofins][0]=inscode;
            instruction[numberofins][1]=variableno;

            // Full create: Set all indices to 0.

            for (i = 0; i < mx_array_dim; i++)
               instruction[numberofins][i+2] = 0;

            instruction[numberofins][4+mx_array_dim]=glno;    
            numberofins++;
	}
    }
    else if (inscode==proc)
    {
        if (Variable.getidnumber(idno)!=-1||QCArray.getidnumber(idno)!=-1)
        {
            outerror(glno, "id exist.");
            return;
        }
        InsVar.insert(0, idno, numberofins, 0);
    }
    else if (inscode==endproc)
    {
        instruction[numberofins][0]=endproc;
        instruction[numberofins][4+mx_array_dim]=glno;
        numberofins++;
    }
    else if (inscode==call)
    {
        variableno=InsVar.getidnumber(idno);
        if (variableno==-1)
        {
            outerror(glno, "call: index can not be found.");
            return;
        }
        else
        {
            instruction[numberofins][0]=call;
            instruction[numberofins][1]=InsVar.getrange1(variableno);
            instruction[numberofins][4+mx_array_dim]=glno;
            numberofins++;
        }
    }
    else if (inscode==exitt||inscode==cycle)
    {
        instruction[numberofins][0]=inscode;
        numberofins++;
    }
    else if (inscode==start)
    {
        instruction[0][0]=jump;
        instruction[0][2]=numberofins;
    }
    
    else
        outerror(glno, "unknown error in instruct.insert #3.");
}



/*do, enddo, execute, create_index (with 2 args) */
void InstructCls::insert(int inscode, int isno, int idno)
{
    int i, variableno;
    if (!HaveEnoughSize())  return;

    if(inscode==doo||inscode==enddo)
    {
        variableno=Variable.getidnumber(isno);
        if (variableno==-1)
        {
            outerror(glno, "do/enddo: index can not be found.");
            return;
        }
        else
        {
            instruction[numberofins][0]=inscode;
            if (idno) variableno |= in_index_mask;
            instruction[numberofins][1]=variableno;
            instruction[numberofins][4+mx_array_dim]=glno;
            numberofins++;
        }
    }
    else if (inscode == execute) 
    {
       instruction[numberofins][3+mx_array_dim]=
           c_get_subroutine_handle((char*)idname[isno].c_str());

       if(instruction[numberofins][3+mx_array_dim]>0)
       {
           instruction[numberofins][0]=inscode;

           if (idno==-1)       /*no argument*/
               instruction[numberofins][1]=-1;
           else
           {
               instruction[numberofins][1]=QCArray.getidnumber(idno);

               if (instruction[numberofins][1]<0)
               {
                   outerror(glno,
                       "array name or scalar name not defined.");
                   return;
               }
           }

           for (i=0; i<mx_array_dim; i++)
               instruction[numberofins][i+2]=-1;
           instruction[numberofins][2+mx_array_dim]=0;
           instruction[numberofins][4+mx_array_dim]=glno;
           instruction[numberofins][5+mx_array_dim]=-1;
           numberofins++;
       }
       else
       {
           outerror(glno, "subroutine not in table");
           printf("c_get_subroutine_handle returned %d.\n", 
            instruction[numberofins][1]);
       }
    }
    else
    {
        outerror(glno, "unknown error in instruct.insert #4.");
        return;
    }
}



/*execute, where*/
void InstructCls::insert(int inscode, int isno, int idno, int idno1)
{
    int i;
    if (!HaveEnoughSize())  return;

    if (inscode!=execute &&
        inscode!=where)
    {
        outerror(glno, "unknown error in instruct.insert #5.");
        return;
    }

    if (inscode == execute) 
    {
       instruction[numberofins][3+mx_array_dim]=
	    c_get_subroutine_handle((char*)idname[isno].c_str());
       if(instruction[numberofins][3+mx_array_dim]>0)
       {
           instruction[numberofins][0]=inscode;
           instruction[numberofins][1]=QCArray.getidnumber(idno);
           if (instruction[numberofins][1]<0)
           {
               outerror(glno,
                   "array name or scalar name not defined.");
               return;
           }

           if (idno1 >=0)
              instruction[numberofins][5+mx_array_dim]=QCArray.getidnumber(idno1);
           else
              instruction[numberofins][5+mx_array_dim]=-idno1;
           if (instruction[numberofins][1]<0)
           {
               outerror(glno,
                   "array name or scalar name not defined.");
               return;
           }

           for (i=0; i<mx_array_dim; i++)
               instruction[numberofins][i+2]=-1;
           instruction[numberofins][2+mx_array_dim]=0;
           instruction[numberofins][4+mx_array_dim]=glno;
           numberofins++;
       }
       else
       {
        outerror(glno, "subroutine not in table");
        printf("c_get_subroutine_handle returned %d.\n", 
            instruction[numberofins][1]);
       }
    }   /* inscode == execute */

    if (inscode == where)
    {
        for (i=0; i<mx_array_dim; i++)
            instruction[numberofins][i+2]=-1;  /* set all indices to -1 */

        instruction[numberofins][0]=inscode;
	instruction[numberofins][1] = Variable.getidnumber(isno)+1;
	instruction[numberofins][2] = idno1;
        instruction[numberofins][3] = Variable.getidnumber(idno)+1;
        instruction[numberofins][4+mx_array_dim]=glno;
        numberofins++;
    }
}



/*execute, request*/
void InstructCls::insert(int inscode, int isno, 
int idno, int nindex, int indexarray[])
{
    int i;
    int temp1;

    if (!HaveEnoughSize())  return;

    if (inscode!=execute && inscode!=request)
    {
        outerror(glno, "unknown error in instruct.insert #6.");
        return;
    }

    if(CheckConvert(idno, nindex, indexarray))
        return;

    if (execute==inscode)
    {
        temp1=c_get_subroutine_handle((char*)idname[isno].c_str());

        /*the subroutine handle is always greater than 0*/
        if (temp1<=0)
        {
            outerror(glno, "subroutine not in table");
            printf("c_get_subroutine_handle returned %d.\n", temp1);
            return;
        }
        instruction[numberofins][5+mx_array_dim]=-1;
    }
    else /*request==inscode*/
    {
        if (QCArray.gettype(idno)!=served)
        {
            outerror(glno, "request has to use served type array.");
            return;
        }

        int find=0;
        temp1=Variable.getidnumber(isno);
        instruction[numberofins][3+mx_array_dim]=temp1;

        if (temp1<0)
        {
            outerror(glno, "variable not defined");
            return;
        }
        for (i=0; i<nindex; i++)
        {
            if (indexarray[i]==temp1)
            {
                find=1;
                break;
            }
        }
        if (find==0)
        {
            outerror(glno, "variable after the array not found in array before it.");
            return;
        }
    }
    instruction[numberofins][3+mx_array_dim]=temp1;
    instruction[numberofins][0]=inscode;
    instruction[numberofins][1]=idno;

    for (i=0; i<nindex; i++)
    {
        instruction[numberofins][i+2]
            =indexarray[i];
    }
    /*add 0 to meet the ll needs*/
    for (i=nindex; i<mx_array_dim; i++)
        instruction[numberofins][i+2]=-1;

    instruction[numberofins][2+mx_array_dim]=nindex;
    instruction[numberofins][4+mx_array_dim]=glno;
    numberofins++;
}



/*fconst*/
int InstructCls::insertsf(int inscode, double result, int op1, int op2)
{
    int rvalue=0;

    if (!HaveEnoughSize())  return 0;

    instruction[numberofins][0]=idf;

    if (fconst == inscode)
    {

        instruction[numberofins][1]=stacknof;
        instruction[numberofins][2]=QCArray.inserttempscalar(result);

        rvalue=stacknof++;
    }

    else
        outerror(glno, "unknown error in instruct.insert fconst.");

    instruction[numberofins][4+mx_array_dim]=glno;
    numberofins++;
    return rvalue;
}



/*iconst, id, idf, symbolic_const*/
int InstructCls::inserts(int inscode, int result, int op1, int op2)
{
    int rvalue=0;

    if (!HaveEnoughSize())  return 0;

    instruction[numberofins][0]=inscode;

    if (inscode == symbolic_const) 
    {
       instruction[numberofins][1]=stackno;
        instruction[numberofins][2]=result;
        rvalue=stackno++;
    }
    else if (iconst==inscode)
    {
        instruction[numberofins][1]=stackno;
        instruction[numberofins][2]=result;
        rvalue=stackno++;
    }
    else if (id== inscode)
    {
        instruction[numberofins][1]=stackno;
        instruction[numberofins][2]=result;
        rvalue=stackno++;
    }
    else if (idf== inscode)
    {
        instruction[numberofins][1]=stacknof;
        instruction[numberofins][2]=result;
        rvalue=stacknof++;
    }
    else if (add==inscode||sub==inscode||mult==inscode||divide==inscode
            ||eq==inscode||ne==inscode||ge==inscode
            ||le==inscode||gt==inscode||lt==inscode)
    {
        instruction[numberofins][1]=result;
        instruction[numberofins][2]=op1;
        instruction[numberofins][3]=op2;
        rvalue=result;
        stackno--;
    }
    else if (addf==inscode||subf==inscode||multf==inscode||dividef==inscode
            ||eqf==inscode||nef==inscode||gef==inscode
            ||lef==inscode||gtf==inscode||ltf==inscode)
    {
        instruction[numberofins][1]=result;
        instruction[numberofins][2]=op1;
        instruction[numberofins][3]=op2;
        rvalue=result;
        stacknof--;
    }
    else
        outerror(glno, "unknown error in instruct.insert iconst.");

    instruction[numberofins][4+mx_array_dim]=glno;
    numberofins++;
    return rvalue;
}



int InstructCls::insertif (int inscode, int condition)
{
    if (!HaveEnoughSize())  return 0;

    instruction[numberofins][0]=inscode;

    switch (inscode)
    {
    case iff:
        instruction[numberofins][1]=condition;
        stackno--;
        break;

    case elsee:
        instruction[numberofins][0]=jump;
        break;

    default:
        outerror(glno, "unknown error in instruct.insert if.");
    }
    instruction[numberofins][4+mx_array_dim]=glno;
    return numberofins++;
}


int InstructCls::reinsert(int structno, int address)
{
    if (!HaveEnoughSize()) return 0;

    instruction[structno][2]=address;
    return 1;
}


int InstructCls::getinsno() const
{
    return numberofins;
}



void InstructCls::output() const
{
    int i, j;
    for (i=0; i<numberofins; i++)
    {
        if (add==instruction[i][0]||sub==instruction[i][0]
            ||mult==instruction[i][0]||divide_array==instruction[i][0]
            ||eq==instruction[i][0]||ne==instruction[i][0]||ge==instruction[i][0]
            ||le==instruction[i][0]||gt==instruction[i][0]||lt==instruction[i][0]
            ||addf==instruction[i][0]||subf==instruction[i][0]
            ||multf==instruction[i][0]||dividef==instruction[i][0]            
            ||eqf==instruction[i][0]||nef==instruction[i][0]||gef==instruction[i][0]
            ||lef==instruction[i][0]||gtf==instruction[i][0]||ltf==instruction[i][0])
        {
        printf("%3d, %10s, %3d, %3d, %3d",i, keywords[instruction[i][0]], 
            instruction[i][1], instruction[i][2], instruction[i][3]);
        }
        else if (iconst==instruction[i][0]||iff==instruction[i][0]
            ||id==instruction[i][0]||idf==instruction[i][0] ||
              instruction[i][0]==symbolic_const)
        {
        printf("%3d, %10s, %3d, %3d",i, keywords[instruction[i][0]], 
            instruction[i][1], instruction[i][2]);
        }
        else if (jump==instruction[i][0])
        {
        printf("%3d, %10s, %3d",i, keywords[instruction[i][0]], 
            instruction[i][2]);
        }
        else if (call==instruction[i][0])
        {
        printf("%3d, %10s, %3d",i, keywords[instruction[i][0]], 
            instruction[i][1]);
        }
        else if (endproc==instruction[i][0])
        {
        printf("%3d, %10s, %3d",i, keywords[instruction[i][0]], 
            instruction[i][1]);
        }
        else if (pardo==instruction[i][0]||endpardo==instruction[i][0])
        {
            printf("%3d, %10s",i, keywords[instruction[i][0]]);
            for (j=0; j<instruction[i][1]; j++)
                printf(", %3d", instruction[i][j+2]);
        }
        else if (where == instruction[i][0])
        {
            printf("%3d, %10s %3d, %3d, %3d",i, keywords[instruction[i][0]],
                 instruction[i][1], instruction[i][2], instruction[i][3]);
        }
        else if (execute==instruction[i][0])
        {
            printf("%3d, %10s, %3d",i, keywords[instruction[i][0]], 
                instruction[i][1]);
            for (j=0; j<instruction[i][6]; j++)
                printf(", %3d", instruction[i][j+2]);
            printf(", %3d", instruction[i][7]);
        }
        else if (request==instruction[i][0])
        {
            printf("%3d, %10s, %3d",i, keywords[instruction[i][0]], 
                instruction[i][1]);
            for (j=0; j<instruction[i][6]; j++)
                printf(", %3d", instruction[i][j+2]);
            printf(", %3d", instruction[i][7]);
        }
        else if (allocatee==instruction[i][0])
        {
            printf("%3d, %10s, %3d",i, keywords[instruction[i][0]], 
                instruction[i][1]);
            for (j=0; j<instruction[i][6]; j++)
                printf(", %3d", instruction[i][j+2]);
        }
        else
        printf("%3d, %10s, ",i, keywords[instruction[i][0]]);
        /*for (int j=1; j<mx_array_dim; j++)
        {
            cout<<instruction[i][j]<<" ";
        }*/
        printf("\n");
    }
    //x();
}



void InstructCls::outputll() const
{
    int i;
    int t1=0, t2=0, t3=0, t5=0, t7=0;
    int type, lineno;
    int *ind_dum;

    ind_dum = new int[mx_array_dim];

    for (i = 0; i < mx_array_dim; i++)
       ind_dum[i] = 0;

    /*int isize=0;*/

    create_optable_c(numberofins);
    for(i=0; i<numberofins; i++)
    {
        lineno = instruction[i][4+mx_array_dim];

        if (contract==instruction[i][0])
        {
            type=101;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, &instruction[i][4], 
                t5, lineno);
        }
        else if (sum==instruction[i][0])
        {
            type=102;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, ind_dum, t5, lineno);
        }

        else if (setindex==instruction[i][0])
        {
            type=103;
            add_optable_c(type, t1, t2, instruction[i][1]+1, 
                &instruction[i][2], 
                t3, lineno);
        }
        else if (doo==instruction[i][0])
        {
            type=104;
            add_optable_c(type, t1, t2, t3, &instruction[i][1], 
                t7, lineno);
        }
        else if (enddo==instruction[i][0])
        {
            type=105;
            add_optable_c(type, t1, t2, t3, &instruction[i][1], 
                t7, lineno);
        }
        else if (get==instruction[i][0])
        {
            type=106;
            add_optable_c(type, t1, t2, instruction[i][1]+1, 
                &instruction[i][2], t3, lineno); 
        }
        else if (execute==instruction[i][0])
        {
            type=107;
            add_optable_c(type, instruction[i][5+mx_array_dim]+1, 
                t2, instruction[i][1]+1, 
                &instruction[i][2], 
                instruction[i][3+mx_array_dim], lineno);
        }
        else if (putadd==instruction[i][0])
        {
            type=108;
            add_optable_c(type, instruction[i][2]+1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (jump==instruction[i][0])
        {
            type=109;
            add_optable_c(type, t1, t2, instruction[i][2]+1, 
                ind_dum, t7, lineno);
        }
        else if (create==instruction[i][0])
        {
            type=110;
            add_optable_c(type, t1, t2, instruction[i][1]+1, 
                &instruction[i][2], t7, lineno);
        }
        else if (deletee==instruction[i][0])
        {
            type=111;
            add_optable_c(type, t1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (destroye == instruction[i][0])
        {
            type=162;
            add_optable_c(type, t1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (call==instruction[i][0])
        {
            type=112;
            add_optable_c(type, t1, t2, instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (endproc==instruction[i][0])                //return
        {
            type=113;
            add_optable_c(type, t1, t2, t3,  
                ind_dum, t7, lineno);
        }
        else if (iff==instruction[i][0])                    //jz
        {
            type=114;
            add_optable_c(type, instruction[i][1]+1, t2, 
                instruction[i][2]+1, 
                ind_dum, t7, lineno);
        }
        else if (add==instruction[i][0])
        {
            type=116;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (sub==instruction[i][0])
        {
            type=117;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (mult==instruction[i][0])
        {
            type=118;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (divide==instruction[i][0])
        {
            type=119;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (eq==instruction[i][0])
        {
            type=120;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (ne==instruction[i][0])
        {
            type=121;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (ge==instruction[i][0])
        {
            type=122;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (le==instruction[i][0])
        {
            type=123;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (gt==instruction[i][0])
        {
            type=124;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (lt==instruction[i][0])
        {
            type=125;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (iconst==instruction[i][0])
        {
            type=126;
            /*                  value,               t1,
                address,             t3, t4, t5, t6, t7*/
            add_optable_c(type, instruction[i][2], t1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (symbolic_const==instruction[i][0])
        {
            type=161;
            /*                  value,               t1,
                address,             t3, t4, t5, t6, t7*/
            add_optable_c(type, instruction[i][2], t1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (id==instruction[i][0])             //load value to
        {
            type=127;
            /*                  from address,        t1,
                to address,          t3, t4, t5, t6, t7*/
            add_optable_c(type, instruction[i][2]+1, t1, 
                instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (pardo==instruction[i][0])
        {
            type=128;
            /*                  from address,        t1,
                to address,          t3, t4, t5, t6, t7*/
            add_optable_c(type, t1, t2, t3, 
                &instruction[i][2], 
                t7, lineno);
        }
        else if (endpardo==instruction[i][0])
        {
            type=129;
            /*                  from address,        t1,
                to address,          t3, t4, t5, t6, t7*/
            add_optable_c(type, t1, t2, t3, 
                &instruction[i][2], 
                t7, lineno);
        }
        else if (exitt==instruction[i][0])
        {
            type=130;
            add_optable_c(type, t1, t2, t3, 
                ind_dum, t7, lineno);
        }
        else if (assign==instruction[i][0])
        {
            type=131;
            add_optable_c(type, instruction[i][2]+1, t1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (cycle==instruction[i][0])
        {
            type=132;
            add_optable_c(type, t1, t2, t3, 
                ind_dum, t7, lineno);
        }
        else if (diff==instruction[i][0])
        {
            type=135;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (collective==instruction[i][0])
        {
            type=136;
            add_optable_c(type, instruction[i][2]+1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (divide_array==instruction[i][0])
        {
            type=137;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (prepare==instruction[i][0])
        {
            type=138;
            add_optable_c(type, instruction[i][2]+1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (request==instruction[i][0])
        {
            type=139;
            add_optable_c(type, t1, t2, instruction[i][1]+1, 
                &instruction[i][2],
                t7,lineno);
        }
        else if (computeint==instruction[i][0])
        {
            type=140;
            add_optable_c(type, t1, t2, instruction[i][1]+1, 
                &instruction[i][2], 
                t3, lineno);
        }
        else if (put==instruction[i][0])
        {
            type=141;
            add_optable_c(type, instruction[i][2]+1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (tensor==instruction[i][0])
        {
            type=142;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, &instruction[i][4], 
                t5, lineno);
        }
        else if (addf==instruction[i][0])
        {
            type=146;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (subf==instruction[i][0])
        {
            type=147;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (multf==instruction[i][0])
        {
            type=148;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (dividef==instruction[i][0])
        {
            type=149;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (eqf==instruction[i][0])
        {
            type=150;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (nef==instruction[i][0])
        {
            type=151;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (gef==instruction[i][0])
        {
            type=152;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (lef==instruction[i][0])
        {
            type=153;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (gtf==instruction[i][0])
        {
            type=154;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (ltf==instruction[i][0])
        {
            type=155;
            add_optable_c(type, instruction[i][2]+1, instruction[i][3]+1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (idf==instruction[i][0])             //load value to
        {
            type=157;
            /*                  from address,        t1,
                to address,          t3, t4, t5, t6, t7*/
            add_optable_c(type, instruction[i][2]+1, t1, 
                instruction[i][1]+1, 
                ind_dum, t7, lineno);
        }
        else if (prepareadd==instruction[i][0])
        {
            type=158;
            add_optable_c(type, instruction[i][2]+1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (allocatee==instruction[i][0])
        {

			type=159;
            /*                  from address,        t1,
                to address,          t3, t4, t5, t6, t7*/
            add_optable_c(type, t1, t1, instruction[i][1]+1, 
			  &instruction[i][2], 
                          t7, lineno);
        }
        else if (deallocatee==instruction[i][0])
        {
            type=160;
            /*                  from address,        t1,
                to address,          t3, t4, t5, t6, t7*/
            add_optable_c(type, t1, t1, instruction[i][1]+1, 
                &instruction[i][2],
                t7, lineno);
        }
        else if (instruction[i][0] == prequest)
        {
           type = 163;
           add_optable_c(type, instruction[i][2]+1, t2, instruction[i][1]+1,
                ind_dum, t7, lineno);
        }
        else if (instruction[i][0] == where)
        {
           type = 164;
           add_optable_c(type, instruction[i][1], instruction[i][2],
                         instruction[i][3], ind_dum, t7, lineno);     
        }
        else
        {
            printf("Error in insert. Invalid code for instruction. ");
            printf("%d, %d\n", i, instruction[i][0]);
        }
    }

    delete [] ind_dum;
}
