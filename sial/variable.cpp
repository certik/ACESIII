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
variable.cpp
Lei Wang
wang@qtp.ufl.edu
May 2004
******************************/
#include <cstdio>
#include <string>
//using namespace std;
#include "variable.h"
#include "keywdcnt.h"
//#include "parsetest.h"
#include "parser_interface.h"
#include "qcarray.h"
#include "errorhl.h"
#include "f77_name.h"
#include "f_types.h"

extern string idname[];
extern int idno;
extern int glno;
extern QCArrayClass QCArray;


VariableClass::VariableClass()
{
    nvars=0;
}


void VariableClass::initializer()
{
    vararray = new string[maxidnumber];
    varrange1 = new int[maxidnumber];
    varrange2 = new int[maxidnumber];
    vartype   = new int[maxidnumber];
}

VariableClass::~VariableClass()
{
    nvars=0;

    if (vararray) delete [] vararray;
    if (varrange1) delete [] varrange1;
    if (varrange2) delete [] varrange2;
    if (vartype) delete [] vartype;

}


void VariableClass::insertconstant()
{
vararray[0]="moindex";
varrange1[0]=1;
varrange2[0]=-norb;
vartype[0]=moindex;

vararray[1]="moaindex";
varrange1[1]=1;
varrange2[1]=-norb;
vartype[1]=moaindex;

vararray[2]="mobindex";
varrange1[2]=1;
varrange2[2]=-norb;
vartype[2]=mobindex;

vararray[3]="aoindex";
varrange1[3]=1;
varrange2[3]=-norb;
vartype[3]=aoindex;

nvars=4;
}

    //skip the first index by now, may need to modify later.

void VariableClass::insert(int type, int idno, int range1, int range2)
{
    int i;

    if (QCArray.getidnumber(idno)!=-1)
    {
        outerror(glno, "id exist.");
        return;
    }

    i=getidnumber(idno);
    if (i!=-1)
    {
        outerror(glno, "id exist.");
        return;
    }
    else        //try insert new
    {
        if (nvars<maxidnumber)
        {
            vararray[nvars]=idname[idno];
            varrange1[nvars]=range1;
            varrange2[nvars]=range2;
            vartype[nvars]=type;
            nvars++;
        }
        else
        {
            outerror (glno, "More id than id container can hold.");
        }
    }
}


/*
int VariableClass::getidnumber(string idname)
{
    int i;
    for (i=0; i<nvars&&idname!=vararray[i]; i++)
        i++;
    if (i<nvars&&idname==vararray[i])
        return i;
    else
        return -1;
}
*/

int VariableClass::getidnumber(int idnumber) const
{
    int i;
    for (i=0; i<nvars&&idname[idnumber]!=vararray[i]; i++);

    if (i<nvars&&idname[idnumber]==vararray[i])
        return i;
    else
        return -1;
}



int VariableClass::getrange1(int idnumber) const
{
    if (idnumber <nvars)
        return varrange1[idnumber];
    else
        return -1;
}



int VariableClass::getidtype(int idnumber) const
{
    if (idnumber <nvars)
        return vartype[idnumber];
    else
        return -1;
}



int VariableClass::getbegrange(int idnumber) const
{
    if (idnumber <nvars)
    {
        return varrange1[idnumber];
    }
    else
        return -1;
}



int VariableClass::getendrange(int idnumber) const
{
    if (idnumber <nvars)
    {
        return varrange2[idnumber];
    }
    else
        return -1;
}



string VariableClass::getidname(int idnumber) const
{
    if (idnumber<nvars)
        return vararray[idnumber];
    else
        return 0;
}



void VariableClass::print() const
{
    int i;
    for (i=0; i<nvars; i++)
    {
        printf("%i, %9s, ", i, vararray[i].c_str());

        if (varrange1[i]<=0)
            printf("%5s,", keywords[-varrange1[i]]);
        else
            printf("%i, ", varrange1[i]);

        if (varrange2[i]<=0)
            printf("%5s\n", keywords[-varrange2[i]]);
        else
            printf("%i\n", varrange2[i]);
    }
}


void VariableClass::outputll() const
{
    f_int type;
    int i;
    f_int nvarstemp=nvars;
    f_int range1temp, range2temp;

    F77_NAME(create_index_table, CREATE_INDEX_TABLE) (&nvarstemp);
    for(i=0; i<nvars; i++)
    {
        range1temp=varrange1[i];
        range2temp=varrange2[i];

        switch (vartype[i])
        {
        case aoindex:
            type=1001;
            break;
        case moindex:
            type=1002;
            break;
        case moaindex:
            type=1003;
            break;
        case mobindex:
            type=1004;
            break;
        case indexx:
            type=1005;
            break;
        case laindex:
            type = 1006;
            break;
         case subindex:
            type = 1007;
            break;
        default:
            outerror(glno, "internal error in variable outputll.");
        }

        F77_NAME(add_index_table, ADD_INDEX_TABLE) (&type,&range1temp, &range2temp
                                                    );
    }
}
