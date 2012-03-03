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
syntax.cpp
Lei Wang
wang@qtp.ufl.edu
May 2004
******************************/
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "variable.h"
#include "instruct.h"
#include "keywdcnt.h"
#include "qcarray.h"
#include "errorhl.h"
#include "f77_name.h"
#include "f_types.h"

extern "C"
{
void F77_NAME(decode_cond_code, DECODE_COND_CODE)(const char *const c, f_int *fi);
};

//using namespace std;

extern int outArray[];
extern int outarray_size;
extern int idno;
extern string idname[];
extern int intconst[];
extern double floatconst[];
extern int intno;
extern int mx_array_dim;

int intnosyn=0;
int floatnosyn=0;
int idnosyn=0;
int endofoutArray=0;
int glno=0;
int pardocount=0;
int returncount=-1;
int docount=0;

extern VariableClass Variable;
extern QCArrayClass QCArray;
extern InstructCls Instruct;
/*ArrayClass QCArray;
*/

void instructionsyntax(int &lexpos, int &lineno);
//int assignsyntax(int &lexpos, int &lineno);
int factorsyntax(int &lexpos, int &datatype, int &lineno);
int primarysyntax(int &lexpos, int &datatype, int &lineno);
int termsyntax(int &lexpos, int &datatype, int &lineno);
int expsyntax(int &lexpos, int &datatype, int &lineno);
int conditionsyntax(int &lexpos, int &lineno);
void returnsilsyntax(int &lexpos, int &lineno);

int lexmovenext(int &lexpos)
{
    if(!endofoutArray)
    {
        if (lexpos<outarray_size)
        {
            lexpos++;
            return 1;
        }
        else //lexpos>=outArray_size
        {
            endofoutArray=1;
            cout<<"end of file unexpected."<<endl;
            return 0;
        }
    }
    else
        return 0;
}



void moveoverreturn (int &lexpos, int &lineno)
{
    if (outArray[lexpos]!=returnn && lexpos <outarray_size)
    {
        outerror(lineno, "@@@  unexpected char, ignored till end of line.");
        printf("Unexpected char: %s\n" , keywords[outArray[lexpos]]);
        if (outArray[lexpos]==id)
            printf("   '%s' ### unexpected char", idname[idnosyn].c_str());
        else
            printf("   '%s' $$$ unexpected char", keywords[outArray[lexpos]]);
        while (outArray[lexpos]!=returnn && !endofoutArray)
        {
            if (outArray[lexpos]==id)
                idnosyn++;
            lexmovenext(lexpos);
        }
    }

    while (outArray[lexpos]==returnn && !endofoutArray)
    {
        lineno++;
        glno=lineno;
        lexmovenext(lexpos);
    }
    if (lexpos>=outarray_size)
        endofoutArray=1;
}



int IsPreConstant(const int x)
{
    return (outArray[x]==norb||outArray[x]==nocc||
        outArray[x]==nvirt||outArray[x]==bocc||
        outArray[x]==eocc||outArray[x]==bvirt||
        outArray[x]==evirt||outArray[x]==naocc||
        outArray[x]==nbocc||outArray[x]==navirt||
        outArray[x]==nbvirt||outArray[x]==baocc||
        outArray[x]==bbocc||outArray[x]==eaocc||
        outArray[x]==ebocc||outArray[x]==bavirt||
        outArray[x]==bbvirt||outArray[x]==eavirt||
        outArray[x]==ebvirt ||
        outArray[x]==noccorb||
        outArray[x]==nvirtorb||outArray[x]==boccorb||
        outArray[x]==eoccorb||outArray[x]==bvirtorb||
        outArray[x]==evirtorb||outArray[x]==naoccorb||
        outArray[x]==nboccorb||outArray[x]==navirtorb||
        outArray[x]==nbvirt||outArray[x]==baocc||
        outArray[x]==bboccorb||outArray[x]==eaoccorb||
        outArray[x]==eboccorb||outArray[x]==bavirtorb||
        outArray[x]==bbvirtorb||outArray[x]==eavirtorb||
        outArray[x]==ebvirtorb ||outArray[x]==cc_iter ||
        outArray[x]==itrips || outArray[x]==itripe ||
        outArray[x]==ihess1 || outArray[x]==ihess2 ||
        outArray[x]==jhess1 || outArray[x]==jhess2 ||
        outArray[x]==cc_hist || outArray[x]==cc_beg ||
        outArray[x]==scf_iter || outArray[x]==scf_hist ||
        outArray[x]==scf_beg || outArray[x]==natoms ||
        outArray[x]==subb || outArray[x]==sube||
        outArray[x]==sip_sub_segsize || outArray[x]==sip_sub_occ_segsize||
        outArray[x]==sip_sub_virt_segsize || outArray[x]==sip_sub_ao_segsize);
}


void subindex_syntax(int &lexpos, int &lineno)
{
    int indexnameno=0;
    int indextype=outArray[lexpos]; //ao, moindex, index
    int superindex = 0;

    /* Syntax of the subindex declaration is: 
       subindex ii of i
       where ii is the name of the subindex
              i is the name of the superindex.
   */

    lexmovenext(lexpos);
    if (id!=outArray[lexpos])
    {
        outerror(lineno, "idname is missing.");
        indexnameno=0;
    }
    else
    {
        indexnameno=idnosyn++;    // subindex name.
        lexmovenext(lexpos);
    }

    /* Next element is the "of" keyword. */

    if (outArray[lexpos] != of_kwd)
        outerror(lineno, "subindex syntax error: of is missing.");
    else
        lexmovenext(lexpos);   // move pointer past "of".

    /* Check next element.  It must be an index of any type except subindex. */
 
    if (id!=outArray[lexpos])
    {
        outerror(lineno, "superindex idname is missing.");
    }
    else
    {
        superindex=Variable.getidnumber(idnosyn);    // superindex number.
        if (superindex == -1)
	{
           outerror(lineno, "superindex must be  a previously declared index.");
	}
        idnosyn++;
        lexmovenext(lexpos);
    }

    Variable.insert(indextype, indexnameno, superindex, 0);
    moveoverreturn(lexpos, lineno);
}

void indexsyntax(int &lexpos, int &lineno)
{
    int indexnameno=0;
    int indextype=outArray[lexpos]; //ao, moindex, index
    int indexrange1=0, indexrange2=0;
    /*int result=0;*/

    if (outArray[lexpos] == subindex)
    {
       subindex_syntax(lexpos, lineno);
       return;
    }

    lexmovenext(lexpos);
    if (id!=outArray[lexpos])
    {
        outerror(lineno, "idname is missing.");
        indexnameno=0;
    }
    else
    {
        indexnameno=idnosyn++;
        lexmovenext(lexpos);
    }

    if (assign!=outArray[lexpos])
        outerror(lineno, "= is missing.");
    else
        lexmovenext(lexpos);

    if (IsPreConstant(lexpos))
    {
        indexrange1=-outArray[lexpos];
        lexmovenext(lexpos);
    }
    else if (outArray[lexpos]==iconst)
    {
        indexrange1=intconst[intnosyn];
        intnosyn++;
        lexmovenext(lexpos);
    }
    else
        outerror(lineno, "start range is missing.");

    if (comma!=outArray[lexpos])
        outerror(lineno, "comma is missing.");
    else
        lexmovenext(lexpos);

    if (IsPreConstant(lexpos))
    {
        indexrange2=-outArray[lexpos];
        lexmovenext(lexpos);
    }
    else if (outArray[lexpos]==iconst)
    {
        lexmovenext(lexpos);
        indexrange2=intconst[intnosyn];
        intnosyn++;
    }
    else
    {
        outerror(lineno, "end range is missing.");
        indexrange2=1;
    }

    Variable.insert(indextype, indexnameno, indexrange1, indexrange2);
    moveoverreturn(lexpos, lineno);
}


int IsAssignOp(const int x)
{
    return (outArray[x]==assign||outArray[x]==addassign
          ||outArray[x]==subassign||outArray[x]==multassign
          ||outArray[x]==divassign);
}



int IsCon(const int x)
{
    return (outArray[x]==ge||outArray[x]==le
          ||outArray[x]==gt||outArray[x]==lt
          ||outArray[x]==eq||outArray[x]==ne);
}



/*====================================
The function IsDeclaration returns true if x is a proc, false otherwise.
====================================*/
int IsProc(const int x)
{
    return (outArray[x]==proc);
}



/*====================================
The function IsDeclaration returns true if x is a proc, false otherwise.
====================================*/
int IsAlsoPardo(const int x)
{
    return (outArray[x]==also);
}


/*====================================
The function IsDeclaration returns true if x is a declaration, false otherwise.
====================================*/
int IsDeclaration(const int x)
{
    return (outArray[x]==aoindex||outArray[x]==moindex||
            outArray[x]==moaindex||outArray[x]==mobindex||
            outArray[x]==laindex || outArray[x]==subindex ||
            outArray[x]==served||outArray[x]==temp||
            outArray[x]==distributed||outArray[x]==staticc||
            outArray[x]==local||
            outArray[x]==scalar||outArray[x]==indexx);
}

/*====================================
The function IsIndex returns true if x is an index type, false otherwise.
====================================*/
int IsIndex(const int x)
{
    return (x==aoindex||x==moindex||
            x==moaindex||x==mobindex||
            x==laindex ||
            x==indexx);
}

int IsScalar(const int x)
{
    return (outArray[x]==iconst||outArray[x]==fconst||
           (x+1<outarray_size&&outArray[x]==sub
            &&(outArray[x+1]==iconst||outArray[x+1]==fconst))||
           (x+1<outarray_size&&outArray[x]==id&&outArray[x+1]!=lparen));
}



/*====================================
The function IsInstruction returns true if x is a instruction, false otherwise.
====================================*/
int IsInstruction(const int x)
{
    return (outArray[x]==setindex||outArray[x]==put||
            outArray[x]==execute ||outArray[x]==doo||
            outArray[x]==call||outArray[x]==iff ||
            outArray[x]==create ||outArray[x]==deletee||
            outArray[x]==pardo ||outArray[x]==get||
            outArray[x]==prepare ||outArray[x]==request||
            outArray[x]==destroye ||outArray[x]==writee||
            outArray[x]==collective||outArray[x]==computeint||
            (outArray[x]==id&&outArray[x+1]==lparen)||
            (outArray[x]==id&&outArray[x+1]==assign)||
            (outArray[x]==id&&outArray[x+1]==addassign)||
            (outArray[x]==id&&outArray[x+1]==subassign)||
            (outArray[x]==id&&outArray[x+1]==multassign)||
            (outArray[x]==id&&outArray[x+1]==divassign)||
            outArray[x]==returnsil||
            outArray[x]==exitt||
            outArray[x]==cycle||
            outArray[x]==allocatee||
            outArray[x]==deallocatee || outArray[x]==prequest ||
            outArray[x]== where);
}



int scalarsyntax(int &lexpos, int &idname)
{
    double value;
    int noerrorlocal=1;
    int negtive=0;

    if (id==outArray[lexpos])
    {
        idname=idnosyn++;
        lexmovenext(lexpos);
    }
    else    //constant
    {
        if (sub==outArray[lexpos])
        {
            negtive=1;
            lexmovenext(lexpos);
        }

        if (iconst==outArray[lexpos])
        {
            value=double(intconst[intnosyn++]);
            if (negtive)
                value=-value;
            lexmovenext(lexpos);
            idname=-QCArray.inserttempscalar(value);
            return noerrorlocal;
        }
        else if (fconst==outArray[lexpos])
        {
            value=floatconst[floatnosyn++];
            if (negtive)
                value=-value;
            lexmovenext(lexpos);
            idname=-QCArray.inserttempscalar(value);
            return noerrorlocal;
        }
    }
    return noerrorlocal;
}



int allocatearraysyntax(int &lexpos, int &lineno, int &arraynameid,
                      int &nindex, int indexarray[])
{
    int noerrorlocal=1;

    if (id==outArray[lexpos])
    {
        arraynameid=idnosyn++;
        lexmovenext(lexpos);
    }
    else
    {
        outerror(lineno, "Arrayname is missing.");
        noerrorlocal=0;
    }

    if (lparen!=outArray[lexpos])
    {
        nindex=0;
        noerrorlocal = 0;
        return noerrorlocal;
    }
    else
        lexmovenext(lexpos);

    nindex=0;
    while ((id==outArray[lexpos] ||mult==outArray[lexpos])&& !endofoutArray)
    {
        if (nindex<mx_array_dim)
        {
            if (id==outArray[lexpos])
            {
                indexarray[nindex]=idnosyn;
                idnosyn++;
            }
            else
                indexarray[nindex]=WILDCARD_INDICATOR;
        }

	nindex++;

        lexmovenext(lexpos);

        if (comma==outArray[lexpos])
        {
            lexmovenext(lexpos);
        }
        else if (rparen==outArray[lexpos])
        {
            lexmovenext(lexpos);
            break;
        }
        else
        {
            outerror(lineno, "error in array declaration.");
            noerrorlocal=0;
            break;
        }
    }

    if (nindex>mx_array_dim)
    {
        outerror(lineno, "Warning: Too many indices detected.");
        printf("Only the first %d are used\n", mx_array_dim);
        printf("%d=nindex, mx_array_dim = %d\n", nindex,mx_array_dim);
        printf("%s\n", idname[arraynameid].c_str());
        for (int i=arraynameid+1; i<idnosyn; i++)
            printf("%s, ", idname[i].c_str());
        printf("\n");
        noerrorlocal=0;
        nindex=mx_array_dim;
    }
    return noerrorlocal;
}

int commonarraysyntax(int &lexpos, int &lineno, int &arraynameid,
                      int &nindex, int indexarray[])
{
    int noerrorlocal=1;

    if (IsScalar(lexpos))
    {
        nindex=0;
        return scalarsyntax(lexpos, arraynameid);
    }

    if (id==outArray[lexpos])
    {
        arraynameid=idnosyn++;
        lexmovenext(lexpos);
    }
    else
    {
        outerror(lineno, "Arrayname is missing.");
        noerrorlocal=0;
    }

    if (lparen!=outArray[lexpos])
    {
        nindex=0;
        return noerrorlocal;
/*        outerror(lineno, "Array ( is missing.");//idname[arraynameid]
        noerrorlocal=0;*/
    }
    else
        lexmovenext(lexpos);

    nindex=0;
    while (id==outArray[lexpos] && !endofoutArray)
    {
        if (nindex<mx_array_dim)
        {
            indexarray[nindex]=idnosyn;
        }

        nindex++;

        idnosyn++;
        lexmovenext(lexpos);

        if (comma==outArray[lexpos])
            lexmovenext(lexpos);
        else if (rparen==outArray[lexpos])
        {
            lexmovenext(lexpos);
            break;
        }
        else
        {
            outerror(lineno, "error in array declaration.");
            noerrorlocal=0;
            break;
        }
    }

    if (nindex>mx_array_dim)
    {
        outerror(lineno, "Warning: Too many indices detected.");
        printf("Only the first %d are used.\n", mx_array_dim);
        printf("%d=nindex\n", nindex);
        printf("%s\n", idname[arraynameid].c_str());
        for (int i=arraynameid+1; i<idnosyn; i++)
            printf("%s, ", idname[i].c_str());
        printf("\n");
        noerrorlocal=0;
        nindex=mx_array_dim;
    }
    return noerrorlocal;
}



void scalardecsyntax(int &lexpos, int &lineno)
{
    int indexnameno=0;
    double value1=0;
    int negtive=0;

    lexmovenext(lexpos);
    if (id!=outArray[lexpos])
        outerror(lineno, "idname is missing.");
    else
    {
        indexnameno=idnosyn++;
        lexmovenext(lexpos);
    }

    if (assign==outArray[lexpos])
    {
        lexmovenext(lexpos);

        if (sub==outArray[lexpos])
        {
            negtive=1;
            lexmovenext(lexpos);
        }

        if (fconst==outArray[lexpos])
        {
            value1=floatconst[floatnosyn];
            if (negtive)
                value1=-value1;
            floatnosyn++;
            lexmovenext(lexpos);
        }
        else if (iconst==outArray[lexpos])
        {
            lexmovenext(lexpos);
            value1=double(intconst[intnosyn]);
            if (negtive)
                value1=-value1;
            intnosyn++;
        }
        else
            outerror(lineno, "value is missing.");
    }

    QCArray.insertscalar(indexnameno, value1);

    moveoverreturn (lexpos, lineno);
}



void arraysyntax(int &lexpos, int &lineno)
/*====================================
The function declarationsyntax anylizes the declaration of token array stored
in outArray.
201 (served integral array)
202 (staticc array)
203 (distributed array)
204 (temp array)
205 (scalar)
206 (local array)
====================================*/
{
    int arraynameid;
    int nindex=0;
    int indexarray[mx_array_dim];

    int arraytype=outArray[lexpos++];

    commonarraysyntax(lexpos, lineno, arraynameid, nindex, indexarray);

    QCArray.insert(arraytype, arraynameid, nindex, indexarray);

    moveoverreturn (lexpos, lineno);
}



void exitcyclesyntax(int &lexpos, int &lineno)
{
    if (docount==0)
        outerror(lineno, "exit or cycle not allowed outside do loop.");
    if (exitt==outArray[lexpos]||cycle==outArray[lexpos])
        Instruct.insert(outArray[lexpos], 0);
    else
        outerror(lineno, "unknow error in exitcyclesyntax.");
    lexmovenext(lexpos);

    moveoverreturn(lexpos, lineno);
}



void callsyntax(int &lexpos, int &lineno)
{
    int idno=0;
    lexmovenext(lexpos);

    if (id!=outArray[lexpos])
        outerror(lineno, "function name is missing.");
    else
    {
        idno=idnosyn++;
        lexmovenext(lexpos);
    }

    Instruct.insert(call, idno);
    moveoverreturn(lexpos, lineno);
}



void executesyntax(int &lexpos, int &lineno)
{
    int arraynameid, nindex;
    int indexarray[mx_array_dim];
    int fnnameid=0;
    int idno=0;
    int idno1=0;

    lexmovenext(lexpos);
    if (id!=outArray[lexpos])
        outerror(lineno, "function name is missing.");
    else
    {
        fnnameid=idnosyn++;
        lexmovenext(lexpos);
        if (id != outArray[lexpos])
        {
            Instruct.insert(execute, fnnameid, -1);
            moveoverreturn(lexpos, lineno);
            return;
        }

    }

    if (id ==outArray[lexpos] && lparen ==outArray[lexpos+1])
    {
        if(!commonarraysyntax(lexpos, lineno, arraynameid, nindex, indexarray))
        {
            outerror(lineno, "Array name is missing.");
            return;
        }

        Instruct.insert(execute, fnnameid, arraynameid, nindex, indexarray);
    }
    else if (id==outArray[lexpos])
    {
        idno=idnosyn++;
        lexmovenext(lexpos);

        if (id==outArray[lexpos])
        {
            idno1=idnosyn++;
            lexmovenext(lexpos);
            /*insert 2 arguements*/
            Instruct.insert(execute, fnnameid, idno, idno1);
        }
        else if (fconst==outArray[lexpos])
        {
           scalarsyntax(lexpos, idno1);
           Instruct.insert(execute, fnnameid, idno, idno1);
         }
        else                            /*insert one argument*/
            Instruct.insert(execute, fnnameid, idno);

    }

    moveoverreturn(lexpos, lineno);
}



void createsyntax(int &lexpos, int &lineno)
{
    int idno=0;
    int save_lexpos;
    int arrayno;
    int nindex;
    int indexarray[mx_array_dim];

    save_lexpos = lexpos;

    // First, check for the "partial create syntax.

    lexmovenext(lexpos);
    if (allocatearraysyntax(lexpos, lineno, arrayno, nindex, indexarray))
    {
       // We have a valid partial create syntax.
       // Generate the proper instruction.

       Instruct.insertallocate(create, arrayno, nindex, indexarray);
    }
    else
    {
       lexpos = save_lexpos;   // back up to the create's lexpos position.
       idnosyn--;
       lexmovenext(lexpos);

       if (id!=outArray[lexpos])
           outerror(lineno, "array name is missing.");
       else
       {
           idno=idnosyn++;
           lexmovenext(lexpos);
       }

       Instruct.insert(create, idno);
    }

    moveoverreturn(lexpos, lineno);
}

void deletesyntax(int &lexpos, int &lineno)
{
    int idno=0;
    int save_lexpos;
    int arrayno;
    int nindex;
    int indexarray[mx_array_dim];

    save_lexpos = lexpos;

    // First, check for the "partial delete syntax.

    lexmovenext(lexpos);
    if (allocatearraysyntax(lexpos, lineno, arrayno, nindex, indexarray))
    {
       // We have a valid partial delete syntax.
       // Generate the proper instruction.

       Instruct.insertallocate(deletee, arrayno, nindex, indexarray);
    }
    else
    {
       lexpos = save_lexpos;
       idnosyn--;
       lexmovenext(lexpos);

       if (id!=outArray[lexpos])
           outerror(lineno, "array name is missing.");
       else
       {
           idno=idnosyn++;
           lexmovenext(lexpos);
       }

       Instruct.insert(deletee, idno);
    }
    moveoverreturn(lexpos, lineno);
}

void destroysyntax(int &lexpos, int &lineno)
{
    int idno=0;

    lexmovenext(lexpos);
    if (id!=outArray[lexpos])
        outerror(lineno, "function name is missing.");
    else
    {
        idno=idnosyn++;
        lexmovenext(lexpos);
    }

    Instruct.insert(destroye, idno);
    moveoverreturn(lexpos, lineno);
}


/*====================================
The function declarationsyntax anylizes the declaration of token array stored
in outArray.
====================================*/
void declarationsyntax(int &lexpos, int &lineno)
{
    if (endofoutArray) return;
    if (aoindex==outArray[lexpos]||moindex==outArray[lexpos]||
        moaindex==outArray[lexpos]||mobindex==outArray[lexpos]||
        outArray[lexpos] == laindex || outArray[lexpos] == subindex ||
        indexx==outArray[lexpos])
        indexsyntax(lexpos, lineno);
    else if (scalar==outArray[lexpos])
        scalardecsyntax(lexpos, lineno);
    else if (served==outArray[lexpos]||temp==outArray[lexpos]
             ||distributed==outArray[lexpos]||staticc==outArray[lexpos]
             ||local==outArray[lexpos])
        arraysyntax(lexpos, lineno);
    else
        outerror(lineno, "error in declaration.");
}



void putpreparesyntax(int &lexpos, int &lineno)
{
    /*int i=0, j=0;*/
    int instruction=0;
    int arraynameid, nindex, indexarray[mx_array_dim];
    int array1nameid, nindex1, indexarray1[mx_array_dim];

    if (put==outArray[lexpos])
        instruction=putadd;
    else if (prepare==outArray[lexpos])
        instruction=prepareadd;
    else
        outerror(lineno, "unknown error.");
    lexmovenext(lexpos);

    if (!commonarraysyntax(lexpos, lineno, arraynameid, nindex, indexarray))
    {
        outerror(lineno, "Array name is missing.");
        moveoverreturn(lexpos, lineno);
        return;
    }

    if (instruction==prepareadd)
    {
        if (addassign==outArray[lexpos])
            lexmovenext(lexpos);
        else if (assign==outArray[lexpos])
        {
            instruction=prepare;
            lexmovenext(lexpos);
        }
        else
            outerror(lineno, "'+=' or '=' expected.");
    }
    else    /*put*/
    {
        if (addassign==outArray[lexpos])
            lexmovenext(lexpos);
        else if (assign==outArray[lexpos])
        {
            instruction=put;
            lexmovenext(lexpos);
        }
        else
            outerror(lineno, "'+=' or '=' expected.");
    }

    if (commonarraysyntax(lexpos, lineno, array1nameid, nindex1, indexarray1))
    {
        if (instruction==put)
        Instruct.insert(put, arraynameid, nindex, indexarray,
                                    array1nameid, nindex1, indexarray1);
        else if (instruction==putadd)
        Instruct.insert(putadd, arraynameid, nindex, indexarray,
                                    array1nameid, nindex1, indexarray1);
        else if (instruction==prepare)
        Instruct.insert(prepare, arraynameid, nindex, indexarray,
                                    array1nameid, nindex1, indexarray1);
        else if (instruction==prepareadd)
        Instruct.insert(prepareadd, arraynameid, nindex, indexarray,
                                    array1nameid, nindex1, indexarray1);
        else
            outerror(lineno, "Internal error in putprepare.");
    }
    else
        outerror(lineno, "second array is missing.");

    moveoverreturn(lexpos, lineno);
}



void collectivesyntax(int &lexpos, int &lineno)
{
    /*int i=0, j=0;*/
    int arraynameid, nindex, indexarray[mx_array_dim];
    int array1nameid, nindex1, indexarray1[mx_array_dim];

    lexmovenext(lexpos);

    if (!commonarraysyntax(lexpos, lineno, arraynameid, nindex, indexarray))
        outerror(lineno, "Array name is missing.");

    if (addassign==outArray[lexpos])
        lexmovenext(lexpos);
    else
        outerror(lineno, "'+=' is missing");

    if (commonarraysyntax(lexpos, lineno, array1nameid, nindex1, indexarray1))
        Instruct.insert(collective, arraynameid, nindex, indexarray,
                                    array1nameid, nindex1, indexarray1);
    else
        outerror(lineno, "second array is missing.");

    moveoverreturn(lexpos, lineno);
}


/*
int assignsyntax(int &lexpos, int &lineno)
{
    int op=0;
    int stackno1, stackno2;

    outerror(lineno, "assign is not available right now.");
    stackno1=expsyntax(lexpos, lineno);


    if (IsAssignOp(lexpos))
    {
        op=outArray[lexpos];
        stackno2=assignsyntax(lexpos, lineno);
        lexmovenext(lexpos);
    }

    Instruct.inserts (assign, op, op, op);
    return stackno1;
}
*/


int factorsyntax(int &lexpos, int &datatype, int &lineno)
{
    int op=0;
    int ivalue;
    double fvalue;
    int stackno=0;
	int temptype=iconst;

    if (outArray[lexpos] >= norb && outArray[lexpos] < indexx)
    {
       datatype = iconst;
       stackno=Instruct.inserts (symbolic_const, outArray[lexpos], 0, 0);
       lexmovenext(lexpos);
    }
    else if (iconst==outArray[lexpos])
    {
		datatype=iconst;
        ivalue=intconst[intnosyn++];
        stackno=Instruct.inserts (iconst, ivalue, 0, 0);
        lexmovenext(lexpos);
    }

    else if (fconst==outArray[lexpos])
    {
		datatype=fconst;
        /*outerror(lineno, "possible float error.");*/
        fvalue=floatconst[floatnosyn++];
		stackno=Instruct.insertsf (fconst, fvalue, 0, 0);
        lexmovenext(lexpos);
    }

    else if (id==outArray[lexpos])
    {
        op=idnosyn;
        stackno=Variable.getidnumber(idnosyn);
        if (stackno<0)
		{
			stackno=QCArray.getidnumber(idnosyn);
			if (stackno<0)
			{
				outerror(lineno, "variable not declared.");
			}
			else if(QCArray.gettype(stackno)!=scalar)
			{
				outerror(lineno, "array is not expected.");
			}
			else
			{
				datatype=fconst;
				stackno=Instruct.inserts (idf, stackno, 0, 0);
			}

		}
		else
		{
			datatype=iconst;
	        stackno=Instruct.inserts(id, stackno, 0, 0);
		}

		idnosyn++;
        lexmovenext(lexpos);
    }

    else if (lparen==outArray[lexpos])
    {
        lexmovenext(lexpos);

        stackno=expsyntax(lexpos, temptype, lineno);
		datatype=temptype;

        if (rparen!=outArray[lexpos])
            outerror(lineno, " %%% unexpected char.");
        else
            lexmovenext(lexpos);
    }
    else
    {
        outerror(lineno, "^^^ unexpected char.");
        printf("Unexpected char: %s\n",keywords[outArray[lexpos]]);
    }

    return stackno;
}



int primarysyntax(int &lexpos, int &datatype, int &lineno)
{
    int op=0;
    int stackno1, stackno2;

    if (sub==outArray[lexpos])                          //-
    {
        op = sub;
        lexmovenext(lexpos);
    }
    stackno1=factorsyntax(lexpos, datatype, lineno);

    if (op==sub)
    {
		if (datatype==iconst)
		{
        stackno2=Instruct.inserts (iconst, -1, 0, 0);
        Instruct.inserts (mult, stackno1, stackno1, stackno2);
		}
		else
		{
			stackno2=Instruct.insertsf (fconst, -1, 0, 0);
			Instruct.inserts (multf, stackno1, stackno1, stackno2);
		}
    }
    return stackno1;
}



int termsyntax(int &lexpos, int &datatype, int &lineno)
{
    int op=0;
    int stackno1, stackno2;
	int datatype2=iconst;

    stackno1=primarysyntax(lexpos, datatype, lineno);

    while (mult==outArray[lexpos]||divide==outArray[lexpos])       //*, /
    {
        op = outArray[lexpos];
        lexmovenext(lexpos);
        stackno2=primarysyntax(lexpos, datatype2, lineno);
		if (datatype2!=datatype)
			outerror(lineno, "data type mismatch.");
		if (datatype==iconst)
			Instruct.inserts (op, stackno1, stackno1, stackno2);
		else
			Instruct.inserts (multf-mult+op, stackno1, stackno1, stackno2);
    }
    return stackno1;
}



int expsyntax(int &lexpos, int &datatype, int &lineno)
{
    int op=0;
    int stackno1, stackno2;
	int datatype2=iconst;

    stackno1=termsyntax(lexpos, datatype, lineno);

    while (add==outArray[lexpos]||sub==outArray[lexpos])       //+, -
    {
        op = outArray[lexpos];
        lexmovenext(lexpos);
        stackno2=termsyntax(lexpos, datatype2, lineno);
		if (datatype2!=datatype)
			outerror(lineno, "data type mismatch.");
		if(datatype==iconst)
			Instruct.inserts (op, stackno1, stackno1, stackno2);
		else
			Instruct.inserts (addf-add+op, stackno1, stackno1, stackno2);
    }
    return stackno1;
}



int conditionsyntax(int &lexpos, int &lineno)
{
    int op=0;
    int stackno1, stackno2;
    int datatype=iconst;
	int datatype2=iconst;

    stackno1=expsyntax(lexpos, datatype, lineno);

    if (IsCon(lexpos))       //==, >=,<=, !=, >,<
    {
        op = outArray[lexpos];
        lexmovenext(lexpos);
        stackno2=expsyntax(lexpos, datatype2, lineno);
		if (datatype!=datatype2)
			outerror(lineno, "data type does not match.");
    }
    else
    {
		if (datatype==iconst)
			stackno2=Instruct.inserts(iconst, 0, 0, 0);
		else
			stackno2=Instruct.inserts(fconst, 0, 0, 0);
        op=ne;
    }
	if (datatype==iconst)
		Instruct.inserts(op, stackno1, stackno1, stackno2);
	else
		Instruct.inserts(gtf-gt+op, stackno1, stackno1, stackno2);

	return stackno1;
}



void ifsyntax(int &lexpos, int &lineno)
{
    int stackno;
    int insif, insfalse, insend;

    lexmovenext(lexpos);                                    //if

    stackno=conditionsyntax(lexpos, lineno);
    insif=Instruct.insertif(iff, stackno);

    moveoverreturn(lexpos, lineno);

    while (IsInstruction(lexpos))
    {
        instructionsyntax(lexpos, lineno);
    }

    if (elsee==outArray[lexpos])                            //else
    {
        lexmovenext(lexpos);
        moveoverreturn(lexpos, lineno);

        insfalse=Instruct.insertif(elsee, 0);  //goto
        Instruct.reinsert(insif, insfalse+1);

        while (IsInstruction(lexpos))
        {
            instructionsyntax(lexpos, lineno);
        }
        insend=Instruct.getinsno();
        Instruct.reinsert(insfalse, insend);
    }
    else
    {
        insend=Instruct.getinsno();
        Instruct.reinsert(insif, insend);
    }

    if (endiff!=outArray[lexpos])                           //elseif
    {
        outerror(lineno, "endif missing.");
    }
    else
    {
        lexmovenext(lexpos);
        moveoverreturn(lexpos, lineno);
    }
}

void in_syntax_check(int &lexpos, int &lineno, int subindex_idnosyn)
{
    int superindex, subindex_varid;
    lexmovenext(lexpos);  // eat the "in" keyword 
    if (outArray[lexpos] == id)  // next element must be the superindex 
    {
       /* Check to make sure the index is the superindex of this 
          subindex. */

       superindex = Variable.getidnumber(idnosyn); 
       subindex_varid   = Variable.getidnumber(subindex_idnosyn);
              
       /* The superindex is stored in the range1 field of the 
          variable table entry of the subindex.  It should match the
          superindex declared in this pardo. */
 
       if (Variable.getrange1(subindex_varid) != superindex)
          outerror(lineno, "subindex and superindex do not correspond.");
                 
    }
    else
    {
       outerror(lineno, "subindex missing superindex");
    }

    /* Done with the "in" syntax check.  Move to next index on pardo. */

    idnosyn++;
    lexmovenext(lexpos);  
    if (outArray[lexpos] == comma) lexmovenext(lexpos);  
}

void pardoline(int &lexpos, int &lineno, int &nindex, int indexarray[],
               int in_flag[])
{
    int i;

    nindex=0;

    for (i = 0; i < mx_array_dim; i++)
      in_flag[i] = 0;
 
    while (!endofoutArray&& id==outArray[lexpos])
    {
        if (nindex<mx_array_dim)
            indexarray[nindex]=idnosyn;

        nindex++;

        idnosyn++;
        lexmovenext(lexpos);

        if (comma==outArray[lexpos])
            lexmovenext(lexpos);
        else if (outArray[lexpos] == in_kwd)   // "pardo ... ii in i" syntax.
        {
           in_syntax_check(lexpos, lineno, indexarray[nindex-1]);
           in_flag[nindex-1] = 1;
        }   
        else
            break;
    }

    if (nindex>mx_array_dim)
    {
        outerror(lineno, "Warning: Too many indices detected.");
        printf("Only the first %d are used.\n",mx_array_dim);
	printf("%d=nindex\n", nindex);
	printf("All the indices are: ");
        for (int ii=indexarray[0]; ii<idnosyn; ii++)
            printf("%s, ", idname[ii].c_str());
        printf("\n");
        nindex=mx_array_dim;
    }

    if (nindex==0)
        outerror(lineno, "id is missing.");
    Instruct.insert(pardo, nindex, indexarray, in_flag);

    moveoverreturn(lexpos, lineno);
}


void pardosyntax(int &lexpos, int &lineno)
{
    int i;
    int nalsopardo=0;
    int **indexa;
    int *nindex;
    //int indexarray[mx_array_dim], indexarray1[mx_array_dim];
    int in_flag[mx_array_dim];

    if (pardocount>=1)
        outerror(lineno, " nested pardo loop is not allowed.");

    pardocount++;

    lexmovenext(lexpos);                            //pardo
    nalsopardo++;

    //count pardo numbers
    i=lexpos;
    while (i+1<outarray_size && endpardo!=outArray[i])
    {
        if (also==outArray[i]&& pardo==outArray[i+1])
            nalsopardo++;
        i++;
    }

    if (i+1==outarray_size)
        outerror(lineno, "no endpardo.");

    indexa=new int * [nalsopardo];
    if (!indexa)
    {
        printf("error in malloc, exit.");
        exit(1);
    }
    nindex=new int [nalsopardo];
    if (!nindex)
    {
        printf("error in malloc, exit.");
        exit(1);
    }

    
    for (i=0; i<nalsopardo; i++)
    {
        indexa[i]=new int [mx_array_dim];
        if (!indexa[i])
        {
            printf("error in malloc, exit.");
            exit(1);
        }
    }

    i=0;
    //first pardo
    pardoline(lexpos, lineno, nindex[i], indexa[i], in_flag);

    while(!endofoutArray && IsInstruction(lexpos) )
    {
        instructionsyntax(lexpos, lineno);
    }

    while(!endofoutArray && IsAlsoPardo(lexpos))
    {
        lexmovenext(lexpos);

        if (pardo!=outArray[lexpos])
        {
            outerror(lineno,  "pardo is missing.");
        }
        else
            lexmovenext(lexpos);
        i++;        //also pardo increment

        pardoline(lexpos, lineno, nindex[i], indexa[i], in_flag);

        while(!endofoutArray && IsInstruction(lexpos) )
        {
            instructionsyntax(lexpos, lineno);
        }
    }


    if (endpardo!=outArray[lexpos])
    {
        outerror(lineno,  "endpardo is missing.");
        printf("%s\n", idname[idnosyn].c_str());
    }
    else
        lexmovenext(lexpos);


    for (i=0; i<nalsopardo; i++)
    {
        for (int j=0; j<nindex[i]; j++)
        {
            if (!endofoutArray&& id==outArray[lexpos])
            {
                lexmovenext(lexpos);
                if (idname[idnosyn++]!=idname[indexa[i][j]])
                    outerror(lineno, "syntax1: index does not match.");
            }
            else
                outerror(lineno, "syntax2: index does not match.");

            in_flag[j] = 0;
            if (comma==outArray[lexpos])
                lexmovenext(lexpos);
            else if (outArray[lexpos] == in_kwd) 
            {
               // This is the "endpardo... ii in i" form.

               // idnosyn has already been bumped when checking the subindex.
               // To get the proper value (i. e. pointing to the subindex, 
               // we must decrement idnosyn.    

               in_syntax_check(lexpos, lineno, idnosyn-1); 
               in_flag[j] = 1;
            }
        }
    }

    Instruct.insert(endpardo, nindex[0], indexa[0], in_flag);

    moveoverreturn(lexpos, lineno);
    pardocount --;

    for (i=0; i<nalsopardo; i++)
    {
        delete [](indexa[i]);
    }
    delete [](indexa);
}

void wheresyntax(int &lexpos, int &lineno)
{
   int index1, index2, icond;
   f_int cond;

   lexmovenext(lexpos);
   if (outArray[lexpos] == id)
   {
      index1 = idnosyn++;
      lexmovenext(lexpos);

      if (IsCon(lexpos))
      {
         icond = outArray[lexpos];
         F77_NAME(decode_cond_code, DECODE_COND_CODE)(keywords[icond], &cond);
         lexmovenext(lexpos);
         if (outArray[lexpos] == id)
         {
            index2 = idnosyn++;
            lexmovenext(lexpos);
            Instruct.insert(where, index1, index2, cond);
         }
         else
            outerror(lineno, "Second operand of WHERE must be an index");
      }
      else
         outerror(lineno, "Invalid conditional in WHERE");
   }
   else
      outerror(lineno, "First operand in WHERE must be an index");

   moveoverreturn(lexpos, lineno);
}


void dosyntax(int &lexpos, int &lineno)
{
    int variableno=0;
    int in_flag = 0;
    lexmovenext(lexpos);

    docount++;

    if (id!=outArray[lexpos])
        outerror(lineno, "id is missing.");
    else
    {
        variableno=idnosyn;
        idnosyn++;
        lexmovenext(lexpos);

        // Check for the "in" syntax of the subindex.

        if (outArray[lexpos] == in_kwd)
        {
           in_syntax_check(lexpos, lineno, idnosyn-1);
           in_flag = 1;
        }

        Instruct.insert(doo, variableno, in_flag);
    }

    moveoverreturn(lexpos, lineno);


    while(!endofoutArray&&enddo!=outArray[lexpos])
    {

        while(IsInstruction(lexpos) && !endofoutArray)
        {
            instructionsyntax(lexpos, lineno);
        }

        if (enddo!=outArray[lexpos])
        {
            moveoverreturn(lexpos, lineno);
        }
    }

    if (enddo!=outArray[lexpos])
            outerror(lineno,  "enddo missing.");
    else
        lexmovenext(lexpos);

    if (id!=outArray[lexpos])
    {
        outerror(lineno,  "id is missing.");
    }
    else
    {
        if (idname[variableno]!=idname[idnosyn])
            outerror(lineno, "index number does not match.");
        else
        {
            idnosyn++;
            lexmovenext(lexpos);

            // Check for "in" syntax.
     
            in_flag = 0;
            if (outArray[lexpos] == in_kwd)
            {
               in_syntax_check(lexpos, lineno, idnosyn-1);
               in_flag = 1; 
            }

            Instruct.insert (enddo, variableno, in_flag);
        }
    }

    moveoverreturn(lexpos, lineno);
    docount--;
}



void getcomputesyntax(int&lexpos, int&lineno)
{
    int arrayno;
    int nindex;
    int indexarray[mx_array_dim];
    int instructionno=0;

    instructionno=outArray[lexpos];
    lexmovenext(lexpos);

    commonarraysyntax(lexpos, lineno, arrayno, nindex, indexarray);

    Instruct.insert(instructionno, arrayno, nindex, indexarray);

    moveoverreturn(lexpos, lineno);
}



void allocatesyntax(int&lexpos, int&lineno)
{
    int arrayno;
    int nindex;
    int indexarray[mx_array_dim];
    int instructionno=0;

    instructionno=outArray[lexpos];
    lexmovenext(lexpos);

    allocatearraysyntax(lexpos, lineno, arrayno, nindex, indexarray);

    Instruct.insertallocate(instructionno, arrayno, nindex, indexarray);

    moveoverreturn(lexpos, lineno);
}



void requestsyntax(int&lexpos, int&lineno)
{
    int arrayno;
    int nindex;
    int indexarray[mx_array_dim];
    int variableno=0;

    lexmovenext(lexpos);

    commonarraysyntax(lexpos, lineno, arrayno, nindex, indexarray);

    if (id!=outArray[lexpos])
        outerror(lineno, "id is missing.");
    else
    {
        variableno=idnosyn;
        idnosyn++;
        lexmovenext(lexpos);
        Instruct.insert(request, variableno, arrayno, nindex, indexarray);
    }


    moveoverreturn(lexpos, lineno);
}


void prequest_syntax(int&lexpos, int&lineno)
{
    int arrayno, arrayno2;
    int nindex, nindex2;
    int indexarray[mx_array_dim], indexarray2[mx_array_dim];

    lexmovenext(lexpos);
    commonarraysyntax(lexpos, lineno, arrayno, nindex, indexarray);

    if (outArray[lexpos] != assign)
    {
       outerror(lineno, "Expected = is missing.");
       printf("outArray[lexpos] = %d\n",outArray[lexpos]);
       moveoverreturn(lexpos, lineno);
       return;
    }

    lexmovenext(lexpos);
    commonarraysyntax(lexpos, lineno, arrayno2, nindex2, indexarray2);

    Instruct.insert(prequest, arrayno, nindex, indexarray,
                     arrayno2, nindex2, indexarray2);
    moveoverreturn(lexpos, lineno);
}

/*====================================
contraction and summation, subtraction, increment, decrement.
====================================*/
void operationidsyntax(int &lexpos, int &lineno)
{
    int array1no, array2no, array3no;
    int nindex1, nindex2, nindex3;
    int indexarray1[mx_array_dim], indexarray2[mx_array_dim], indexarray3[mx_array_dim];
    int op;
    int i;

    commonarraysyntax(lexpos, lineno, array1no, nindex1, indexarray1);

    /*    array1no=idnosyn;
    idnosyn++;
    lexmovenext(lexpos);
*/

    if (addassign==outArray[lexpos]||
        subassign==outArray[lexpos]||
        multassign==outArray[lexpos]||
        divassign==outArray[lexpos])
    {
        if (addassign==outArray[lexpos])
            op=add;
        else if (subassign==outArray[lexpos])
            op=sub;
        else if (multassign==outArray[lexpos])
            op=mult;
        else
            op=divide;
        lexmovenext(lexpos);
        commonarraysyntax(lexpos, lineno, array2no, nindex2, indexarray2);

        array3no=array1no;
        nindex3=nindex1;
        for(i=0;i<nindex3; i++)
            indexarray3[i]=indexarray1[i];
        Instruct.insert(op, array1no, nindex1, indexarray1,
            array3no, nindex3, indexarray3,
            array2no, nindex2, indexarray2);
    }
    else if (assign==outArray[lexpos])
    {
        lexmovenext(lexpos);

        commonarraysyntax(lexpos, lineno, array2no, nindex2, indexarray2);

        if (mult==outArray[lexpos]||
            add==outArray[lexpos]||
            sub==outArray[lexpos]||
            divide==outArray[lexpos]||
            tensor==outArray[lexpos])
        {
            op=outArray[lexpos];
            if (op == divide) op = divide_array;
            lexmovenext(lexpos);

            commonarraysyntax(lexpos, lineno, array3no, nindex3, indexarray3);

            Instruct.insert(op, array1no, nindex1, indexarray1,
                array2no, nindex2, indexarray2,
                array3no, nindex3, indexarray3);
        }
        else /* = */
        {
            Instruct.insert(assign, array1no, nindex1, indexarray1,
                array2no, nindex2, indexarray2);
        }

    }
    else
        outerror(lineno,  "= is missing.");

    moveoverreturn(lexpos, lineno);
}



void setindexsyntax(int&lexpos, int &lineno)
{
    int arrayno, indexarray[mx_array_dim], nindex=0;
    int i;

    lexmovenext(lexpos);                                   //setindex

    if (id!=outArray[lexpos])
        outerror(lineno, "id is missing.");
    else
    {
        arrayno=idnosyn++;
        lexmovenext(lexpos);
    }

    for (i = 0; i < mx_array_dim; i++)
       indexarray[i] = 0;

    while (id==outArray[lexpos]&& !endofoutArray)
    {
        if (nindex<mx_array_dim)
        {
            indexarray[nindex]=idnosyn;
        }
        nindex++;
        idnosyn++;
        lexmovenext(lexpos);
    }
    if(nindex>mx_array_dim)
        outerror(lineno,
        "resetindex: Too many indices detected.");

/*    Instruct.insert(setindex, arrayno, nindex, indexarray);*/
    outerror(lineno, "resetindex not supported anymore.");

    moveoverreturn(lexpos, lineno);

}



void procsyntax(int &lexpos, int &lineno)
{
    int variableno=0;
    returncount=0;
    lexmovenext(lexpos);                                   //proc
    if (id!=outArray[lexpos])
        outerror(lineno, "id is missing.");
    else
    {
        variableno=idnosyn;
        Instruct.insert(proc, variableno);
        idnosyn++;
        lexmovenext(lexpos);
    }

    moveoverreturn(lexpos, lineno);

    while(!endofoutArray&&(IsInstruction(lexpos) ))
    {
            instructionsyntax(lexpos, lineno);
    }
/*
	while(!endofoutArray&&(IsInstruction(lexpos) ))//||returnsil==outArray[lexpos]))
    {
        while(IsInstruction(lexpos) && !endofoutArray)
        {
            instructionsyntax(lexpos, lineno);
        }
        if (returnsil==outArray[lexpos])
        {
            Instruct.insert (endproc, variableno);
            lexmovenext(lexpos);
            moveoverreturn(lexpos, lineno);

        }
    }
*/
    if (endproc!=outArray[lexpos])
        outerror(lineno,  "endproc is missing.");
    else
        lexmovenext(lexpos);

    if (idname[variableno]!=idname[idnosyn])
        outerror(lineno, "proc name does not match.");
    else
    {
        Instruct.insert (endproc, variableno);
        idnosyn++;
        lexmovenext(lexpos);
    }
    returncount=-1;
    moveoverreturn(lexpos, lineno);
}



void returnsilsyntax(int &lexpos, int &lineno)
{
    int variableno=0;
    if (returncount==-1)
        outerror(lineno, "return cannot be used in the main part of sial.");

    if (returnsil==outArray[lexpos])
    {
        Instruct.insert (endproc, variableno);
        lexmovenext(lexpos);
        moveoverreturn(lexpos, lineno);

    }
    else
        outerror(lineno, "unknown error in returnsilsyntax.");
}



void readwritesyntax(int &lexpos, int &lineno)
{
    lexpos=0;
    outerror(lineno, "read / write instruction has not implemented.");
}



/*====================================
The function instructionsyntax anylizes the instruction in function.
There are two parts: assign and semicolon. Assign will call assignsyntax function.
====================================*/
void instructionsyntax(int &lexpos, int &lineno)
{
    glno=lineno;
    if (endofoutArray) return;
    if (setindex==outArray[lexpos])
        setindexsyntax(lexpos, lineno);

    else if (put==outArray[lexpos])
        putpreparesyntax(lexpos, lineno);

    else if (collective==outArray[lexpos])
        collectivesyntax(lexpos, lineno);

    else if (execute==outArray[lexpos])
        executesyntax(lexpos, lineno);

    else if (id==outArray[lexpos]&&lexpos+1<outarray_size&&
        (lparen==outArray[lexpos+1]||assign==outArray[lexpos+1]||
            outArray[lexpos+1]==addassign||
            outArray[lexpos+1]==subassign||
            outArray[lexpos+1]==multassign||
            outArray[lexpos+1]==divassign))
        operationidsyntax(lexpos, lineno);

    else if (doo==outArray[lexpos])
        dosyntax(lexpos, lineno);

    else if (call==outArray[lexpos])
        callsyntax(lexpos, lineno);

    else if (create==outArray[lexpos])
        createsyntax(lexpos, lineno);

    else if (deletee==outArray[lexpos])
        deletesyntax(lexpos, lineno);

    else if (iff==outArray[lexpos])
        ifsyntax(lexpos, lineno);

    else if (pardo==outArray[lexpos])
        pardosyntax(lexpos, lineno);

    else if (where==outArray[lexpos])
        wheresyntax(lexpos, lineno);

    else if (get==outArray[lexpos])
        getcomputesyntax(lexpos, lineno);

    else if (prepare==outArray[lexpos])
        putpreparesyntax(lexpos, lineno);

    else if (request==outArray[lexpos])
        requestsyntax(lexpos, lineno);

    else if (computeint==outArray[lexpos])
        getcomputesyntax(lexpos, lineno);

    else if (destroye==outArray[lexpos])
        destroysyntax(lexpos, lineno);

    else if (writee==outArray[lexpos])
        readwritesyntax(lexpos, lineno);

    else if (returnsil==outArray[lexpos])
        returnsilsyntax(lexpos, lineno);

    else if (exitt==outArray[lexpos]||cycle==outArray[lexpos])
        exitcyclesyntax(lexpos, lineno);

    else if (allocatee==outArray[lexpos]||deallocatee==outArray[lexpos])
        allocatesyntax(lexpos, lineno);

    else if (outArray[lexpos] == prequest)
        prequest_syntax(lexpos, lineno);
    else
        outerror(lineno, "error in instruction syntax");
}



/*====================================
The function functionsyntax anylyzes the token array stored in outArray.
variable lexpos keeps tracking the position of the outArray.
string temp is the char array of lexical meanings.
====================================*/
void functionsyntax(int &lexpos, int& lineno)
{
    int sialnameno=0;
    if (outArray[lexpos]!=sial)
        outerror(lineno, "Data type is missing.");
    else
    {
        if(!lexmovenext(lexpos)) return;
    }

    if (id!=outArray[lexpos])               //id
        outerror(lineno, "<---Identifier missing.");
    else
    {
        sialnameno=idnosyn++;
        if(!lexmovenext(lexpos)) return;
    }
    moveoverreturn(lexpos, lineno);


    while (IsDeclaration (lexpos) && !endofoutArray)        //declaration
    {
        declarationsyntax(lexpos, lineno);
    }

    while (IsProc (lexpos)&& !endofoutArray)                //proc
    {
        procsyntax(lexpos, lineno);
    }

    Instruct.insert(start, 0);

    while(!endofoutArray &&endsial!=outArray[lexpos])
    {
        while(IsInstruction (lexpos) && !endofoutArray )    //instruction
        {
            instructionsyntax(lexpos, lineno);
        }
        if (!endofoutArray&&endsial!=outArray[lexpos])
            moveoverreturn(lexpos, lineno);
    }

    if (endsial!=outArray[lexpos])                         //endsial
    {
        outerror(lineno, "not instructions. endsial is missing.");
        return;
    }
    else
        lexmovenext(lexpos);

    if (idname[sialnameno]!=idname[idnosyn])
    {
        outerror(lineno, "endsial name does not match.");
    }
    else
    {
        idnosyn++;
        lexmovenext(lexpos);
    }
}
