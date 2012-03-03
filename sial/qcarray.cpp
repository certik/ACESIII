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
qcarray.cpp
Lei Wang
wang@qtp.ufl.edu
May 2004
******************************/
#include <cstdio>
#include <string>
#include "qcarray.h"
#include "variable.h"
/*#include "parsetest.h"*/
#include "parser_interface.h"
#include "keywdcnt.h"
#include "errorhl.h"
#include "instruct.h"
//using namespace std;

extern int mx_array_dim;
extern int glno;
extern string idname[];
extern int idno;
extern VariableClass Variable;
extern InstructCls Instruct;

QCArrayClass::QCArrayClass()
{
    nvars=22;
}

void QCArrayClass::initializer()
{
    vararray    = new string[maxarraynumber];
    arrayindex  = new int*[maxarraynumber];
    for (int i = 0; i < maxarraynumber; i++)
       arrayindex[i] = new int[mx_array_dim];
    arraynindex = new int[maxarraynumber];
    arraytype   = new int[maxarraynumber];
    fvalue      = new double[maxarraynumber];

    insertconstant();
}



QCArrayClass::~QCArrayClass()
{
    nvars=0;
    if (vararray) delete [] vararray;
    if (arrayindex) delete [] arrayindex;
    if (arraynindex) delete [] arraynindex;
    if (arraytype) delete [] arraytype;
    if (fvalue) delete [] fvalue;
}


void QCArrayClass::insertconstant()
{
    vararray[0]="c";
    arrayindex[0][0]=3;
    arrayindex[0][1]=0;
    arraynindex[0]=2;
    arraytype[0]=staticc;

    vararray[1]="ca";
    arrayindex[1][0]=3;
    arrayindex[1][1]=1;
    arraynindex[1]=2;
    arraytype[1]=staticc;

    vararray[2]="cb";
    arrayindex[2][0]=3;
    arrayindex[2][1]=2;
    arraynindex[2]=2;
    arraytype[2]=staticc;

    vararray[3]="scfeneg";
    arraytype[3]=scalar;
    fvalue[3]=0;

    vararray[4]="totenerg";
    arraytype[4]=scalar;
    fvalue[4]=0;

    vararray[5]="e";
    arrayindex[5][0]=3;
    arrayindex[5][1]=0;
    arraynindex[5]=2;
    arraytype[5]=staticc;

    vararray[6]="ea";
    arrayindex[6][0]=3;
    arrayindex[6][1]=1;
    arraynindex[6]=2;
    arraytype[6]=staticc;

    vararray[7]="eb";
    arrayindex[7][0]=3;
    arrayindex[7][1]=2;
    arraynindex[7]=2;
    arraytype[7]=staticc;

    vararray[8] = "fock_a";
    arrayindex[8][0]=1;
    arrayindex[8][1] = 1;
    arraynindex[8] = 2;
    arraytype[8] = staticc;

    vararray[9] = "fock_b";
    arrayindex[9][0]= 2;
    arrayindex[9][1] = 2;
    arraynindex[9] = 2;
    arraytype[9] = staticc;

    vararray[10]="oed_nai";
    arrayindex[10][0]=3;
    arrayindex[10][1]=3;
    arraynindex[10]=2;
    arraytype[10]=staticc;

    vararray[11]="oed_kin";
    arrayindex[11][0]=3;
    arrayindex[11][1]=3;
    arraynindex[11]=2;
    arraytype[11]=staticc;

    vararray[12]="oed_ovl";
    arrayindex[12][0]=3;
    arrayindex[12][1]=3;
    arraynindex[12]=2;
    arraytype[12]=staticc;

    vararray[13]="damp";
    arraytype[13]=scalar;
    fvalue[13]=0;

    vararray[14]="cc_conv";
    arraytype[14]=scalar;
    fvalue[14]=0;

    vararray[15]="scf_conv";
    arraytype[15]=scalar;
    fvalue[15]=0;

    vararray[16]="fockrohf_a";
    arrayindex[16][0]=3;
    arrayindex[16][1]=3;
    arraynindex[16]=2;
    arraytype[16]=staticc;

    vararray[17]="fockrohf_b";
    arrayindex[17][0]=3;
    arrayindex[17][1]=3;
    arraynindex[17]=2;
    arraytype[17]=staticc;

    vararray[18]="stabvalue";
    arraytype[18]=scalar;
    fvalue[18]=0;

    vararray[19]="eom_tol";
    arraytype[19]=scalar;
    fvalue[19]=0;

    vararray[20]="eom_roots";
    arraytype[20]=scalar;
    fvalue[20]=0;

    vararray[21]="excite";
    arraytype[21]=scalar;
    fvalue[21]=0;

    nvars = 22;
}



void QCArrayClass::insert(int type, int idno, int nindex, int*indexarray)
{
    int i, j;
    int variableid;
    int found;

    if (Variable.getidnumber(idno)!=-1)
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

    i=0;

    if (nvars<maxidnumber)
    {
        found=1;
        while(i<nindex&&found)
        {
            variableid=Variable.getidnumber(indexarray[i]);

            if(variableid==-1)
                found=0;
            else
            {
                arrayindex[nvars][i]=variableid;
                i++;
            }
        }

        if (!found)
        {
            outerror(glno, "index not defined.");
            printf("%s, \n",idname[indexarray[i]].c_str());
            return;
        }

        /*add 0 to meet the ll needs*/
        for (j=i; j<mx_array_dim; j++)
            arrayindex[nvars][j]=0;

        vararray[nvars]=idname[idno];
        arraynindex[nvars]=nindex;
        arraytype[nvars]=type;
        nvars++;
    }
    else
    {
        outerror(glno, 
            "More id than id container can hold, please contact writer.\n");
    }
}



void QCArrayClass::insertscalar(int idno, double value)
{
    if (nvars>=maxidnumber)
    {
        outerror(glno, 
            "More id than id container can hold, please contact writer.\n");
        return;
    }
    if (Variable.getidnumber(idno)!=-1||getidnumber(idno)!=-1)
    {
        outerror(glno, "id exist.");
        return;
    }

    vararray[nvars]=idname[idno];
    arraytype[nvars]=scalar;
    arraynindex[nvars]=0;
    arrayindex[nvars][0]=0;
    fvalue[nvars]=value;
    nvars++;
}


int QCArrayClass::inserttempscalar(double value)
{
    int i;

    if (nvars>=maxidnumber)
    {
        outerror(glno, 
            "More id than id container can hold, please contact writer.\n");
        return -1;
    }

    for (i = 0; i < nvars; i++)
    {
       if (vararray[i] == "temp" && fvalue[i] == value) 
          return i;
    }

    vararray[nvars]="temp";
    arraytype[nvars]=scalar;
    arraynindex[nvars]=0;
    arrayindex[nvars][0]=temp;
    fvalue[nvars]=value;
    return nvars++;
}




/*
int QCArrayClass::getidnumber(string idname)
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



int QCArrayClass::getidnumber(int idnumber) const
{
    int i;
    for (i=0; i<nvars&&idname[idnumber]!=vararray[i]; i++);
    if (i<nvars&&idname[idnumber]==vararray[i])
        return i;
    else
        return -1;
}



int QCArrayClass::gettype(int idnumber) const
{
    if (idnumber<nvars)
        return arraytype[idnumber];
    else
        return 0;
}



string QCArrayClass::getidname(int idnumber) const
{
    if (idnumber<nvars)
        return vararray[idnumber];
    else
    {
        printf("Error in getidname: idnumber %d is larger than nvars %d\n",
               idnumber, nvars);
        return 0;
    }
}



int QCArrayClass::istempscalar(int idnumber) const
{
    if (idnumber<nvars)
        return (arraytype[idnumber]==scalar&&
            arraynindex[idnumber]==0&&arrayindex[idnumber][0]==temp);
    else
        return -1;
}


int QCArrayClass::isscalar(int idnumber) const
{
    if (idnumber<nvars)
        return (arraytype[idnumber]==scalar);
    else
        return -1;
}



/*
string QCArrayClass::getidname(int idnumber) const
{
    if (idnumber<nvars)
        return vararray[idnumber];
    else
        return 0;
}
*/

/*
int QCArrayClass::CheckExeIndex(int idno, int nindex, int indexarray[]) const
{
    int i;
    if (nindex!=arraynindex[idno])
    {
        outerror(glno, "CheckExeIndex: nindex doesn't match.");
        return -3;
    }
    for (i=0; i<nindex; i++)
    {
        if (Variable.getidnumber(indexarray[i])!=arrayindex[idno][i])
        {
            outerror(glno, "index definition does not match before +=.");
            return -2;
        }
    }
    return 0;
}



int QCArrayClass::CheckAccIndex(int array1, int array2,
                      int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[]) const
{
    int i;
    if (arraynindex[array1]!=arraynindex[array2])
    {
        outerror(glno, "nindex does not match.");
        return -3;
    }

    for (i=0; i<arraynindex[array1]; i++)
    {
        if (Variable.getidnumber(indexarray1[i])!=arrayindex[array1][i])
        {
            outerror(glno, "index definition does not match before +=.");
            return -2;
        }
        if (Variable.getidnumber(indexarray2[i])!=arrayindex[array2][i])
        {
            outerror(glno, "index definition does not match after +=.");
            return -2;
        }

        if (arrayindex[array1][i]!=arrayindex[array2][i])
        {
            outerror(glno, "CheckAccIndex: index does not match.");
            return -1;
        }
    }

    return 0;
}
*/

int IsParameter(int x)
{
   int i, y;

   i = 0;
   if (x < 0)
   {
      y = -x;
      if (y == cc_iter ||
          y == cc_hist ||
          y == cc_beg  ||
          y == scf_iter ||
          y == scf_hist ||
          y == scf_beg  ||
          y == natoms   ||
          y == itrips ||
          y == itripe ||
          y == ihess1 ||
          y == ihess2 ||
          y == jhess1 ||
          y == jhess2 ||
          y == subb   ||
          y == sube   ||
          y == sip_sub_segsize ||
          y == sip_sub_occ_segsize ||
          y == sip_sub_virt_segsize ||
          y == sip_sub_ao_segsize) i = 1;
   }
   return i;
}


int QCArrayClass::CheckArray(int &arrayid, int nindex, int indexarray[]) const
{
    int i;
    int number;
    int indextype, number_type;
    int b1, b2, e1, e2;
    int rtnvalue=0;
    int superindex, subindex_flag; 

    if (arrayid<=0)
        arrayid=-arrayid;
    else
        arrayid=getidnumber(arrayid);

    if (arrayid==-1)
    {
        outerror(glno, "array or scalar does not exist.");
        return -2;
    }

    if (nindex!=arraynindex[arrayid])
    {
        outerror(glno, "CheckArray: nindex doesn't match.");
        return -2;
    }
    for (i=0; i<nindex; i++)
    {
        number=Variable.getidnumber(indexarray[i]);
        number_type = Variable.getidtype(number);
        indexarray[i]=number;

        if (number==-1)
        {
            outerror(glno, "index not defined.");
            return -1;
        }

        if (number_type == subindex)
        {
           //  Look up the type of the superindex instead of the subindex.

           superindex = Variable.getrange1(indexarray[i]);
           number_type = Variable.getidtype(superindex);
           number      = superindex;
        } 

        indextype = Variable.getidtype(arrayindex[arrayid][i]);
        if (indextype == subindex)
        {
                      //  Look up the type of the superindex instead of the subindex.

           superindex = Variable.getrange1(arrayindex[arrayid][i]);
           indextype  = Variable.getidtype(superindex);
           subindex_flag = 1; 
        }
        else
        {
           subindex_flag = 0;  
        }

           if (number_type!= indextype)
           {
               printf("Variable.getidtype(number) = %d\n",Variable.getidtype(number));
               printf("Variable.getidtype(arrayindex[arrayid][i]) = %d\n",
                       Variable.getidtype(arrayindex[arrayid][i]));
               outerror(glno, "!!! index types do not match.");
               return -1;
           }
/*continue;*/
           b1=Variable.getbegrange(number);
           e1=Variable.getendrange(number);
           
           if (subindex_flag) 
           {
              b2 = Variable.getbegrange(superindex);
              e2 = Variable.getendrange(superindex);
           }
           else
           {  
              b2=Variable.getbegrange(arrayindex[arrayid][i]);
              e2=Variable.getendrange(arrayindex[arrayid][i]);
           }

           if (!IsParameter(b2) && !IsParameter(e2))
           {
              if (b1!=b2)
              {
                  if ((b2==-bocc &&b1==-bvirt)||
                      (b2==-baocc&&b1==-bavirt)||
                      (b2==-bbocc&&b1==-bbvirt)||
                      b2==1)
                  {
                      /*should be ok, but warning*/
                  }
                  else
                  {
                      outerror(glno, "QCArrayClass::CheckArray #1");
                      outerror(glno, "index beginning ranges do not match.");
                      printf("  In array '%s', index %d, using variable '%s', beginning value: ",
                      getidname(arrayid).c_str(), i, Variable.getidname(number).c_str());
                      if (b1<=0)
                          printf("%s; ", keywords[-b1]);
                      else
                          printf("%1d; ", b1);
                      printf("\n  By definition: variable '%s', beginning value: ",
                      Variable.getidname(arrayindex[arrayid][i]).c_str());
                      if (b2<=0)
                          printf("%s.", keywords[-b2]);
                      else
                          printf("%1d.", b2);
                      printf("\n");
                      rtnvalue=-1;
                      continue;
                  }
              }   /* b1 != b2 */
              if (e1!=e2)
              {
                  if ((e1==-eocc &&e2==- evirt)||
                      (e1==-eaocc&&e2==-eavirt)||
                      (e1==-ebocc&&e2==-ebvirt)||
                      (e1==-eocc &&e2==-norb)||
                      (e1==-eaocc&&e2==-norb)||
                      (e1==-ebocc&&e2==-norb)||
                      (e1==-evirt&&e2==-norb)||
                      (e1==-eavirt&&e2==-norb)||
                      (e1==-ebvirt&&e2==-norb)||
                      /*question about it, but should be ok*/
                      (e1==e2)
                     )
                  {
                      /*should be ok, but warning*/
                  }
                  else
                  {
                      outerror(glno, "QCArray::CheckArray #2");
                      outerror(glno, "index ending range does not match.");
                      printf("  In array '%s', index %d, using variable '%s', ending value: ",
                       getidname(arrayid).c_str(), i, Variable.getidname(number).c_str());
                      if (e1<=0)
                          printf("%s; ", keywords[-e1]);
                      else
                          printf("%1d; ", e1);
                      printf("\n  By definition: variable '%s', ending value: ",
                       Variable.getidname(arrayindex[arrayid][i]).c_str());
                      if (e2<=0)
                          printf("%s.", keywords[-e2]);
                      else
                          printf("%1d.", e2);
                      printf("\n");
                      rtnvalue=-1;
                      continue;
      /*                printf("number=%d, arrayid=%d, i=%d, ",number, arrayid, i);
                      printf("arrayindex[arrayid][i]=%d, %d, %d\n", 
                          arrayindex[arrayid][i],
                                  e1, e2);
                      outerror(glno, "index ranges do not match.");
                      return -1;
                      */
                  }
              }   /* e1 != e2 */
           }      /* !IsParameter... */
    }
    return rtnvalue;
}



int QCArrayClass::CheckArrayAllocate(int &arrayid, int nindex, int indexarray[]) const
{
    int i;
    int number;
    int indextype, number_type;
    int b1, b2, e1, e2;
    int rtnvalue=0;

    if (arrayid<=0)
        arrayid=-arrayid;
    else
        arrayid=getidnumber(arrayid);

    if (arrayid==-1)
    {
        outerror(glno, "array or scalar does not exist.");
        return -2;
    }

    if (nindex!=arraynindex[arrayid])
    {
        outerror(glno, "CheckArrayAllocate: nindex doesn't match.");
        return -2;
    }
    for (i=0; i<nindex; i++)
    {
        if (indexarray[i]!=WILDCARD_INDICATOR)  /* * */
        {
            number=Variable.getidnumber(indexarray[i]);
            number_type = Variable.getidtype(number); 
            indexarray[i] = number+1;  // Return (Fortran) index rather than id.

            if (number==-1)
            {
                outerror(glno, "index not defined.");
                return -1;
            }

               indextype = Variable.getidtype(arrayindex[arrayid][i]);
                  if (number_type != indextype)
                  {
                      printf("Variable.getidtype(number) = %d\n",
                           Variable.getidtype(number));
                       printf("Variable.getidtype(arrayindex[arrayid][i]) = %d\n",
                               Variable.getidtype(arrayindex[arrayid][i])); 
                      outerror(glno, "index types do not match.");
                      return -1;
                  }

    /*continue;*/

               b1=Variable.getbegrange(number);
               e1=Variable.getendrange(number);
               b2=Variable.getbegrange(arrayindex[arrayid][i]);
               e2=Variable.getendrange(arrayindex[arrayid][i]);

               if (!IsParameter(b2) && !IsParameter(e2))
               {
                  if (b1!=b2)
                  {
                      if ((b2==-bocc &&b1==-bvirt)||
                          (b2==-baocc&&b1==-bavirt)||
                          (b2==-bbocc&&b1==-bbvirt)||
                          b2==1)
                      {
                          /*should be ok, but warning*/
                      }
                      else
                      {
                          outerror(glno, "QCArrayClass::CheckArrayAllocate");
                          outerror(glno, "index beginning ranges do not match.");
                          printf("  In array '%s', index %d, using variable '%s', beginning value: ",
                          getidname(arrayid).c_str(), i, Variable.getidname(number).c_str());
                          if (b1<=0)
                              printf("%s; ", keywords[-b1]);
                          else
                              printf("%1d; ", b1);
                          printf("\n  By definition: variable '%s', beginning value: ",
                          Variable.getidname(arrayindex[arrayid][i]).c_str());
                          if (b2<=0)
                              printf("%s.", keywords[-b2]);
                          else
                              printf("%1d.", b2);
                          printf("\n");
                          rtnvalue=-1;
                          continue;
                      }
                  }    /* b1 != b2 */
                  if (e1!=e2)
                  {
                      if ((e1==-eocc &&e2==- evirt)||
                          (e1==-eaocc&&e2==-eavirt)||
                          (e1==-ebocc&&e2==-ebvirt)||
                          (e1==-eocc &&e2==-norb)||
                          (e1==-eaocc&&e2==-norb)||
                          (e1==-ebocc&&e2==-norb)||
                          (e1==-evirt&&e2==-norb)||
                          (e1==-eavirt&&e2==-norb)||
                          (e1==-ebvirt&&e2==-norb)||
                          /*question about it, but should be ok*/
                          (e1==e2)
                      )
                      {
                          /*should be ok, but warning*/
                      }
                      else
                      {
                          outerror(glno, "index ending ranges do not match.");
                          printf("  In array '%s', index %d, using variable '%s', ending value: ",
                              getidname(arrayid).c_str(), i, Variable.getidname(number).c_str());
                          if (e1<=0)
                              printf("%s; ", keywords[-e1]);
                          else
                              printf("%1d; ", e1);
                          printf("\n  By definition: variable '%s', ending value: ",
                           Variable.getidname(arrayindex[arrayid][i]).c_str());
                          if (e2<=0)
                              printf("%s.", keywords[-e2]);
                          else
                              printf("%1d.", e2);
                          printf("\n");
                          rtnvalue=-1;
                          continue;
          /*                printf("number=%d, arrayid=%d, i=%d, ",number, arrayid, i);
                          printf("arrayindex[arrayid][i]=%d, %d, %d\n", 
                              arrayindex[arrayid][i],
                                      e1, e2);
                          outerror(glno, "index ranges do not match.");
                          return -1;
                          */
                      }
		  }   /* e1 != e2 */
	       }      /* !IsParameter */
        }  /* else indexarray[i] != WILDCARD_INDICATOR */
    } /*end of for*/
    return rtnvalue;
}




int QCArrayClass::CheckSumDiff( int nindex1, int indexarray1[],
                                int nindex2, int indexarray2[],
                                int nindex3, int indexarray3[]) const
{
    int i;

    if (nindex1!=nindex2 || nindex1 != nindex3)
    {
        outerror(glno, "nindex does not match.");
        return -3;
    }

    for (i=0; i<nindex1; i++)
    {
        if (indexarray1[i] != indexarray2[i] ||
            indexarray1[i] != indexarray3[i])
        {
            outerror(glno, "CheckSumDiff: index does not match.");
            return -1;
        }
    }
    return 0;
}



void selectionSort(const int oldarray[], int numbers[], int array_size)
{
  int i, j;
  int min, temp;

    for (i=0; i<array_size; i++)
        numbers[i]=oldarray[i];

    for (i = 0; i < array_size-1; i++)
    {
        min = i;
        for (j = i+1; j < array_size; j++)
        {
            if (numbers[j] < numbers[min])
                min = j;
        }

        temp = numbers[i];
        numbers[i] = numbers[min];
        numbers[min] = temp;
    }
}


/*
int QCArrayClass::resetindex(int idno, int nindex, const int indexarray[])
{
    int i;
    if (idno>=nvars)
    {
        outerror(glno, "error in resetindex, please contact writer.");
        return -1;
    }
    if (nindex!=arraynindex[idno])
    {
        outerror(glno, "QCArrayClass: nindex doesn't match.");
        return -2;
    }
    
    for (i=0; i<nindex; i++)
    {
        arrayindex[idno][i]=indexarray[i];
    }
    return 0;
}
*/

/*
void QCArrayClass::implicitsetindex(int idno, int nindex, int * indexarray)
{
    Instruct.insert(setindex, idno, nindex, indexarray);
}
*/

/*
int QCArrayClass::CheckAssignIndex(int array1, int array2,
                      int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[]) const
{
    int i;
    int variableindexarray1[mx_array_dim], variableindexarray2[mx_array_dim];
    int sortstorearray1[mx_array_dim], sortstorearray2[mx_array_dim];
    int var_index_type, id_array1, id_array2;

    if (nindex1!=arraynindex[array1] ||
        nindex2!=arraynindex[array2])
    {
        outerror(glno, "array dimension doesn't match.");
        return -4;
    }


    if (nindex1 == nindex2)
    {
       for (i=0; i<arraynindex[array1]; i++)
       {
           variableindexarray1[i]=Variable.getidnumber(indexarray1[i]);
           variableindexarray2[i]=Variable.getidnumber(indexarray2[i]);
           if (variableindexarray1[i]!=arrayindex[array1][i])
           {
               outerror(glno, "index definition does not match before =.");
               return -2;
           }
           if (Variable.getidnumber(indexarray2[i])!=arrayindex[array2][i])
           {
               outerror(glno, "index definition does not match after =.");
               return -2;
           }
       }

       selectionSort(arrayindex[array1], sortstorearray1, nindex1);

       selectionSort(arrayindex[array2], sortstorearray2, nindex2);
       for (i=0; i<nindex1; i++)
       {
           if (sortstorearray1[i]!=sortstorearray1[i])
           {
               outerror(glno, "two array indices do not match.");
               return -1;
           }
       }
    }
    
    if (nindex1 < nindex2) 
    {
       for (i = 0; i < nindex1; i++)
       {
          if (arrayindex[array1][i] != arrayindex[array2][i])
          {
             outerror(glno, "Array indices do not match");
             return -1;
          }
       }  

       id_array2 = Variable.getidnumber(array2);

       for (i = nindex1; i < nindex2; i++)
       {
          var_index_type = Variable.getidtype(id_array2);
          if (var_index_type != indexx) 
          {
             outerror(glno, "Final indices in assign must all be type index");
             return -1;
          }
       }
    }

    if (nindex1 > nindex2)
    {
       for i = 0; i < nindex2; i++)
       {
          if (arrayindex[array1][i] != arrayindex[array2][i])
          {
             outerror(glno, "Array indices do not match");
             return -1;
          }
       }

       id_array1 = Variable.getidnumber(array1);

       for (i = nindex2; i < nindex1; i++)
       {
          var_index_type = Variable.getidtype(id_array1);
          if (var_index_type != indexx)
          {
             outerror(glno, "Final indices in assign must all be type index");
             return -1;
          }
       }
    }

    return 0;
}
*/

int QCArrayClass::CheckMult(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[],
                      int nindex3, int indexarray3[],
                      int &nindexr, int indexarrayr[]) const
{
    int i, i2, i3;
    nindexr=0;
    int sortindexarray1[mx_array_dim], sortindexarray2[mx_array_dim], 
        sortindexarray3[mx_array_dim];
    int contractionarray[2*mx_array_dim], contractnindex;

    selectionSort(indexarray1, sortindexarray1, nindex1);
    selectionSort(indexarray2, sortindexarray2, nindex2);
    selectionSort(indexarray3, sortindexarray3, nindex3);

/*check contraction index of array*/
    /*merge two arrays*/
    i2=i3=0;
    contractnindex=0;
    while (i2<nindex2 && i3<nindex3)
    {
        if (sortindexarray2[i2]>sortindexarray3[i3])
            contractionarray[contractnindex++]=sortindexarray3[i3++];
        else if (sortindexarray2[i2]<sortindexarray3[i3])
            contractionarray[contractnindex++]=sortindexarray2[i2++];
        else /*==*/
        {
            indexarrayr[nindexr++]=sortindexarray2[i2];
            i2++;
            i3++;
        }
    }
    while (i2<nindex2)
            contractionarray[contractnindex++]=sortindexarray2[i2++];
    while (i3<nindex3)
            contractionarray[contractnindex++]=sortindexarray3[i3++];

    /*compare array1 and calculated contraction array*/
    if (nindex1!=contractnindex)
    {
        outerror(glno, "error in contraction index.");
        return -2;
    }
    
    for (i=0; i<nindex1; i++)
    {
        if (contractionarray[i]!=sortindexarray1[i])
        {
            outerror(glno, "contraction array doesn't match.");
            return -1;
        }
    }
    return 0;
}



int QCArrayClass::CheckTensor(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[],
                      int nindex3, int indexarray3[],
                      int &nindexr, int indexarrayr[]) const
{
    int i, i2, i3;
    nindexr=0;
    int sortindexarray1[mx_array_dim], sortindexarray2[mx_array_dim], 
        sortindexarray3[mx_array_dim];
    int contractionarray[2*mx_array_dim], contractnindex;

    selectionSort(indexarray1, sortindexarray1, nindex1);
    selectionSort(indexarray2, sortindexarray2, nindex2);
    selectionSort(indexarray3, sortindexarray3, nindex3);

/*check tensor index of array*/
    /*merge two arrays*/
    i2=i3=0;
    contractnindex=0;
    while (i2<nindex2 && i3<nindex3)
    {
        if (sortindexarray2[i2]>sortindexarray3[i3])
            contractionarray[contractnindex++]=sortindexarray3[i3++];
        else if (sortindexarray2[i2]<sortindexarray3[i3])
            contractionarray[contractnindex++]=sortindexarray2[i2++];
        else /*==*/
        {
            indexarrayr[nindexr++]=sortindexarray2[i2];
            i2++;
            i3++;
        }
    }
    while (i2<nindex2)
            contractionarray[contractnindex++]=sortindexarray2[i2++];
    while (i3<nindex3)
            contractionarray[contractnindex++]=sortindexarray3[i3++];

    /*compare array1 and calculated contraction array*/
    if (nindex1!=contractnindex)
    {
        outerror(glno, "error in tensor index.");
        return -2;
    }
    
    for (i=0; i<nindex1; i++)
    {
        if (contractionarray[i]!=sortindexarray1[i])
        {
            outerror(glno, "tensor array doesn't match.");
            return -1;
        }
    }
    return 0;
}



int QCArrayClass::CheckDivide(int nindex1, int indexarray1[],
                    int nindex2, int indexarray2[],
                    int nindex3, int indexarray3[]) const
{
    if (nindex1!=0||nindex2!=0||nindex3!=0)
    {
        outerror(glno, "only scalar is allowed for divide.");
        return -1;
    }
    return 0;
}



int QCArrayClass::Checkcollective(int nindex1, int indexarray1[],
                    int nindex2, int indexarray2[]) const
{
    if (nindex1!=0||nindex2!=0)
    {
        outerror(glno, "only scalar is allowed for Scalar Collective Sum.");
        return -1;
    }
    return 0;
}



int QCArrayClass::CheckAssign(int nindex1, int indexarray1[],
                    int nindex2, int indexarray2[]) const
{
	return 0;
    if (nindex1 != 0 &&nindex2!=0)
    {
        outerror(glno, "Error in array assign, array dimension mismatch.");
        return -1;
    }
    return 0;
}



void QCArrayClass::print()
{
    int i, j;
    for (i=0; i<nvars; i++)
    {
        if (arraytype[i]==scalar)
            printf("%i, %10s, %f", i, vararray[i].c_str(), fvalue[i]);
        else
        {
            printf("%i, %10s, %d, ", i, vararray[i].c_str(), arraynindex[i]);
            for (j=0; j<arraynindex[i]; j++)
                printf("%6s ", Variable.getidname(arrayindex[i][j]).c_str());
        }
        printf("\n");
    }
}


void add_array_table_c(int arraynindex, int type, int isize, int *indarray, 
            double fvalue)
{
    f_int f_arraynindex = arraynindex;
    f_int f_type = type;
    f_int f_isize = isize;
    f_int *findarray;

    int i;

    findarray = new f_int[mx_array_dim];

    for (i = 0; i < mx_array_dim; i++)
    {
       if (i < arraynindex) 
          findarray[i] = indarray[i] + 1;   // Convert to Fortran index
       else
          findarray[i] = 0;
    }

    F77_NAME(add_array_table, ADD_ARRAY_TABLE)
            (&f_arraynindex, &f_type, &f_isize, findarray, 
             &fvalue);

    delete [] findarray;
}



void QCArrayClass::outputll()
{
    int type=0;
    int isize=0;
    int i;
    f_int f_nvars = nvars;

    F77_NAME(create_array_table, CREATE_ARRAY_TABLE)(&f_nvars);
    for(i=0; i<nvars; i++)
    {
        if (served==arraytype[i])
            type=201;
        else if (staticc==arraytype[i])
            type=202;
        else if (distributed==arraytype[i])
            type=203;
        else if (temp==arraytype[i])
            type=204;
        else if (scalar==arraytype[i])
            type=205;
        else if (local==arraytype[i])
            type=206;
        else
            outerror(glno, "error in create_array_table");

        add_array_table_c
            (arraynindex[i], type, isize, &arrayindex[i][0], 
            fvalue[i]);

    }
}
