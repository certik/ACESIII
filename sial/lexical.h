/******************************
lexical.h
Lei Wang 
wang@qtp.ufl.edu
May 2004
******************************/

#ifndef lexical_h
#define lexical_h

using namespace std;

/*====================================
lexical function lexicalizes the characters stored in sourceArray,
variable sAsize tells the function the size of sourceArray is. 
====================================*/
void lexical(char sourceArray[], int sAsize, 
             int &idno, string idname[], int &rnumber);


#endif
