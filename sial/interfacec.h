#ifndef _interfacec_h
#define _interfacec_h
#include "f_types.h"
#include "f77_name.h"

#ifdef __cplusplus /* Tell the compiler not mangle the name of foo1*/
extern "C"{
#endif


int create_optable_c(int nentries);

#ifdef __cplusplus
}
#endif



#ifdef __cplusplus /* Tell the compiler not mangle the name of foo1*/
extern "C"{
#endif

int add_optable_c(int opcode, 
                    int op1_array,
                    int op2_array,
                    int result_array,
                    int *indarray,
                    int user_sub_index,
                    int lineno);
#ifdef __cplusplus
}
#endif

/*
int create_index_table(int *nentries);

int create_array_table(int *nentries);

int add_index_table (int *index_size,
                                                  int *nsegments,
                                                  int *virtual_index,
                                                  int *equiv_index);

int add_array_table(int *nindex,
                                                  int *array_type,
                                                  int *numblks,
                                                  int *index1,
                                                  int *index2,
                                                  int *index3,
                                                  int *index4);
*/

#endif
