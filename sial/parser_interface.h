#include <stdio.h>
#include "f77_name.h"
#include "f_types.h"
#ifdef __cplusplus /* Tell the compiler not mangle the name of foo1*/
extern "C"{
#endif


f_int F77_NAME(create_optable, CREATE_OPTABLE) (f_int *nentries);


f_int F77_NAME(create_index_table, CREATE_INDEX_TABLE) (f_int *nentries);


f_int F77_NAME(create_array_table, CREATE_ARRAY_TABLE) (f_int *nentries);


f_int F77_NAME(add_optable, ADD_OPTABLE) (f_int *opcode, 
                                          f_int *op1_array,
                                          f_int *op2_array,
                                          f_int *result_array,
                                          f_int *indarray,
                                          f_int *user_sub_index,
                                          f_int *lineno);



f_int F77_NAME(add_index_table, ADD_INDEX_TABLE) (f_int *index_type,
                                                  f_int *beginning_seg,
                                                  f_int *ending_seg);


f_int F77_NAME(add_array_table, ADD_ARRAY_TABLE) (f_int *nindex,
                                                  f_int *array_type,
                                                  f_int *numblks,
                                                  f_int *indarray,
                                                  double *scalar_value);

void F77_NAME(dump_sip_tables, DUMP_SIP_TABLES)();

void F77_NAME(write_tables, WRITE_TABLES) (char *obj_filename);

int c_get_subroutine_handle(char *sub_name);
#ifdef __cplusplus
}
#endif
