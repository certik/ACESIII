#ifdef __cplusplus /* Tell the compiler not mangle the name of foo1*/
extern "C"{
#endif

void allocidname();


int interpret (char* infilename, int lext, int indext, int arrayt, int inst);

#ifdef __cplusplus
}
#endif
