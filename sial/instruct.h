#ifndef _instruct_h
#define _instruct_h
const int maxinstructionnumber=1000000;
const int in_index_mask = 65536;

class InstructCls
{
    int numberofins;
    int stackno;
	int stacknof;
    int **instruction;
    int CheckConvert(int &idno, int nindex, int indexarray[]);
    int CheckConvertAllocate(int &idno, int nindex, int indexarray[]);
    int HaveEnoughSize();

public:
    InstructCls();
    ~InstructCls();

    void insert(int type, int idno, int nindex1, int indexarray1[],
                         int idno2, int nindex2, int indexarray2[],
                         int idno3, int nindex3, int indexarray3[]);
    void insert(int type, int idno, int nindex, int*indexarray,
                          int idno1, int nindex1, int* indexarray1);
    void insert(int inscode, int isno, 
                int idno, int nindex, int indexarray[]);

    void insert(int type, int idno, int nindex, int*indexarray);
    void insertallocate(int type, int idno, int nindex, int*indexarray);
    void insert(int iscode, int nindex, int*indexarray);
    void insert(int iscode, int nindex, int*indexarray, int flagarray[]);
    void insert(int inscode, int isno, int idno, int idno1);
    void insert(int inscode, int isno, int idno);
    void insert(int inscode, int idno);
    int inserts(int inscode, int result, int op1, int op2);
    int insertsf(int inscode, double result, int op1, int op2);
    int insertif(int inscode, int address);
    int reinsert(int structno, int address);
    int getinsno() const;


    void output() const;
    void outputll() const;
    void initializer();
};
#endif
