
using namespace std;

const int maxarraynumber=10000;

class QCArrayClass
{
    int nvars;

    string *vararray;
    int **arrayindex;
    int *arraynindex;
    int *arraytype;
    double *fvalue;
    void insertconstant();
/*    void implicitsetindex(int idno, int nindex, int * indexarray);
*/

public:

    QCArrayClass();

    ~QCArrayClass();

    void insert(int type, int idno, int nindex, int*indexarray);

    void insertscalar(int idno, double value);

    int inserttempscalar(double value);

    int getidnumber(int idnumber) const;

    int gettype(int idnumber) const;

    string getidname(int idnumber) const;

    int istempscalar(int idnumber) const;

    int isscalar(int idnumber) const;
    /*string getidname(int idnumber) const;*/

    int CheckArray(int &arrayid, int nindex, int indexarray[]) const;
    int CheckArrayAllocate(int &arrayid, int nindex, int indexarray[]) const;
/*
    int resetindex(int idno, int nindex, const int indexarray[]);
    int CheckExeIndex(int array, int nindex, int indexarray[]) const;
    int CheckAccIndex(int array1, int array2,
                      int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[]) const;
    int CheckAssignIndex(int array1, int array2,
                      int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[]) const;
*/
    int CheckSumDiff(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[],
                      int nindex3, int indexarray3[]) const;
    int CheckMult(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[],
                      int nindex3, int indexarray3[],
                      int &nindexr, int indexarrayr[]) const;
    int CheckTensor(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[],
                      int nindex3, int indexarray3[],
                      int &nindexr, int indexarrayr[]) const;
    int CheckDivide(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[],
                      int nindex3, int indexarray3[]) const;
    int Checkcollective(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[]) const;
    int CheckAssign(int nindex1, int indexarray1[],
                      int nindex2, int indexarray2[]) const;
    void print();
    void outputll();
    void initializer();
};
