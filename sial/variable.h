#ifndef _variable_h
#define _variable_h

using namespace std;

const int maxidnumber=10000;

class VariableClass
{
    int nvars;

    string *vararray;
    int *varrange1;
    int *varrange2;
    int *vartype;


public:

    VariableClass();

    ~VariableClass();

    void insertconstant();

    void insert(int type, int idno, int range1, int range2);

    int getidnumber(int idnumber) const;

    string getidname(int idnumber) const;

    int getrange1(int idnumber) const;
    int getidtype(int idnumber) const;
    int getbegrange(int idnumber) const;
    int getendrange(int idnumber) const;
//    int getidnumber(string idname);
    void print() const;
    void outputll() const;
    void initializer();
};
#endif
