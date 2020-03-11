#ifndef _CCSUBMATLIB_HPP_
#define _CCSUBMATLIB_HPP_

#include "nummatrix.hpp"

class CCSubMatLib {
public:
    class Entry {
    public:
        NumMatrix& SM() { return _SM; }  //subdivision matrix
        NumMatrix& IV() { return _IV; }  //inverse of eigen vectors
        NumMatrix& V()  { return _V;  }  //eigen vectors
        NumVector& D()  { return _D;  }  //eigen values
    protected:
        NumMatrix _SM;
        NumMatrix _IV;
        NumMatrix _V;
        NumVector _D;
    };
  
    enum {
        MAXVALENCE = 50
    };
  
    enum {
        INTERIOR      = 1,
        BOUNDARY      = 2,
        CORNERCONVEX  = 3,
        CORNERCONCAVE = 4
    };
  
public:
    CCSubMatLib() {;}
    ~CCSubMatLib() {;}
    //ops
    int setup(istream& sin);
    Entry& chooseMatrix(int flag, int k);
  
protected:
    vector<Entry> _interior;
    vector<Entry> _boundary;
    vector<Entry> _cornerconvex;
    vector<Entry> _cornerconcave;
};

#endif
