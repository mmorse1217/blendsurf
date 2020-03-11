#ifndef _BDULIB_HPP_
#define _BDULIB_HPP_

#include "nummatrix.hpp"


class BdULib {
public:
    class Entry {
    protected:
        NumMatrix kU, kIU;
    public:
        NumMatrix& U() { return kU; }
        NumMatrix& IU() { return kIU; }
    };
    enum {
        MAXVALENCE = 13
    };
    BdULib() {;}
    ~BdULib() {;}
    //ops
    int setup(istream& sin);
    Entry& chooseEntry(int k) { return _evec[k]; }
protected:
    vector<Entry> _evec;
};


#endif
