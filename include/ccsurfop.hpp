#ifndef _CCSURFOP_HPP_
#define _CCSURFOP_HPP_

#include "ccsurf.hpp"


class CCSurfOp {
public:
    typedef pair<int,int> intpair;
    typedef CCSurf::DataLevel DataLevel;
  
    static int share(    GpMesh&, int lvl, DataLevel& wk);
  
    static int midsub(   GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
    static int subdivide(GpMesh&, int lvl, double fs, DataLevel& fm, DataLevel& to, CCSubMatLib&);
    static int limit(    GpMesh&, int lvl, DataLevel& fm, DataLevel& to, CCSubMatLib&);
    static int normal(   GpMesh&, int lvl, DataLevel& fm, DataLevel& to, CCSubMatLib&);
    static int du(       GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
    static int dv(       GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
    static int duu(      GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
    static int duv(      GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
    static int dvv(      GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
    //todo: addto, scale, assign
    static int rotate(   GpMesh&, int lvl, DataLevel& wk, Quaternion r);
    static int shift(    GpMesh&, int lvl, DataLevel& wk, Point3 s);
  
    static int eno2index(int lvl, int eno, int t, int& i, int& j);
  
    static int getcenter( GpMesh&, int lvl, DataLevel& wk, int V, Point3&);
    static int putcenter( GpMesh&, int lvl, DataLevel& wk, int V, Point3&);
    static int getonering(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
    static int putonering(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
    static int gethalfring(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
    static int puthalfring(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
};


#endif
