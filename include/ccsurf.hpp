#ifndef _CCSURF_HPP_
#define _CCSURF_HPP_

#include "gpmesh.hpp"
#include "quaternion.hpp"
#include "nummatrix.hpp"
#include "ccsubmatlib.hpp"


//-------------------------------------------------
template <class F>
class CCRect {
    //each face is represented with a ccrect in a data level.
    //the indeces go from -1,-1 to m, n so the real size is 
    //(m+2)X(n+2) although looks like mxn to the user.this is
    //to allow sharing.
public:
    int _m, _n;
    F* _data;
public:
    CCRect(int m=0, int n=0): _m(m), _n(n) {
        _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
    }
    CCRect(const CCRect& C): _m(C._m), _n(C._n) {
        _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
        memcpy( _data, C._data, (_m+2)*(_n+2)*sizeof(F) );
    }
    ~CCRect() {
        delete[] _data; _data = NULL;
    }
    CCRect& operator=(const CCRect& C) {
        delete[] _data; _data = NULL;
        _m = C._m; _n=C._n;
        _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
        memcpy( _data, C._data, (_m+2)*(_n+2)*sizeof(F) );
        return *this;
    }
    void resize(int m, int n)  {
        if(_m!=m || _n!=n) {
            delete[] _data; _data = NULL;
            _m = m; _n = n;
            _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
        }
    }
    const F& operator()(int i, int j) const  { 
        assert( i>=-1 && i<_m+1 && j>=-1 && j<_n+1 );
        return _data[(i+1)+(j+1)*(_m+2)];
    }
    F& operator()(int i, int j)  { 
        assert( i>=-1 && i<_m+1 && j>=-1 && j<_n+1 );
        return _data[(i+1)+(j+1)*(_m+2)];
    }
    int m() const { return _m; }
    int n() const { return _n; }
    F* data() const { return _data; }
};

//-------------------------------------------------
//cc subdivision surface
class CCSurf {
public:
    enum { EVAL_VALUE=1,     //evaluate the values
           EVAL_1ST_DERIV=2,     //the first derivative
           EVAL_2ND_DERIV=4 };   //the second derivative
    enum { MAXLEVELNUM=7 };
    typedef pair<int,int> intpair;
    typedef vector< CCRect<Point3> > DataLevel; //vector of faces
  
protected:
    CCSubMatLib* _subdivMatrices;   /// Catmull-clark subdivision surface library, 
                                    /// stores information about subdivision matrices
    double _flatS;      /// flatness parameter used at concave corners
                        /// Biermann H., Levin A., Zorin D., Piecewise Smooth
                        /// Subdivision Surfaces with Normal Control
    GpMesh _gpmesh;     /// Mesh object
    int _inlvl;         /// input level
    vector<DataLevel> _pos; /// vector of levels
public:
    CCSurf() {;}
    ~CCSurf() {;}
    
    /**
     * Construct Catmull-Clark surface
     * @param double    _iniFlatS           flatness parameter used at concave corners
     * @param istream&  sin                 input stream to read data from
     * @param int       flipNormal          flag to flip the normal direction of faces
     * @param Point3    scale               scaling factor
     * @param void*     _iniSubdivMatrices  subdivision surface matrices
     */
    int setup(double _iniflatS, istream& sin, int flipNormal, Point3 scale, void* _iniSubdivMatrices);
    /**
     * Outputs CC surface informations to given stream.
     * @param ostream&  sout    output stream
     */
    int dump(ostream& out, int lvl);
    /**
     * Evaluate CC surface using cubic spline interpolation.
     * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
     * @param int       fid     Global face number (F)
     * @param double*   uv      cd poistion in cd coordinates of face F.
     * @param Point3*   ret     surface position|derivatives
     * 
     * The expected dimension of uv is 2.
     * The expected dimension of ret depends on the flags (1 | 3 | 6)
     */
    int eval(int flags, int fid, double* uv, Point3*);
    /**
     * Access method for ctr method of gpmesh.
     * Returns the averaged center location of the mesh.
     * @return Point3 center location of mesh
     */
    Point3 ctr() { return _gpmesh.ctr(); }
    /**
     * Access method for bbx method of gpmesh.
     * Calculates the bounding box of mesh.
     * @param Point3&   bbmin   lower corner of the bounding box.
     * @param Point3&   bbmax   upper corner of the bounding box.
     */
    void bbox(Point3& bbmin, Point3& bbmax) { _gpmesh.bbox(bbmin,bbmax); }//TODO: change
    /**
     * Rotates the mesh using the given quaternion matrix.
     * @param Quaternion& qr    quaternion matrix for transformation.
     */
    int rotate(const Quaternion& qr);
    /**
     * Shifts the mesh by the given amount in x,y,z directions.
     * @param Point3&   sh      amount of translation.
     */
    int shift(const Point3& sh);
    /**
     * Subdivide CC surface
     * @param int           lvl     subdivision level
     * @param CCSubMatlib&  sb      CC subdivision matrices
     */  
    int subdivide(int lvl, CCSubMatLib&);
    //access
    /**
     * Returns number of faces in the original mesh (before subdivision).
     * @return int number of faces
     */
    int numFaces() { return _gpmesh.numFaces(); }
   /**
     * Returns number of vertices in the original mesh (before subdivision).
     * @return int number of vertices
     */
    int numVertices() { return _gpmesh.numVertices(); }
    /**
     * Access method for mesh representation
     * @return GpMesh& reference to mesh representation
     */
    GpMesh& gpmesh() { return _gpmesh; }
    /**
     * Access method for subdivision matrices.
     * @return CCSubMatlib* subdivision matrices.
     */
    CCSubMatLib* subdivMatrices() { return _subdivMatrices; }
    /**
     * Returns initial subdivision level
     * @return int  _inlvl
     */
    int inlvl() { return _inlvl; }
    /**
     * Returns datalavel corresponding to given lvl
     * @param int       lvl     level for which data is required
     */
    DataLevel& pos(int lvl) { return _pos[lvl]; }
    /**
     * returns true if CC surface at given level exists
     */
    bool validpos(int lvl) { return _pos[lvl][0].m()==pow2(lvl)+1; }
    /**
     * returns the last subdivision level
     */
    int lastposlvl() {	 int la = inlvl();	 while(la+1<MAXLEVELNUM && validpos(la+1)) la++;	 return la;  }
};

#endif
