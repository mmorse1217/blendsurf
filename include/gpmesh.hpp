#ifndef _GPMESH_HPP_
#define _GPMESH_HPP_

#include "common.hpp"
#include "vec3t.hpp"

//-------------------------------------------------
//general purpose mesh
class GpMesh {
public:
    typedef pair<int,int> intpair;
    typedef pair<double,double> dblpair;
  
    GpMesh() {;}
    ~GpMesh() {;}
    /**
     * Setup the mesh structure
     * @param istream&      sin         input stream to read data from
     * @param int           flipNormal  flag to flip the normal direction of faces
     * @param Point3        scale       scaling factor
     */
    int setup(istream& sin, int flipNormal, Point3 scale);
    /**
     * Ouputs mesh information to given output stream.
     * @param ostream&  sout    output stream
     */
    int dump(ostream& sout);
    /**
     * Returns the number of faces contained in mesh.
     */
    int numFaces() { return _Fe2Fe.size(); }
    /**
     * Returns the number of vertices contained in mesh.
     */
    int numVertices() { return _Vf2Fv.size(); }
    /**
     * Returns the valence of a vertex V.
     * @param   int   V   Global vertex number
     * @return  int   valence of V.
     */
    int valence(int V) { return _Vf2Fv[V].size(); }
    /**
     * Returns the averaged center location of the mesh.
     * @return Point3 center location of mesh
     */
    Point3 ctr();
    /**
     * Calculates the bounding box of mesh.
     * @param Point3&   bbmin   lower corner of the bounding box.
     * @param Point3&   bbmax   upper corner of the bounding box.
     */
    void bbox(Point3&, Point3&);
    //access
    /**
     * Returns a reference to coordinate of vertex V.
     * @param int       V       Global vertex number.
     */
    Point3& vpoint(int V)        { assert(V>=0&&V<numVertices()); return _vpoint[V]; }
    /**
     * returns the pair containing neighbor global face number F and local edege e
     * for a given Global face F and local edge e.
     */
    intpair& Fe2Fe(int F, int e) { assert(F>=0&&F<numFaces()); return _Fe2Fe[F][e]; }
    /**
     * returns the pair containing corresponding Global face F and local vertex v
     * number for a given global vertex V, local face f.
     */
    intpair& Vf2Fv(int V, int f) { assert(V>=0&&V<numVertices()); return _Vf2Fv[V][f]; } 
    /**
     * returns the pair containing corresponding global vertex V and local face f
     * number for a given global face F and local vertex v.
     */
    intpair& Fv2Vf(int F, int v) { assert(F>=0&&F<numFaces()); return _Fv2Vf[F][v]; }
    ///  INTERIOR_VERTEX | BOUNDARY_VERTEX
    int& boun(int V)             { assert(V>=0&&V<numVertices()); return _boun[V]; }
    /// CREASE_VERTEX | CONVEX_VERTEX | CONCAVE_VERTEX
    int& corner(int V)           { assert(V>=0&&V<numVertices()); return _corner[V]; }
    /**
     * Flags for vertex type.
     * Interior or boundary vertex
     */
    enum { INTERIOR_VERTEX = 0, BOUNDARY_VERTEX = 1 };
    /**
     * Flags for boundary vertex type
     * crease(smoth boundary) | convex | concave
     */
    enum { CREASE_VERTEX = 0, CONVEX_VERTEX = 1, CONCAVE_VERTEX = 2 };
public:
    //e: 0,1,2,3
    //v: 0,1,2,3
    /// Returns next vertex number
    int nextvno(int vno) { return (vno+1)%4; }
    /// Returns previous vertex number
    int prevvno(int vno) { return (vno+3)%4; }
    /// Returns next edge number
    int nexteno(int eno) { return (eno+1)%4; }
    /// Returns previous edge number
    int preveno(int eno) { return (eno+3)%4; }
    /// Returns the vertex number that given edge starts from.
    int fromvno(int eno) { return eno; }
    /// Returns the vertex number that given edge goes to
    int gotovno(int eno) { return (eno+1)%4; }
    /// Returns the edge number which has this vertex at its end
    int enteeno(int vno) { return (vno+3)%4; }
    /// Returns the edge number whihc has this vertex at its start
    int leaveno(int vno) { return vno; }

    // Return vector of interior points
    vector<Point3> get_interior_points(){
        return _interior_points;
    }

    // Return vector of interior point orientations
    vector<int> get_intpt_orientation(){
        return _intpt_orientation;
    }
    // Access vector of group ids, read in from options files 
    vector<int> get_group_ids(){
        return _gids;
    }
    vector<int> get_face_group_ids(){
        return _face_group_ids;
    }
protected:
    vector< Point3 > _vpoint;
    vector< vector<intpair> > _Fe2Fe;
    vector< vector<intpair> > _Vf2Fv;
    vector< vector<intpair> > _Fv2Vf;
    vector< int > _boun;
    vector< int > _corner;
    vector<Point3> _interior_points;
    vector<int> _intpt_orientation;
    vector<int> _gids; //TODO: remove this
    vector<int> _face_group_ids;
};

#endif
