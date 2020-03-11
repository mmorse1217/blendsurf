#ifndef _BDSURF_HPP_
#define _BDSURF_HPP_

#include "common.hpp"
#include "quaternion.hpp"
#include "nummatrix.hpp"
#include "ccsurf.hpp"
#include "vec2t.hpp"
#include "bdulib.hpp"

/**
 * Manifold surface definition
 */

class BdSurf {
public:
  /**
   * Flags for evaluation type
   * EVAL_VALUE : Evaluate only value of a point on the surface
   * EVAL_1ST_DERIV : Evaluate 1st derivative
   * EVAL_2ND_DERIV : Evaluate 2nd derivative
   * EVAL_1ST_DERIV and EVAL_2ND_DERIV can not be used without EVAL_VALUE
   */
  enum { EVAL_VALUE=1,  EVAL_1ST_DERIV=2,    EVAL_2ND_DERIV=4  };
  
  typedef pair<int,int> intpair;
  typedef complex<double> dblcomplex;
  typedef vector< CCRect<Point3> > DataLevel;
  
  DataLevel  _ccLimitPos;  
  
protected:

  /**
   * Type of partition of unity function flag
   * ANALYTIC_BLEND_FUNC : Analytic function
   * BSPLINE_BLEND_FUNC : Spline basis functions as pou
   */
  enum { ANALYTIC_BLEND_FUNC=0,   BSPLINE_BLEND_FUNC=1  , OPTIM_BLEND_FUNC=2};
  /// Type of chart flag
  enum { FULLY_COMPLEX_CHART = 0,
	 CHARACTERISTIC_MAP_CHART = 1, 
	 ISODISTANCE_CHART = 2  };
  /// Basis Function Type
  enum { POLYNOMIAL_BASIS = 0, 
	 SPLINE_BASIS = 1 };
  
  ///Type of matrix construction for independent boundary
  enum { SMOOTH_MAT = 0, 
	 CONVEX_MAT = 1,
	 CONCAVE2_MAT = 2,
	 CONCAVE3_MAT = 3};
  
  CCSubMatLib* _subdivMatrices;   /// Catmull-clark subdivision surface library, 
                                  /// stores information about subdivision matrices
  BdULib* _matricesU;             /// Manifold surface library, stores information about 
                                  /// matrices (U) used by manifold surface.
                                  /// Least squares fit : min||Ua-s||^2
                                  /// Ying L., Zorin D. A Simple Manifold Based Construction 
                                  /// of Surfaces of Arbitrary Smoothness, Siggraph'04,p.3
  CCSurf _ccSurf;                 /// Catmull-clark subdivision surface object
  vector< Vector<Point3>* > _controls; /// this stores for each vertex a vector of 
                                  /// polynomial coefficients for the corresponding 
                                  /// chart, or the spline control pts for the chart
  int _subdivCtrlLevel;           /// control level, where points are used to approximate 
                                  /// around extraordinary vertex
  int _blendFuncType;             /// what partition of unity to use, C^\infty or Spline based
  int _chartType;                 /// chart type: FULLYCOMPLEX | CHARACTERISTIC | ISODISTANCE
  int _basisType;                 /// basis type: polynomial | spline with extra ring | 
                                  /// spline with constant derivatives
  bool _preprocess;               /// setup phase preprocesing
  bool _indepBoundaryFit;         /// do a global least squares or separate boundary
  double _flatS;                  /// flatness parameter used at concave corners
                                  /// Biermann H., Levin A., Zorin D., Piecewise Smooth
                                  /// Subdivision Surfaces with Normal Control
  int _max_degree_valence;        /// controls degree for high order 
  int _splineDeg;                 /// degree of spline
  int _splineBlendFuncDeg;        /// degree of spline function used in POU
  double _evalLowerBnd;           /// Lower bound of cut-off value for evaluation [0,1]
  double _evalUpperBnd;           /// Upper bound of cut-off value for evaluation [0,1]
  int _refinement_factor;         /// number of levels to subdivide the original control mesh
  double _interpolation_spacing;     /// sampling rate of surface interpolant
  bool _evaluate_via_interpolant;     /// sampling rate of surface interpolant
public:
  /**
   * Constructor method
   */
  BdSurf(): _subdivMatrices(NULL) {;}
  /**
   * Destructor method
   */
  ~BdSurf();

  double EVAL_UB() {
      return _evalUpperBnd;
  }
  
  //KK : The interface to this method is very heavy I think. It may be better
  // to use a double array for this one or a struct, but array would be cleaner I guess
  // for this case.

  /**
   * Upload the input parameters for manifold surface
   * @param int       _iniSubdivCtrlLevel      control level where points are used to 
   *                                           approximate around extraordinary vertex
   * @param int       _iniBlendFuncType        type of partition of unity function 
   *                                           (ANALYTIC_BLEND_FUNC | BSPLINE_BLEND_FUNC)
   * @param int       _iniChartType            chart type 
   *                                           (FULLYCOMPLEX | CHARACTERISTIC | EQUALDISTANCE)
   * @param int       _iniBasisType            basis function type 
   *                                           (Polynomial | Spline with extra ring | Spline)
   * @param int       _iniPreprocess           setup phase preprocessing
   * @param double    _iniEvalLowerBnd         lower bound for evaluation
   * @param double    _iniEvalUpperBnd         upper bound for evaluation
   * @param int       _iniIndepBoundaryFit      
   * @param double    _iniFlatS                flatness parameter used at concave corners
   * @param int       _inimax_degree_valence  
   * @param int       _iniSplineDeg            degree of spline curve
   * @param int       _iniSplineBlendFuncDeg   degree of spline function used in POU
   */
  int setParams(int _iniSubdivCtrlLevel, int _iniBlendFuncType,
		int _iniChartType, int _iniBasisType, int _iniPreprocess, 
		double _iniEvalLowerBnd, double _iniEvalUpperBnd, int _iniIndepBoundaryFit, 
		double _iniFlatS, int _inimax_degree_valence, int _iniSplineDeg, 
		int _iniSplineBlendFuncDeg, 
        double interpolation_spacing,
        int evaluate_via_interpolant,
        int _refinement_factor=0,
		int loadW=0); /*** for loading precomputed W matrices DZ **/
  
  /**
   * Construct the manifold surface
   * @param istream&  sin                 input stream to read data from
   * @param int       flipNormal          flag to flip the normal direction of faces
   * @param Point3    scale               scaling factor
   * @param void*     _iniSubdivMatrices  subdivision surface matrices
   * @param void*     _iniMatricesU       manifold surface matrices
   */
  int setup(istream&, int flipNormal, Point3 scale, void* _iniSubdivMatrices, void* _iniMatricesU);
  
  /**
   * Access method for underlying catmull-clark surface
   * @return CCSurf& reference to catmull-clark surface representation
   */
  CCSurf& ccSurf() { return _ccSurf; }
  
  /**
   * Access method for mesh representation
   * @return GpMesh& reference to mesh representation
   */
  GpMesh& gpmesh() { return _ccSurf.gpmesh(); }

  /**
   * Returns number of faces in the original mesh (before subdivision).
   * @return int number of faces
   */
  int numFaces() { return gpmesh().numFaces(); }

  /**
   * Returns number of vertices in the original mesh (before subdivision).
   * @return int number of vertices
   */
  int numVertices() { return gpmesh().numVertices(); }

  /**
   * Number of subdivision steps for constructing catmull-clark 
   * surface representation.
   * @return int subdivision control level
   */
  int subdivCtrlLevel() { return _subdivCtrlLevel;}

  /**
   * Returns valence of a vertex V.
   * @param int V vertex number
   * @return int valence of V
   */
  int valence(int V) { return gpmesh().valence(V); }

  /**
   * Access method for ctr method of gpmesh.
   * Returns the averaged center location of the mesh.
   * @return Point3 center location of mesh
   */
  Point3 ctr() { return gpmesh().ctr(); }

  /**
   * Access method for bbx method of gpmesh.
   * Calculates the bounding box of mesh.
   * @param Point3&   bbmin   lower corner of the bounding box.
   * @param Point3&   bbmax   upper corner of the bounding box.
   */
  void bbox(Point3& bbmin, Point3& bbmax) { gpmesh().bbox(bbmin,bbmax); }

  /**
   * Outputs vertex locations to given stream.
   * @param ostream&  sout    output stream to direct vertex locations
   */
  int dump(ostream& sout);

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
   * Evaluate (U wrt C) * (C wrt X) to get U wrt X = R
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       udof    number of degrees of freedom of U, C, R
   * @param double*   U       position|derivatives wrt C
   * @param double*   C       position|derivatives wrt X
   * @param double*   R       position|derivatives wrt X
   * 
   * The expected dimensions for U,C,R are as follows.
   * U,R shall contain variable number of Point3, Point2 or Point1 depending on
   * the dof specified. If dof=3, They should contain Point3. The number of
   * points contained depends on the flags variable.
   *      * If only evaluation at a point is required (EVAL_VALUE) then
   *      U and R must contain a single point.
   *      * If value and 1st derivative are required then
   *      U and R must contain 3 points (u, u_c, u_d) and (u, u_x, u_y)
   *      * If value, 1st derivative and 2nd derivatrive are required then
   *      U and R must contain 6 points (u, u_c, u_d, u_cc, u_cd, u_dd) and
   *      (u, u_x, u_y, u_xx, u_xy, u_yy)
   * Since PointX is actually a double array with X elements, for dof=3 and 
   * evaluation, 1st deriv., 2nd deriv, U,R must contain 6 Point3 that is double[18].
   * C contains partial derivatives of c,d wrt x,y. Therefore it contains
   * c, d, c_x, c_y, d_y, c_xx, d_xx, c_xy, d_xy, c_yy, d_yy for value, 1st deriv
   * and 2nd deriv.
   */
  int compose(int flags, int udof, double* U, double* C, double* R);

  /**
   * Transformation from xy coordinate system of global vertex V to cd 
   * coordinate system in face f of global vertex V.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates.
   *                          Coordinate system in domain D of Figure 1 as 
   *                          curved star shape, composed of K curved wedges,
   *                          where K is valence of V. The origin is the vertex
   *                          which corresponds to the chart center. The xy
   *                          coordinate system is defined on 360deg curved star.
   * @param int&      f       face number [0,valence(V)]
   * @param double*   cd      position|derivatives in cd coordinates.
   *                          Coordinate system in domain S of Figure 1 as
   *                          regular star shape composed of K regular wedges,
   *                          where K is valence of V. The origin is the vertex
   *                          which corresponds to the chart center. The cd
   *                          coordinate system is defined in a single wedge,
   *                          the position in S domain is defined by the local
   *                          face number f.
   * 
   * In this method transformation from Vxy domain to Vfcd domain is done. 
   * Vxy is xy coordinate system centered on chart number V. Since charts are 
   * defined per vertex, chart number V = vertex number V. Vfcd if cd 
   * coordinate system centered on vertex(chart) V.
   * 
   * The expected dimension of xy is 2.
   * The expected dimension of cd depends on the flags (2 | 6 | 12)
   */
  int Vxy2Vfcd(int flags, int V, double* xy, int& f,  double* cd);

  /**
   * Transformation from cd coordinate system in face f of global vertex V to 
   * xy coordinate system of global vertex V.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param int&      f       face number [0,valence(V)]
   * @param double*   cd      position in cd coordinates.
   *                          Coordinate system in domain S of Figure 1 as
   *                          regular star shape composed of K regular wedges,
   *                          where K is valence of V. The origin is the vertex
   *                          which corresponds to the chart center. The cd
   *                          coordinate system is defined in a single wedge,
   *                          the position in S domain is defined by the local
   *                          face number f.
   * @param double*   xy      position|derivatives in xy coordinates.
   *                          Coordinate system in domain D of Figure 1 as 
   *                          curved star shape, composed of K curved wedges,
   *                          where K is valence of V. The origin is the vertex
   *                          which corresponds to the chart center. The xy
   *                          coordinate system is defined on 360deg curved star.
   * 
   * This method is opposite if Vxy2Vfcd
   * 
   * The expected dimension of cd is 2.
   * The expected dimension of xy depends on the flags (2 | 6 | 12)
   */
  int Vfcd2Vxy(int flags, int V, int f, double* cd, double* xy);


  int Vfcd2Vxy_val(int flags, int K,int bountype, int contype,  
			   int f, double* cd, double* xy);
  /**
   * Transformation from cd coordinates in face f of global Vertex V to cd
   * coordinates of global face F.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param int&      f       face number [0,valence(V)]
   * @param double*   cd      position in cd coordinates.
   *                          Coordinate system in domain S of Figure 1 as
   *                          regular star shape composed of K regular wedges,
   *                          where K is valence of V. The origin is the vertex
   *                          which corresponds to the chart center. The cd
   *                          coordinate system is defined in a single wedge,
   *                          the position in S domain is defined by the local
   *                          face number f.
   * @param int&      F       Global face number
   * @param double*   st      position|derivatives in cd coordinates.
   *                          Coordinate system in a single regular wedge. The
   *                          origin of the coordinate system is arbitrary and
   *                          independent of chart choice.
   * 
   * In this method transformation from Vfcd domain to Fcd domain is done.
   * Vfcd is the cd coordinate system centered on vertex(chart) V and Fcd is
   * the cd coordinate system defined in face F and centered independent of 
   * the chart choice.
   * 
   * The expected dimension of cd of Vfcd is 2.
   * The expected dimension of cd of Fcd depends on the flags (2 | 6 | 12)
   */
  int Vfcd2Fcd(int flags, int V, int f, double* cd, int& F, double* st);

  /**
   * Transformation from cd coordinate system of global face F to cd 
   * coordinate system in face f of global vertex V corresponding to vertex 
   * 'v' of F.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int&      F       Global face number
   * @param double*   cd      position|derivatives in cd coordinates.
   *                          Coordinate system in a single regular wedge. The
   *                          origin of the coordinate system is arbitrary and
   *                          independent of chart choice.
   * @param int       v       vertex number [0,3] in face F.
   * @param int&      V       Global vertex number
   * @param int&      f       face number [0,valence(V)]
   * @param double*   st      position|derivatives in cd coordinates.
   *                          Coordinate system in domain S of Figure 1 as
   *                          regular star shape composed of K regular wedges,
   *                          where K is valence of V. The origin is the vertex
   *                          which corresponds to the chart center. The cd
   *                          coordinate system is defined in a single wedge,
   *                          the position in S domain is defined by the local
   *                          face number f.
   * 
   * This method is the opposite of Vfcd2Fcd.
   * 
   * The expected dimension of cd of Fcd is 2.
   * The expected dimension of cd of Vfcd depends on the flags (2 | 6 | 12)
   */
  int Fcd2Vfcd(int flags, int F, double* cd, int v, int& V, int& f, double* st);

  /**
   * Application of Partition Of Unity Algorithm on cd coordinate system 
   * in face f of a global Vertex V.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param double*   cd      position in cd coordinates
   * @param double    lb      Lower bound for POU algorithm
   * @param double    ub      Upper bound for POU algorithm
   * @param double*   pou     position|derivatives after POU
   * 
   * The expected dimension of cd is 1.
   * The expected dimension of pou depends on the flags (1 | 3| 6)
   */    
  int blendFuncEval(int flags, double* cd, double lb, double ub, double* pou);    

  /**
   * Given x,y=ub in cd coordinates, calculates the corresponding x and y 
   * coordinates in a face of vertex V (the result is the same for all faces) 
   * and updates rad with the absolute value of maximum of calculated x and y.
   * @param int       V       Global vertex number
   * @param double    ub      x,y position in cd coordinates (x=y=ub)
   * @param double&   rad     abs(max. value of position) in xy coordinates
   */
  int paramBound(int V, double ub, double& rad);

  /**
   * Evaluate surface (that is calculate the 3D surface point location and 
   * derivatives) from given position in xy coordinate system
   * of global vertex V using contirbutions of all charts that share the face
   * at which xy is located.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates
   * @param Point3*   ret     surface position|derivatives in R3.
   * 
   * The expected dimension of xy is 2.
   * the expected dimension of ret depends on the flags (1 | 3 | 6)
   */
  int eval(int flags, int V, double* xy, Point3* ret, bool skip_interp=false);
    int eval2(int flags, int V, double* xy, Point3* ret);
  /**
   * Evaluates the surface position and derivatives using given position
   * in x,y coordinate system of a global vertex V usign only the
   * contribution from chart of vertex V. Makes call to charteval.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates
   * @param Point3*   ret     surface position|derivatives in R3.
   * 
   * The expected dimension of xy is 2.
   * The expected dimension of ret depends on flags (1 | 3 | 6)
   */
  int singleChartEval(int flags, int V, double* xy, Point3* ret);

  /**
   * Evaluates the surface at all points. used for timing purposes only.
   * @param           lvl     number of subdivision steps
   * @param           gen     rendering flag (full | 1/2 charts | 
   *                          extraordinary full)
   */
  int evalall(int lvl, int gen);
  
protected:
  /**
   * Carries out the partition of unity calculation in a single dimension.
   * Higher level partition of unity calculators make calls to this method.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param double*   u       position (x or y) in cd coordinates
   * @param double    lb      Lower bound for POU algorithm
   * @param double    ub      Upper bound for POU algorithm
   * @param double*   res     position|derivatives after POU
   * 
   * The expected dimension of res depends on flags (1 | 2 | 3)
   * (val, 1st deriv, 2nd deriv)
   */
  int blendFunc1D(int flags, double u, double lb, double ub, double* res);
  
  /**
   * Partition of unity calculation using analytic function.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param double    t       position in cd coordinates normalized by actual 
   *                          length of basis function (ub - lb)
   * @param double    lb      Lower bound for POU algorithm
   * @param double    ub      Upper bound for POU algorithm
   * @param double*   res     position|derivatives after POU
   * 
   * The expected dimension of res depends on flags (1 | 2 | 3)
   * (val, 1st deriv, 2nd deriv)
   */
  int  analyticBlendFunc1D(int flags, double t, double lb, double ub, double* res);


  /**
   * Partition of unity calculation using the optimal blending function constructed using the lower bound proof.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param double    t       position in cd coordinates normalized by actual 
   *                          length of basis function (ub - lb)
   * @param double    deg     deg of blending function
   * @param double*   res     position|derivatives after POU
   * 
   * The expected dimension of res depends on flags (1 | 2 | 3)
   * (val, 1st deriv, 2nd deriv)
   */
  int  optimBlendFunc1D(int flags, int deg, double t, double scale, double* res);

  /**
   * Partition of unity calculation using spline basis functions
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       deg     degree of spline basis functions
   * @param double    t       position in cd coordinates normalized by actual 
   *                          length of basis function (ub - lb)
   * @param double    scale
   * @param double*   res     Coordinates after POU
   * 
   * The expected dimension of res depends on flags (1 | 2 | 3)
   * (val, 1st deriv, 2nd deriv)
   */
  int splineBlendFunc1D(int flags, int deg, double t, double scale, double* res);
  
  /**
   * Evaluates the surface position and derivatives using given position 
   * in x,y coordinate system of a global vertex V using only the
   * contribution from chart of vertex V. Makes call to evalC2 or evalSp or
   * evalPoly.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates
   * @param Point3*   ret     surface position|derivatives in R3.
   * 
   * The expected dimension of xy is 2.
   * The expected dimension of ret depends on flags (1 | 3 | 6)
   */
  int chartEval(int flags, int V, double* xy, Point3* ret, int flag);
  
  /**
   * Evaluate surface using polynomial basis functions.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates
   * @param Point3*   ret     surface position|derivatives in R3.
   * 
   * The expected dimension of xy is 2.
   * The expected dimension of ret depends on flags (1 | 3 | 6)
   */
  int evalPoly(int flags, int V, double* xy, Point3* ret, int flag);
  
  /**
   * Evaluate surface using 3rd degree spline basis functions. Used only for
   * C2 and closed surfaces.
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates
   * @param Point3*   ret     surface position|derivatives in R3.
   * 
   * The expected dimension of xy is 2.
   * The expected dimension of ret is 3.
   */
  inline int evalCubicBSpline(int flags, int V, double* xy, Point3* ret);
  
  /**
   * Evaluate surface using spline basis functions.
   * @param int       flags   EVAL_VALUE | EVAL_1ST_DERIV | EVAL_2ND_DERIV
   * @param int       V       Global vertex number
   * @param double*   xy      position in xy coordinates
   * @param Point3*   ret     surface position|derivatives in R3.
   * 
   * The expected dimension of xy is 2.
   * The expected dimension of ret depends on flags (1 | 3 | 6)
   */
  int evalBSpline(int flags, int V, double* xy, Point3* ret);
  
  /**
   * Calculate cd, xy and surface limit positions of control mesh.
   * @param int               V       Global vertex number
   * @param vector<Point2>&   tmpcd   cd coordinates in Vfcd coordinate system
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   */
  void gatherData(int V, vector<Point2>& tmpcd, vector<Point2>& tmpxy, 
		  vector<Point3>& tmpval);

  /**
   * Project surface points using _matricesU of appropriate valence.
   * @param int               K       valence
   * @param vector<Point3>&   tmpval  limit position from CC surface
   */
  void preProcessPoly(int K, vector<Point3>& tmpval);
  
  /**
   * Calculate limit position of CC mesh into _ccLimitPos and initialize
   * _controls vector with numv number of NULL pointers.
   * @param int       numv    number of vertices in gpmesh.
   */
  void constructControls(int numv);
  
  /**
   * Solve Least Squares Problem using polynomial basis functions for a vertex
   * @param int               V       Global vertex number
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   * @param vector<Point2>&   tmsqD   weight for each point
   */
  void constructPoly(int V, vector<Point2>& tmpxy, vector<Point3>& tmpval);

    
  /**
   * Solve Least Squares Problem using polynomial basis functions for an
   *  interior vertex
   * @param int               V       Global vertex number
   * @param int               DEG     degree of monomials used in fitting
   * @param int               cm      number of vertices (tmpxy.size())
   * @param int               cn      number of monomials used in fitting
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   * @param vector<Point2>&   tmsqD   weight for each point
   */
  void constructPolyInterior(int V, int DEG, int cm, int cn, 
			     vector<Point2>& tmpxy, vector<Point3>& tmpval); 
   
  /**
   * Solve Least Squares Problem using polynomial basis functions for a
   * boundary vertex
   * @param int               V       Global vertex number
   * @param int               DEG     degree of monomials used in fitting
   * @param int               cm      number of vertices (tmpxy.size())
   * @param int               cn      number of monomials used in fitting
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   * @param vector<Point2>&   tmsqD   weight for each point
   */
  void constructPolyBoundary(int V, int DEG, int cm, int cn, 
			     vector<Point2>& tmpxy, vector<Point3>& tmpval); 

  /**
   * Solve Least Squares Problem using polynomial basis functions for a
   * corner vertex
   * @param int               V       Global vertex number
   * @param int               DEG     degree of monomials used in fitting
   * @param int               cm      number of vertices (tmpxy.size())
   * @param int               cn      number of monomials used in fitting
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   * @param vector<Point2>&   tmsqD   weight for each point
   */
  void constructPolyCorner(int V, int DEG, int cm, int cn, 
			   vector<Point2>& tmpxy, vector<Point3>& tmpval);

  

  /**
   * Initialize x and y values for monomial calculation.
   * @param int           n       degree of monomial
   * @param double        x       x coordinate of tmpxy of a vertex
   * @param double        y       y coordinate of tmpxy of a vertex
   * @param double*       xs      
   * @param double*       ys
   */
  void initPolyXsYs(int n, double x, double y, double* xs, double* ys);
  
  /**
   * Calculates the control points for vertex V using pinv(Monomial Matrix).
   * @param int               V       Global vertex number
   * @param int               cn      number of monomials used in fitting
   * @param int               cm      number of vertices (tmpxy.size())
   * @param NumMatrix&        M       stores evaluation of all monomials
   * @param vector<Point3>&   Val     limit poisition from CC surface
   * @param vector<Point2>&   D       weight for each point
   */
  void updatePolyControlPoint(int V, /*int cn, int cm,*/ NumMatrix& M, 
			  vector<Point3>& Val);
  
  
  /**
   * Solve Least Squares Problem using spline basis functions which results
   * in a set of control points that define a surface patch per chart.
   * Calls different construction routines based on the properties of the vertex
   * @param int               V       Global vertex number
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   */
  void constructSpline(int V, vector<Point2>& tmpxy,  vector<Point3>& tmpval);
  
  /**
   * Creates nth degree spline basis functions in one dimension.
   * @param double    u       point where the function is evaluated
   * @param int       k       grid point that marks the end of the segment pt falls in
   * @param int       n       degree of spline 
   * @param int       m       number of gridpts 
   * @param double*   N       array to store the values at each grid pt
   */
  inline int SplineBasis(double u, int k, int n, int m,double* N);
  
  /**
   * Initialize variables related to spline fitting
   * @param int               V       Global vertex number
   * @param int&              degk    spline degree   
   * @param int&              cn      number of data points on chart
   * @param int&              cm      total number of grid points used in fitting
   * @param double&           xmax    max x-coord of the grid  
   * @param double&           ymax    max y-coord of the grid
   * @param double&           totalx  length of grid in x direction
   * @param double&           totaly  length of grid in y direction
   * @param double&           xmin    min x-coord of the grid
   * @param double&           ymin    min y-coord of the grid
   * @param int&              nrgrid  nr of grid pts
   */
  void initSplineVars(int V, int & degk, int &cn, int &cm, double &xmax, 
		      double &ymax, double &totalx, double &totaly,
		      double &xmin, double &ymin, int &ex, int &nrgrid);
 
   /**
   * Compute control points that define the chart of an interior vertex
   * or a boundary vertex where the boundary of the chart is dependent on 
   * the interior. Regular spline surface fitting procedure  based on
   * least squares with no extra constraints.
   * @param int               V       Global vertex number
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   */
  void constructSplineInterior(int V, vector<Point2>& tmpxy, 
			       vector<Point3>& tmpval);
  
  /**
   * Compute control points that define the chart of an boundary vertex
   * where the interior of the chart is dependent on the boundary.
   * The curve defining the boundary is solved for first and the rest
   * of the surface is fit to this boundary. Each step requires a least squares
   * fit. Calls other subroutines for this purpose.
   * @param int               V       Global vertex number
   * @param int               bpts    Nr of data points on the boundary
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit position from CC surface
   */
  void constructSplineIndepSmooth(int V, int bpts, vector<Point2>& tmpxy, 
				  vector<Point3>& tmpval);
  
  /**
   * Compute control points that define the chart of a convex corner
   * where the interior of the chart is dependent on the boundary.
   * Two curves defining the two piecewise smooth boundaries of the corner are
   * solved for first. The rest of the surface is fit to these boundaries.
   * Each step requires a least squares fit. Calls other subroutines for this purpose.
   * @param int               V       Global vertex number
   * @param int               bpts    Nr of data points on the boundary
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   */
  void constructSplineIndepConvex(int V, int bpts, vector<Point2>& tmpxy, 
				  vector<Point3>& tmpval);
  
  /**
   * Compute control points that define the chart of a concave corner
   * where the interior of the chart is dependent on the boundary.
   * Two curves defining the two piecewise smooth boundaries of the corner are
   * solved for first. The rest of the surface is fit to these boundaries.
   * Each step requires a least squares fit. Calls other subroutines for this purpose.
   * @param int               V       Global vertex number
   * @param int               bpts    Nr of data points on the boundary
   * @param vector<Point2>&   tmpxy   xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval  limit poisition from CC surface
   */
  void constructSplineIndepConcave(int V, int bpts, vector<Point2>& tmpxy, 
				   vector<Point3>& tmpval);

  /**
   * Compute control points that define an any degree spline boundary curve
   * for a boundary chart where the boundary is independent
   * @param int               degk
   * @param int               nrgrid    nr of grid pts
   * @param int               bpts      nr of data points on the boundary
   * @param int               cn        nr of data points on chart
   * @param double            totalx    length of grid in x direction
   * @param double            xmin      min x-coord of the grid
   * @param vector<Point2>&   tmpxy     xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval    limit poisition from CC surface
   * @param vector<Point3>&   bcoef     vector that stores control points 
   *                                    for the boundary(ies)
   */
  void computeSpBounCoefSmooth(int degk, int nrgrid, int bpts, int cn, 
			       double totalx, double xmin, vector<Point2> &tmpxy, 
			       vector<Point3> &tmpval, vector<Point3> &bcoef);

 /**
   * Compute control points that define 2 any degree spline boundary curves
   * for a convex corner  chart where the boundary is independent. Then reorder the 
   * boundary coefficients in one vector so that it can be used in matrix computations
   * @param int               degk
   * @param int               nrgrid     nr of grid pts
   * @param int               bpts       nr of data points on the boundary
   * @param int               cn         nr of data points on chart
   * @param double            totalx     length of grid in x direction
   * @param double            xmin       min x-coord of the grid
   * @param double            totaly     length of grid in y direction
   * @param double            ymin       min y-coord of the grid
   * @param vector<Point2>&   tmpxy      xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval     limit poisition from CC surface
   * @param vector<Point3>&   bcoef      vector that stores control points for the boundary(ies)
   */
  void computeSpBounCoefConvex(int degk, int nrgrid, int bpts, int cn, 
			       double totalx, double xmin, double totaly, 
			       double ymin, vector<Point2> &tmpxy, 
			       vector<Point3> &tmpval, vector<Point3> &bcoef);

 /**
   * Compute control points that define 2 any degree spline boundary curves
   * for a cocave corner  chart where the boundary is independent. Then reorder the 
   * boundary coefficients in one vector so that it can be used in matrix computations
   * @param int               degk
   * @param int               nrgrid     nr of grid pts
   * @param int               bpts       nr of data points on the boundary
   * @param int               cn         nr of data points on chart
   * @param double            totalx     length of grid in x direction
   * @param double            xmin       min x-coord of the grid
   * @param double            totaly     length of grid in y direction
   * @param double            ymin       min y-coord of the grid
   * @param vector<Point2>&   tmpxy      xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval     limit poisition from CC surface
   * @param vector<Point3>&   bcoef      vector that stores control points for the boundary(ies) 
   */
  void computeSpBounCoefConcave(int degk, int nrgrid, int bpts, int cn, 
				double totalx, double xmin, double totaly,double ymin,
				vector<Point2>&tmpxy,vector<Point3>&  tmpval,
				vector<Point3>& bcoef);
  
 /**
   * Compute control points that define a cubic boundary curve
   * for a boundary chart where the boundary is independent. 
   * This is necessary when the number of data points on boundary 
   * is not enough for a prescribed degree spline curve
   * to define the boundary. In this case, compute a cubic curve and use
   * degree elevation to match the prescribed degree.
   * @param int               bpts     nr of data points on the boundary
   * @param int               cn       nr of data points on chart
   * @param int               nrgrid   nr of grid pts
   * @param double            totalx   length of grid in x direction
   * @param double            xmin     min x-coord of the grid
   * @param vector<Point2>&   tmpxy    xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval   limit poisition from CC surface
   * @param vector<Point3>&   cubic    control points of the cubic curve stored here
   */
  void computeSpCubicControlsSmooth(int bpts, int cn, int nrgrid, double totalx, 
				    double xmin, vector<Point2>& tmpxy, 
				    vector<Point3>& tmpval, vector<Point3>& cubic);
  
  /**
   * Compute control points that define the 2 cubic boundary curves
   * for a convex corner chart where the boundary is independent
   * This is necessary when the number of data points on boundary 
   * is not enough for a prescribed degree spline curve
   * to define the boundary. In this case, compute a cubic curve and use
   * degree elevation to match the prescribed degree.
   * @param int               bpts     nr of data points on the boundary
   * @param int               cn       nr of data points on chart
   * @param int               nrgrid   nr of grid pts
   * @param double            totalx   length of grid in x direction
   * @param double            xmin     min x-coord of the grid
   * @param double            totaly   length of grid in y direction
   * @param double            ymin     min y-coord of the grid
   * @param vector<Point2>&   tmpxy    xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval   limit poisition from CC surface
   * @param vector<Point3>&   cubicx   control points of the cubic curve on x-direction
   * @param vector<Point3>&   cubicy   control points of the cubic curve on y-direction
   */ 
  void computeSpCubicControlsConvex(int bpts, int cn, int nrgrid, double totalx, 
				    double xmin, double totaly, double ymin,
				    vector<Point2>& tmpxy, vector<Point3>& tmpval, 
				    vector<Point3>& cubicx, vector<Point3>& cubicy);

 /**
   * Compute control points that define the 2 cubic boundary curves
   * This is necessary when the number of data points on boundary 
   * is not enough for a prescribed degree spline curve
   * to define the boundary. In this case, compute a cubic curve and use
   * degree elevation to match the prescribed degree.
   * for a concave corner chart where the boundary is independent
   * @param int               bpts     nr of data points on the boundary
   * @param int               cn       nr of data points on chart
   * @param int               nrgrid   nr of grid pts
   * @param double            totalx   length of grid in x direction
   * @param double            xmin     min x-coord of the grid
   * @param double            totaly   length of grid in y direction
   * @param double            ymin     min y-coord of the grid
   * @param vector<Point2>&   tmpxy    xy coordinates in Vxy coordinate system
   * @param vector<Point3>&   tmpval   limit poisition from CC surface
   * @param vector<Point3>&   cubicx   control points of the cubic curve on x-direction
   * @param vector<Point3>&   cubicy   control points of the cubic curve on y-direction
   */
  void computeSpCubicControlsConcave(int bpts, int cn, int nrgrid, double totalx, 
				     double xmin, double totaly, double ymin,
				     vector<Point2>& tmpxy, vector<Point3>& tmpval, 
				     vector<Point3>& cubicx, vector<Point3>& cubicy);

 /**
   * Create the matrices needed in the construction of spline
   * control points for independent boundary
   * The idea is to force the surface control points to match the boundary
   * curve computed in computeBounCoef<> functions near the boundary. 
   * We fix a number of control points s.t. they satisfy these constraints exactly.
   * The rest are fit using least squares.
   * Let nre be the number of constraints
   * A2 is cn by nre, A1 is cn by cm-nre.
   * B2 is nre by nre, B1 is nre by cm-nre. 
   *
   * @param int               flag      flag to determine whether the setup is for smoooth/convex/concave
   * @param int               degk      prescribed degree of spline
   * @param int               nrgrid    nr of grid pts
   * @param vector<Point2>&   tmpxy     xy coordinates in Vxy coordinate system
   * @param double            xmin      min x-coord of the grid
   * @param double            ymin      min y-coord of the grid
   * @param double            totalx    length of grid in x direction
   * @param double            totaly    length of grid in y direction
   * @param NumMatrix&        A1        stores least squares equations wrt other control points        
   * @param NumMatrix&        A2        stores least squares equations wrt fixed control points
   * @param NumMatrix&        B1        stores the constraint equations wrt other control points
   * @param NumMatrix&        B2        stores the constraint equations wrt fixed control points
   */
  void createSpIndepMatrices(int flag, int degk, int nrgrid, vector<Point2> &tmpxy,
			     double xmin, double ymin, double totalx, double totaly,
			     NumMatrix &A1, NumMatrix &A2, NumMatrix &B1, NumMatrix &B2);
   
 /**
   * Execute matrix computations used in the construction of
   * spline control points for independent boundary
   *
   * #constraint equations
   * [B1 B2] * [p1] = [bcoef]
   *           [p2]   
   * p2 =  inv(B2)*(bcoef-B1p1)
   *
   * #least squares
   * [A1 A2] * [p1]  = tmpval
   *           [p2]  
   * (A1-A2*inv(B2)*B1)*p1 = tmpval - A2*inv(B2)*bcoef
   * p1 = pinv(A1-A2*inv(B2)*B1)*(tmpval - A2*inv(B2)*bcoef)
   *
   * @param NumMatrix&       B1       stores the constraint equations wrt other control points
   * @param NumMatrix&       B2       stores the constraint equations wrt fixed control points
   * @param NumMatrix&       A1       stores least squares equations wrt other control points        
   * @param NumMatrix&       A2       stores least squares equations wrt fixed control points
   * @param vector<Point3>&  bcoef    vector that stores control points for the boundary(ies)
   * @param vector<Point3>&  tmpval   limit poisition from CC surface
   * @param vector<Point3>&  p1       other control points
   * @param vector<Point3>&  p2       control points that satisfy constraints exactly 
   */
  void computeSpIndepControls(NumMatrix& B1, NumMatrix& B2, NumMatrix& A1, NumMatrix& A2, 
			      vector<Point3> &bcoef,vector<Point3> &tmpval, 
			      vector<Point3>& p1, vector<Point3>& p2);
   
  /**
   * Degree elevation routine from a cubic Bezier curve to any order Bezier Curve
   * @param int               deg        degree to elevate to
   * @param vector<Point3>&   cubic      vector that stores the cubic bezier controls
   * @param vector<Point3>&   high       vector to store the high order bezier controls
   */  
  void elevateDegree(int deg, vector<Point3>& cubic, vector<Point3> & high);
  
  /**
   * Routine to switch from a BSpline representation of a curve to Bezier representation 
   * @param vector<Point3>&   cubic      vector that stores the input bspline rep
   *                                     and has the resulting bezier rep at the end
   */  
  void switchToBezier(vector<Point3>&cubic);

  /**
   * Routine to switch from a Bezier representation of a curve to BSpline representation 
   * @param vector<Point3>&   high       vector that stores the input bezier rep and
   *                                     has the resulting bspline rep at the end
   */  
  void switchToBSpline(vector<Point3> &high);

  class DenseSampleData {
      public:
          double _init;
          double _step;
          int _num_samples;
          BolVector _can_interpolate_at_sample;
          NumMatrix _xy_values;
          NumMatrix _sample_positions;
          NumMatrix _sample_du;
          NumMatrix _sample_dv;
          double _bary_weights_x[4];
          //NumVector _bary_weights_x;
          //NumMatrix bary_weightx_x;
          
          bool can_interpolate(Point2 xy);
          Index2 to_index(Point2 xy);
          void get_interpolation_data(int i, Point2 xy, NumMatrix& interpolation_nodes, NumMatrix& function_values);
          DenseSampleData& operator=(const DenseSampleData& d){
            _init = d._init;
            _step = d._step;
            _num_samples = d._num_samples;
            _can_interpolate_at_sample = d._can_interpolate_at_sample;
            _xy_values       = d._xy_values       ;
            _sample_positions= d._sample_positions;
            _sample_du       = d._sample_du       ;
            _sample_dv       = d._sample_dv       ;
            memcpy(_bary_weights_x,
                    d._bary_weights_x,
                    sizeof(double)*4);
            //_bary_weights_x  = d._bary_weights_x  ;
            return *this;
          }

  };
  vector<DenseSampleData*> _dense_samples_on_patch;
  void setup_interpolation_samples(int V, double spacing);


};

extern int LoadW;


#endif

