#include "bdsurf.hpp"
#include "ccsubmatlib.hpp"
#include "ccsurf.hpp"
#include "ccsurfop.hpp"
#include "vecmatop.hpp"
#include "vec2t.hpp"
#include "mat2t.hpp"
#include "mat3t.hpp"
#include <iomanip>
#include "evalfromw.hpp"
#include "utils.hpp"
/*** DZ ***/

// global variables added not to mess with the main data structures to 
// implement loading W matrices mapping sample points to spline control points

int LoadW = 0;  // use W matrices loaded from files (default = no)

// convert gpmesh flags to chart ids used in evalfromW
int flags2chartId(int boun_flag, int corner_flag) { 
  if (boun_flag != GpMesh::BOUNDARY_VERTEX) return 0; 
  else if (corner_flag == GpMesh::CREASE_VERTEX) return 1; 
  else if (corner_flag == GpMesh::CONVEX_VERTEX) return 2; 
  else if (corner_flag == GpMesh::CONCAVE_VERTEX) return 3; 
  else { cerr << "bad corner" << endl;} return -1; 
}

/**** end DZ ***/

// ---------------------------------------------------------------------- 
BdSurf::~BdSurf(){
  for(int V=0; V<numVertices(); V++) {
    Vector<Point3>* cur = _controls[V];
    if(cur!=NULL) {
      delete cur;
    } 
    else { }
  }
}
// ---------------------------------------------------------------------- 
int BdSurf::setParams(int _iniSubdivCtrlLevel, int _iniBlendFuncType, 
		      int _iniChartType, int _iniBasisType, int _iniPreprocess, 
                      double _iniEvalLowerBnd, double _iniEvalUpperBnd, 
		      int _iniIndepBoundaryFit,  double _iniFlatS, 
		      int _inimax_degree_valence, int _iniSplineDeg, 
		      int _iniSplineBlendFuncDeg,double interpolation_spacing,
              int evaluate_via_interpolant,
              int refinement_factor,
               int _loadW){
    _interpolation_spacing = interpolation_spacing;
    _evaluate_via_interpolant = evaluate_via_interpolant;
  _subdivCtrlLevel = _iniSubdivCtrlLevel;
  assert(_iniSubdivCtrlLevel>=0 && _iniSubdivCtrlLevel<=2);
  
  _blendFuncType = _iniBlendFuncType;   //analytic vs b-spline partition of unity
  assert(_iniBlendFuncType==ANALYTIC_BLEND_FUNC || _iniBlendFuncType==BSPLINE_BLEND_FUNC || _iniBlendFuncType == OPTIM_BLEND_FUNC); 

  //BASIC
  _chartType = _iniChartType;
  assert(_iniChartType>=0 && _iniChartType<=2);
    
  _basisType = _iniBasisType;
  assert(_iniBasisType>=0 && _iniBasisType<=1);
  
  //SETUP
  _preprocess = (bool)_iniPreprocess;
  assert(_iniPreprocess==0 || _iniPreprocess==1);
  
  _evalLowerBnd = _iniEvalLowerBnd; 
  _evalUpperBnd = _iniEvalUpperBnd; 
  assert(_evalLowerBnd+_evalUpperBnd == 1);

  _max_degree_valence = _inimax_degree_valence;
  assert(_inimax_degree_valence == 0 || _inimax_degree_valence >=12);       
  //to avoid high degree polynomials, new deg=_max_degree_valence+2
  
  //CORNER RELATED
  //when set to 1 boundary depends on interior points
  _indepBoundaryFit = (bool)_iniIndepBoundaryFit; 
  assert(_indepBoundaryFit == 0 || _indepBoundaryFit == 1);
  _flatS = _iniFlatS;           //flatness parameter used in concave corners
  assert(_flatS == -1 || (_flatS>=0 && _flatS <= 1));
  //SPLINE FIT RELATED
  _splineDeg = _iniSplineDeg;   //order of spline curve
  assert(_splineDeg >=3);
  _splineBlendFuncDeg = _iniSplineBlendFuncDeg;    //degree of the spline function used as partition of unity
  //assert(_splineBlendFuncDeg >=3);

  _refinement_factor = refinement_factor;
  /**** DZ ***/
  ::LoadW = _loadW;
  /**** end DZ ***/

  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::setup(istream& in, int flipNormal, Point3 scale, void* _iniSubdivMatrices, void* _inimatricesU){

  /*** DZ ***/
  if(::LoadW) loadAllMatrices();
  /**** end DZ ***/
  //----------------------------------------------------------------------
  //0. set _subdivMatrices
  _subdivMatrices = (CCSubMatLib*)_iniSubdivMatrices;
  _matricesU = (BdULib*)_inimatricesU;
  

  //----------------------------------------------------------------------
  //1. _ccSurf preparation based on the subdivision surface
  _ccSurf.setup(_flatS, in,flipNormal,scale,_iniSubdivMatrices);
  assert(_ccSurf.inlvl()==0); //LEXING: no multires version
  for(int i=0; i<5; i++) {
    _ccSurf.subdivide(i, *_subdivMatrices);
  }
  //cerr<<"writing file"<<endl;
  //ofstream os("test.wrl");
  //_ccSurf.dump(os, 1);
  
  //----------------------------------------------------------------------
  //2. construction of _controls
  int numv = gpmesh().numVertices();
  constructControls(numv);
  
  //----------------------------------------------------------------------
  // loop over all vertices
  for(int V=0; V<numv; V++) {
    int K = valence(V);
    if(K==0) continue; //unused vertices
    
    vector<Point2> tmpcd;  //cd poistion
    vector<Point2> tmpxy;  //xy position   
    vector<Point3> tmpval; //surface position
    
    //----------------------------------------------------------------------
    //3.gather
    gatherData(V, tmpcd, tmpxy, tmpval);
    
    /*** DZ ***/
    //        if( valence(V) == 7 && gpmesh().boun(V) == GpMesh::BOUNDARY_VERTEX 
    //    	&& gpmesh().corner(V) == GpMesh::CREASE_VERTEX) {
    //          cout << "param locs" << valence(V) << endl;
    //          for(int i = 0; i < tmpxy.size(); i++)
    //	    cout << "[" <<tmpxy[i](0) << "," << tmpxy[i](1) << "]," << endl;
    //        }
    /*** end DZ ***/
    //----------------------------------------------------------------------
    //preprocessing
    if(_preprocess) { // only possible for closed surfaces for now.
      preProcessPoly(K, tmpval);
    }
    
    //----------------------------------------------------------------------
    //3.3. solve least square system
    if(       _basisType == POLYNOMIAL_BASIS) { //POLYNOMIAL
      //--------------------------------
      constructPoly(V, tmpxy, tmpval);
    }
    else if (_basisType == SPLINE_BASIS) {//SPLINE CONSTRUCTION
      constructSpline(V, tmpxy, tmpval);
    }
    else {
      assert(0);
    }
  }
    
  _dense_samples_on_patch.resize(numv);
  if(_evaluate_via_interpolant){
#pragma omp parallel for
  for(int V=0; V<numv; V++) {
      setup_interpolation_samples(V, _interpolation_spacing);
  } 
  }
  return 0;
}
void BdSurf::setup_interpolation_samples(int V, double spacing){
    DenseSampleData* sample_data = new DenseSampleData();
    
    // Estimate number of samples based on provided spacing
    double param_bound;
    paramBound(V, EVAL_UB(), param_bound); // max allowed xy values
    
    // compute position and derivatives
    Point2 xy(0.);
    vector<Point3> position_and_derivs(3,Point3());
    eval(EVAL_VALUE|EVAL_1ST_DERIV, V, xy.array(), position_and_derivs.data(), true);

    //mystery calculation
    double temp = spacing/max(position_and_derivs[1].l2(), position_and_derivs[2].l2()); 
    //cout << spacing << ", " << max(position_and_derivs[1].l2(), position_and_derivs[2].l2()) << endl;
    //cout << position_and_derivs[1] << ", "<< position_and_derivs[2] << endl;
    int tnum = int(ceil(param_bound/temp));
    //cout << temp << ", " << tnum << endl;
    double step_size= param_bound/tnum;
    //cout << "V: " << V << ", step: " << step_size << endl;
    double min_param_value= (-tnum+1)*step_size;
    double num_samples = 2*tnum+3;
    cout << "num_samples: " << num_samples << endl;
    
    NumMatrix& sample_positions = sample_data->_sample_positions;
    NumMatrix& xy_values= sample_data->_xy_values;
    NumMatrix& sample_du = sample_data->_sample_du;
    NumMatrix& sample_dv = sample_data->_sample_dv;
    BolVector& can_interpolate_at_sample= sample_data->_can_interpolate_at_sample;
    
    xy_values.resize(2, num_samples*num_samples);
    sample_positions.resize(3, num_samples*num_samples);
    sample_du.resize(3, num_samples*num_samples);
    sample_dv.resize(3, num_samples*num_samples);
    can_interpolate_at_sample.resize((num_samples-1)*(num_samples-1));
    clear(can_interpolate_at_sample);
    clear(xy_values);
    clear(sample_positions);
    clear(sample_du);
    clear(sample_dv);
    //setvalue(can_interpolate_at_sample, true);

    BolVector samples_in_interior(num_samples*num_samples);
    clear(samples_in_interior);
    
    Point2 init(min_param_value, min_param_value);
//#pragma omp parallel for
    for(int i=0; i < num_samples; i++){
        for(int j=0; j < num_samples; j++){
            int index = num_samples*i + j;

            // uniformly sample patch in xy-domain
            xy = init + step_size*Point2(i,j);
            
            // transform to cd coordinates
            Point2 cd; int f;
            Vxy2Vfcd(EVAL_VALUE, V, xy.array(), f, cd.array());

            // if the sample is in the patch interior, compute and save its
            // position on surface
            bool is_sample_in_interior =  cd.x()<.95*EVAL_UB() && cd.y()<.95*EVAL_UB();
            // TODO use memcpy instead
            if(is_sample_in_interior) { 
                vector<Point3> position_and_derivs(3,Point3());
                eval(EVAL_VALUE|EVAL_1ST_DERIV, V, xy.array(), position_and_derivs.data(), true);
                for(int d = 0; d < 2; d++)
                    xy_values(d,index) = xy(d);

                for(int d = 0; d < 3; d++){
                    sample_positions(d,index) = position_and_derivs[0](d);
                    sample_du(d,index) = position_and_derivs[1](d);
                    sample_dv(d,index) = position_and_derivs[2](d);
                }

            }
            samples_in_interior(index) = is_sample_in_interior;
        }
    }

    // Check if a 4x4 neighborhood around each sample is in the interior; if so,
    // we can interpolate the surface here.
    int num_interp = 0; 
    for(int i=1; i < num_samples-2; i++){
        for(int j=1; j < num_samples-2; j++){ // for each sample computed above...
            int idx = num_samples*i + j;
            bool can_interp = true;
            for(int pi=i-1; pi < i+3; pi++){ // for a 4x4 patch around the sample...
                for(int pj=j-1; pj < j+3; pj++){
                    int pidx = num_samples*pi + pj;
                    //cout << can_interpolate_at_sample(idx) << endl;
                    can_interp = can_interp &&
                         samples_in_interior(pidx);
                }
            }
            can_interpolate_at_sample(idx) = can_interp;
            num_interp += can_interp ? 1 : 0;
        }
    }
    cerr << "ninterp: " << num_interp << ", numsample: " << num_samples << endl;

    sample_data->_init = min_param_value;
    sample_data->_step = step_size;
    sample_data->_num_samples = num_samples;
    /*cout << "setup" << endl;
cout << sample_data->_init  << ", " 
 << sample_data->_step  << endl;*/
    int idx = 0;
    int it =0;
    for(int i =0; i < can_interpolate_at_sample.length(); i++){
        if(can_interpolate_at_sample(i) && it > 20){
            idx = i; it++;break;
        }
    }
    NumMatrix interp_nodes(2,16);
    NumMatrix func_values(3,16);
    sample_data->get_interpolation_data(0, 
            Point2(xy_values.clmdata(idx)), interp_nodes, func_values);

    NumVector temp_mat(4);
    for (int i = 0; i < 4; i++) {
       temp_mat(i) = interp_nodes(1,i); 
    }
    NumVector temp_weights =  bary_weights_1d(temp_mat);
    for (int i = 0; i < 4; i++)
        sample_data->_bary_weights_x[i] = temp_weights(i);
   _dense_samples_on_patch[V] = sample_data; 
}

Index2 BdSurf::DenseSampleData::to_index(Point2 xy){
    /*double xya[] = {xy.x(), xy.y()};
    double hashed_xy_arr[] = {(xya[0] - _init)/_step,
                            (xya[1] - _init)/_step};
    int index_xy_arr[] = {int(floor(hashed_xy_arr[0])), int(floor(hashed_xy_arr[1]))};
    Index2 array_based(
            min(max(index_xy_arr[0],1),_num_samples-3),
            min(max(index_xy_arr[1],1),_num_samples-3)
            );
    return array_based;*/
    return Index2(
            min(max(int(floor((xy.x() - _init)/_step)), 1), _num_samples - 3),
            min(max(int(floor((xy.y() - _init)/_step)), 1), _num_samples - 3)

            );
}

bool BdSurf::DenseSampleData::can_interpolate(Point2 xy){
    Index2 index_xy = to_index(xy);
    //cout << "xy: " << xy << endl;
    //cout << "idx: " << index_xy << endl;
    //cout << "interp ok: " << _can_interpolate_at_sample(index_xy.x()*_num_samples + index_xy.y()) << endl;
    return _can_interpolate_at_sample(index_xy.x()*_num_samples + index_xy.y());
}
void BdSurf::DenseSampleData::get_interpolation_data(int data_index, Point2 xy, 
        NumMatrix& interpolation_nodes, NumMatrix& function_values){
    assert(interpolation_nodes.n() == function_values.n());
    assert(interpolation_nodes.n() == 16);
    assert(interpolation_nodes.m() == 2);
    assert(function_values.m() == 3);
    
    Index2 index_xy = to_index(xy);
    //assert(index_xy.x() >= 1 && index_xy.x() <= _num_samples-2);
    //assert(index_xy.y() >= 1 && index_xy.y() <= _num_samples-2);
    // TODO FIGURE OUT WHY THIS FAILS
    assert(can_interpolate(xy));
    int i = index_xy.x();
    int j = index_xy.y();
    for(int pi=0; pi < 4; pi++){ // for a 4x4 patch around the sample...
        //for(int pj=0; pj < 4; pj++){
        int pj = 0;
        int pidx = 4*pi + pj;
        int index = _num_samples*(i-1+pi) + (j-1+pj);
        memcpy(
                interpolation_nodes.clmdata(pidx),_xy_values.clmdata(index), 
                sizeof(double)*2*4);
        switch(data_index){
            case 0:
                memcpy(
                        function_values.clmdata(pidx),_sample_positions.clmdata(index), 
                        sizeof(double)*3*4);
                break;
            case 1:
                memcpy(
                        function_values.clmdata(pidx),_sample_du.clmdata(index), 
                        sizeof(double)*3*4);
                break;
            case 2:
                memcpy(
                        function_values.clmdata(pidx),_sample_dv.clmdata(index), 
                        sizeof(double)*3*4);
                break;
            default:
                assert(0); 
        }

        //}
    }
    /*for(int pi=0; pi < 4; pi++){ // for a 4x4 patch around the sample...
        for(int pj=0; pj < 4; pj++){
            //int pj = 0;
            int pidx = 4*pi + pj;
            int index = _num_samples*(i-1+pi) + (j-1+pj);
            
            for (int d = 0; d < 2; d++) {
               interpolation_nodes(d,pidx) = _xy_values(d,index); 
            }

            for (int d = 0; d < 3; d++) {
                switch(data_index){
                    case 0:
                        function_values(d,pidx) = _sample_positions(d,index); 
                        break;
                    case 1:
                        function_values(d,pidx) = _sample_du(d,index); 
                        break;
                    case 2:
                        function_values(d,pidx) = _sample_dv(d,index); 
                        break;
                    default:
                        assert(0); 
                }
            }
            
        }
    }*/
}
// ----------------------------------------------------------------------  
int BdSurf::rotate(const Quaternion& qr){
  
  //1. ccsurf
  _ccSurf.rotate(qr);
  //2. cgr
  NumMatrix tmp(3,3);  qr.toRotationMatrix(tmp);
  Matrix3 rot(tmp.data());
  for(int V=0; V<numVertices(); V++){
    Vector<Point3>* cur = _controls[V];
    if(cur!=NULL) {
      Vector<Point3>& tmp = *cur;
      for(int i=0; i<tmp.length(); i++)
	tmp(i) = rot * tmp(i);
    } 
        else { }
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::shift(const Point3& sh){
  
  //1. ccsurf
  _ccSurf.shift(sh);
    //2. cgr
  for(int V=0; V<numVertices(); V++) {
    Vector<Point3>* cur = _controls[V];
    if(cur!=NULL) {
      Vector<Point3>& tmp = *cur;
      tmp(0) = tmp(0) + sh;
    } 
    else { }
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::paramBound(int V, double ub, double& rad){
  
  double bbd = 0;
  int K = valence(V);
  NumVector D;
  double p = 0; double q = 0;  //constants

  if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX) {//not boundary    
    D = _subdivMatrices->chooseMatrix(CCSubMatLib::INTERIOR, K).D();
    if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
      p = 4.0/K;	 q = 4.0/K;
    } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
      p = log(1.0/D(1))/log(2.0);	 q = 4.0/K;
    } else if(_chartType == ISODISTANCE_CHART) {
      p = 1.0;	 q = 4.0/K;
    }
  } 
  
  else if(gpmesh().corner(V) == GpMesh::CREASE_VERTEX){
    D = _subdivMatrices->chooseMatrix(CCSubMatLib::BOUNDARY, K).D();
    if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
      p = 2.0/K;	 q = 2.0/K;
    } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
      p = log(1.0/D(1))/log(2.0);	 q = 2.0/K;
    } else if(_chartType == ISODISTANCE_CHART) {
      p = 1.0;	 q = 2.0/K;
    }
  }
  else if(gpmesh().corner(V) == GpMesh::CONVEX_VERTEX){
    D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONVEX, K).D();
    if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
      p = 1.0/K;	 q = 1.0/K;
    } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
      p = log(1.0/D(1))/log(2.0);	 q = 1.0/K;
    } else if(_chartType == ISODISTANCE_CHART) {
      p = 1.0;	 q = 1.0/K;
    }
  }
  else if(gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX){
    D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).D();
    if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
      p = 3.0/K;	 q = 3.0/K;
    } else {//CHARACTERISTIC = ISODISTANCE
      p = 1.0;	 q = 3.0/K;
    }
  }
  
  
  
  //evaluate
  double r = sqrt(2.0) * ub;   double a = PI/4.0;
  double rp = pow(r,p);        double aq = a*q;
  double xn = rp*cos(aq);      double yn = rp*sin(aq);
 
  double ang = 0.0;

  for(int f=0; f<K; f++) {
    if(gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX)
      ang = (2*PI*f)/K;
    else if(gpmesh().corner(V) == GpMesh::CREASE_VERTEX)
      ang = (PI*f)/K;
    else if(gpmesh().corner(V) == GpMesh::CONVEX_VERTEX)
      ang = ((1/2)*PI*f)/K;
    else if(gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX)
      ang = ((3/2)*PI*f)/K;
    
    double xt = xn*cos(ang)-yn*sin(ang);            double yt = xn*sin(ang)+yn*cos(ang);
    bbd = max(max(bbd,abs(xt)),abs(yt));
  }
  rad = bbd;
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::eval(int flags, int V, double* xy, Point3* ret, bool skip_interp){
    assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
    DenseSampleData* sample_data = _dense_samples_on_patch[V];
    if(_evaluate_via_interpolant && !skip_interp  &&
            ( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) ) &&
           sample_data->can_interpolate(Point2(xy)) ){
            Index2 index_xy = sample_data->to_index(xy);

            Point2 eval_point(xy);
            //NumVector& bary_weights = sample_data->_bary_weights_x;
            //double* bary_weights_arr = bary_weights.data();
            double* bary_weights_arr = sample_data->_bary_weights_x;

            if(flags & EVAL_VALUE){
                bary2d_2(index_xy, 3, sample_data->_xy_values, sample_data->_sample_positions, 
                        bary_weights_arr,  eval_point, ret[0]);
            }
            if(flags & EVAL_1ST_DERIV){
                // interpolate du
                bary2d_2(index_xy, 3, sample_data->_xy_values, sample_data->_sample_du,  
                        bary_weights_arr,  eval_point, ret[1]);

                // interpolate dv 
                bary2d_2(index_xy, 3, sample_data->_xy_values, sample_data->_sample_dv, 
                        bary_weights_arr,  eval_point, ret[2]);
                
                /*if(ret[0].l2() > 20 || ret[1].l2() > 20 ||ret[2].l2() > 20){
                    Point2 cd; int f;
                    Vxy2Vfcd(EVAL_VALUE, V, xy, f, cd.array());
                    cout << "xy: " << Point2(xy) << "cd: " << cd <<", f: "<< f << endl;
                }*/
            }

    } else {

        double reprLB = _evalLowerBnd;
        double reprUB = _evalUpperBnd;
        int Vi = V;
        double* xyi = xy;
        int fi;
        double cdi[12];   Vxy2Vfcd(flags, Vi, xyi, fi, cdi);    

        //LEXING: IMPORTANT check valid range for evaluation
        //assert(abs(cdi[0])<=reprUB+SCL_EPS && abs(cdi[1])<=reprUB+SCL_EPS);
        if(abs(cdi[0])<=reprLB && abs(cdi[1])<=reprLB)  // ONE CHART
            chartEval(flags, Vi, xyi, ret,0);
        else { //MULTIPLE CHARTS
            int F;
            double st[12]; Vfcd2Fcd(flags, Vi, fi, cdi, F, st);

            //left weight
            double lb[6];   lb[0] = 1.;	 
            for(int g=1; g<6; g++) lb[g] = 0.;	 

            //accumulation of results
            Point3 ap[6];	 
            for(int i=0; i<6; i++)  ap[i]=Point3(0,0,0);	      


            for(int v=0; v<4; v++) {
                //get cdj
                int Vj;
                int fj;
                double cdj[12];      Fcd2Vfcd(flags, F, st, v, Vj, fj, cdj);      

                if( abs(cdj[0])<reprUB && abs(cdj[1])<reprUB) {
                    //within evaluation range 
                    double xyj[12]; Vfcd2Vxy(flags, Vj, fj, cdj, xyj);         
                    //evaluate bj wrt cdj
                    double bj[6];   blendFuncEval(flags, cdj, reprLB, reprUB, bj);       
                    //evaluate pj wrt xyj
                    Point3 pj[6];   chartEval(flags, Vj, xyj, pj,0);   
                    //evaluate pj wrt cdj
                    Point3 pjn[6];  compose(flags, 3, (double*)pj, xyj, (double*)pjn); 
                    //evaluate pj wrt st
                    Point3 pjs[6];  compose(flags, 3, (double*)pjn, cdj, (double*)pjs);      
                    //evaluate bj wrt st
                    double bjs[6];  compose(flags, 1, (double*)bj,  cdj, (double*)bjs);  

                    //combine bjs+pjs into ap and lb
                    if(flags & EVAL_VALUE) {
                        ap[0] += bjs[0] * pjs[0];
                        lb[0] -= bjs[0];
                    }
                    if(flags & EVAL_1ST_DERIV) {
                        ap[1] += bjs[1] * pjs[0] + bjs[0] * pjs[1];
                        ap[2] += bjs[2] * pjs[0] + bjs[0] * pjs[2];
                        lb[1] -= bjs[1];
                        lb[2] -= bjs[2];
                    }
                    if(flags & EVAL_2ND_DERIV) {
                        ap[3] += bjs[3] * pjs[0] + 2*bjs[1]*pjs[1] + bjs[0] * pjs[3];
                        ap[4] += bjs[4] * pjs[0] + bjs[2]*pjs[1] + bjs[1]*pjs[2] + bjs[0] * pjs[4];
                        ap[5] += bjs[5] * pjs[0] + 2*bjs[2]*pjs[2] + bjs[0] * pjs[5];
                        lb[3] -= bjs[3];
                        lb[4] -= bjs[4];
                        lb[5] -= bjs[5];
                    }
                }
            }
            assert(lb[0] < 1e-10);
            Point3 apn[6];
            compose(flags, 3, (double*)ap, st, (double*)apn);
            compose(flags, 3, (double*)apn, cdi, (double*)ret);  
        }
        return 0;
    }
}



// ---------------------------------------------------------------------- 
int BdSurf::eval2(int flags, int V, double* xy, Point3* ret){
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );

  double reprLB = _evalLowerBnd;
  double reprUB = _evalUpperBnd;
  int Vi = V;
  double* xyi = xy;
  int fi;
  /*       
      NumMatrix& IV = _subdivMatrices->chooseMatrix(CCSubMatLib::INTERIOR, 10).IV();
     cout<<1.0/IV(0,0)<<endl;
     exit(1);
  */
  double cdi[12];   Vxy2Vfcd(flags, Vi, xyi, fi, cdi);    
//  if(V == 12) cout<<"cdi = "<<cdi[0]<<" "<<cdi[1]<<endl;
  //LEXING: IMPORTANT check valid range for evaluation
//  assert(abs(cdi[0])<=reprUB+SCL_EPS && abs(cdi[1])<=reprUB+SCL_EPS);
  if(abs(cdi[0])<=reprLB && abs(cdi[1])<=reprLB)  // ONE CHART
    chartEval(flags, Vi, xyi, ret,0);
  else { //MULTIPLE CHARTS
    int F;
    double st[12]; Vfcd2Fcd(flags, Vi, fi, cdi, F, st);
    //  if(V==12) cout<<"st = "<<st[0]<<" "<<st[1]<<endl;
    //left weight
    double lb[6];   lb[0] = 1.;	 
    for(int g=1; g<6; g++) lb[g] = 0.;	 
    
    //accumulation of results
    Point3 ap[6];	 
    for(int i=0; i<6; i++)  ap[i]=Point3(0,0,0);	      
    
    for(int v=0; v<4; v++) {
        
        //get cdj
        int Vj;
        int fj;
        double cdj[12];      Fcd2Vfcd(flags, F, st, v, Vj, fj, cdj);      
        //  if(V==12)
        //{ cout<<"Vj = "<<Vj<<endl;
        //cout<<"cdj = "<<cdj[0]<<" "<<cdj[1]<<endl;
        // }
        //      if(Vj == 12){
        
        if( abs(cdj[0])<reprUB && abs(cdj[1])<reprUB) {
            //within evaluation range 
            double xyj[12]; Vfcd2Vxy(flags, Vj, fj, cdj, xyj);         
            //evaluate bj wrt cdj
            double bj[6];   blendFuncEval(flags, cdj, reprLB, reprUB, bj);       
            //double bj[6]; bj[0] = 1.; bj[1]=bj[2] =bj[3] =bj[4]=bj[5] =0.;
            //evaluate pj wrt xyj
            Point3 pj[6];   chartEval(flags, Vj, xyj, pj,0);   
            //assert(pj[0](2) == 1);
//evaluate pj wrt cdj
            Point3 pjn[6];  compose(flags, 3, (double*)pj, xyj, (double*)pjn); 
            //evaluate pj wrt st
            Point3 pjs[6];  compose(flags, 3, (double*)pjn, cdj, (double*)pjs);      
            //evaluate bj wrt st
            double bjs[6];  compose(flags, 1, (double*)bj,  cdj, (double*)bjs);  
            
/*
  if(V==12){
  cout<<"bj = "<<bj[0]<<" "<<bj[1]<<" "<<bj[2]<<endl;
  cout<<"pj = "<<pj[0]<<endl;
  cout<<pj[1]<<endl;
  cout<<pj[2]<<endl;
  cout<<"bjs = "<<bjs[0]<<" "<<bjs[1]<<" "<<bjs[2]<<endl;
  cout<<"pjs = "<<pjs[0]<<endl;
  cout<<pjs[1]<<endl;
  cout<<pjs[2]<<endl;
  exit(1);
  }
*/
            
            //combine bjs+pjs into ap and lb
            if(flags & EVAL_VALUE) {
                //if(V==12)
                
                //cout<<"in val = "<<bjs[0]*pjs[0]<<endl;
                ap[0] += bjs[0] * pjs[0];
                lb[0] -= bjs[0];
            }
            if(flags & EVAL_1ST_DERIV) {
                //if(V==12){
                //cout<<"in v1 = "<<bjs[1] * pjs[0] + bjs[0] * pjs[1]<<endl;
                //cout<<"in v2 = "<<bjs[2] * pjs[0] + bjs[0] * pjs[2]<<endl;
                // }
                ap[1] += bjs[1] * pjs[0] + bjs[0] * pjs[1];
                ap[2] += bjs[2] * pjs[0] + bjs[0] * pjs[2];
                lb[1] -= bjs[1];
                lb[2] -= bjs[2];
            }
            if(flags & EVAL_2ND_DERIV) {
                ap[3] += bjs[3] * pjs[0] + 2*bjs[1]*pjs[1] + bjs[0] * pjs[3];
                ap[4] += bjs[4] * pjs[0] + bjs[2]*pjs[1] + bjs[1]*pjs[2] + bjs[0] * pjs[4];
                ap[5] += bjs[5] * pjs[0] + 2*bjs[2]*pjs[2] + bjs[0] * pjs[5];
                lb[3] -= bjs[3];
                lb[4] -= bjs[4];
                lb[5] -= bjs[5];
            }
            
//                if(V==12)
            //cout<<"v = "<<v<<" before comp1 = "<<ap[0]<<" "<<ap[1]<<" "<<ap[2]<<endl;
            // }
        }
    }
    Point3 apn[6];
    compose(flags, 3, (double*)ap, st, (double*)apn);
    //if(V==12)
    // cout<<"after comp1 = "<<ret[0]<<" "<<ret[1]<<" "<<ret[2]<<endl;
    compose(flags, 3, (double*)apn, cdi, (double*)ret);  
    //if(V==12)
    //cout<<"after comp2 = "<<ret[0]<<" "<<ret[1]<<" "<<ret[2]<<endl;
    
  }
//  if(V==12)
//      exit(1);
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::singleChartEval(int flags, int V, double* xy, Point3* ret){
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
  int Vi = V;
  double* xyi = xy;
  int fi;
  double cdi[12];  Vxy2Vfcd(flags, Vi, xyi, fi, cdi);  
  double TMPEPS;
  TMPEPS = 1e-5;  
  
//  assert(cdi[0]<=1+TMPEPS && cdi[1]<=1+TMPEPS);  
  chartEval(flags, Vi, xyi, ret,0);
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::dump(ostream& fout){
  for(int V=0; V<numVertices(); V++) {
    fout<<"VERTEX "<<V<<endl;
    Vector<Point3>* cur = _controls[V];
    if(cur!=NULL) {
      Vector<Point3>& tmp = *cur;
      fout<<tmp<<endl;
    }
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::compose(int flags, int dof, double* U, double* C, double* R){
  //u--u(c,d), c--c(x,y), r--u(x,y)
  double* u  =&(U[0]);  
  double* uc =&(U[dof]);    double* ud =&(U[2*dof]);  
  double* ucc=&(U[3*dof]);  double* ucd=&(U[4*dof]);  double* udd=&(U[5*dof]);
  
  //  double& c  = C[0];  
  double& cx = C[2];  double& cy = C[4];  
  double& cxx= C[6];  double& cxy= C[8];  double& cyy= C[10];
  
  //  double& d  = C[1];  
  double& dx = C[3];  double& dy = C[5];  
  double& dxx= C[7];  double& dxy= C[9];  double& dyy= C[11];
  
  double* r  =&(R[0]);  
  double* ux =&(R[dof]);    double* uy =&(R[2*dof]);  
  double* uxx=&(R[3*dof]);  double* uxy=&(R[4*dof]);  
  double* uyy=&(R[5*dof]);
  
  //0 order
  if(flags & EVAL_VALUE) {
    for(int d=0; d<dof; d++) r[d] = u[d];
  }
  if(flags & EVAL_1ST_DERIV) {
    for(int i=0; i<dof; i++) {
      ux[i] = uc[i]*cx + ud[i]*dx;
      uy[i] = uc[i]*cy + ud[i]*dy;
    }
  }
  if(flags & EVAL_2ND_DERIV) {
    for(int i=0; i<dof; i++) {
      uxx[i] = ucc[i]*cx*cx + 2*ucd[i]*cx*dx + udd[i]*dx*dx       + uc[i]*cxx + ud[i]*dxx;
      uxy[i] = ucc[i]*cx*cy + ucd[i]*(cx*dy+cy*dx) + udd[i]*dx*dy + uc[i]*cxy + ud[i]*dxy;
      uyy[i] = ucc[i]*cy*cy + 2*ucd[i]*cy*dy + udd[i]*dy*dy       + uc[i]*cyy + ud[i]*dyy;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::Vxy2Vfcd(int flags, int V, double* xy, int& f, double* cd){
    assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );

    int K = valence(V);    
    NumVector D;
    double p = 0; double q = 0;
    if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX){
      D = _subdivMatrices->chooseMatrix(CCSubMatLib::INTERIOR, K).D();
      if(_chartType == FULLY_COMPLEX_CHART) { //COMPLEX
	p = K/4.0;	 q = K/4.0;
      } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
	p =log(2.0)/log(1.0/D(1)); 	 q = K/4.0;
      } else if(_chartType == ISODISTANCE_CHART) { //ISODISTANCE
	p = 1.0;	 q = K/4.0;
      }
    }
    else{
      if (gpmesh().corner(V) == GpMesh::CREASE_VERTEX) {//crease vertex
	D = _subdivMatrices->chooseMatrix(CCSubMatLib::BOUNDARY, K).D();
	if(_chartType == FULLY_COMPLEX_CHART) { //COMPLEX
	  p = K/2.0;	 q = K/2.0;
	} else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
	  p = log(2.0)/log(1.0/D(1));	 q = K/2.0;
	} else if(_chartType == ISODISTANCE_CHART) { //ISODISTANCE
	  p = 1.0;	 q = K/2.0;
	}
      }
      else if (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX){ //convex vertex
	D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONVEX, K).D();
	if(_chartType == FULLY_COMPLEX_CHART) { //COMPLEX
	  p = K;	 q = K;
	} else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
	  p = log(2.0)/log(1.0/D(1));	 q = K;
	} else if(_chartType == ISODISTANCE_CHART) { //ISODISTANCE
	  p = 1.0;	 q = K;
	}
      }
      else if (gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX){ //concave vertex
	D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).D();
	if(_chartType == FULLY_COMPLEX_CHART) { //COMPLEX
	  p = K/3.0;	 q = K/3.0;
	}
	else {//CHARACTERISTIC OR ISODISTANCE
	  p = 1.0;	 q = K/3.0;
	}	
      } 
    }
    
    // constants
    //----------------------------
    if(abs(xy[0])<=SCL_EPS && abs(xy[1])<=SCL_EPS) {
      //----------------------------
      //map to face 0
      double* val = cd;	 double* ext = cd+2;	 double* fnl = cd+6;
      if(flags & EVAL_VALUE) {
	val[0] = 0; val[1] = 0;
      }
      if(flags & EVAL_1ST_DERIV) { //ONLY VALID FOR VALENCE 4 
	ext[0] = 1; ext[1] = 0;
	ext[2] = 0; ext[3] = 1;
      }
      if(flags & EVAL_2ND_DERIV) { //ONLY VALID FOR VALENCE 4
	for(int g=0; g<6; g++)
	  fnl[g] = 0;
      }
      f = 0;
    } 
    else {      
      double x = xy[0];  
      double y = xy[1];
      double ang = atan2(y,x);
      int i = 0; double R = 0;
      double overK = (1.0/double(K));
      double overPI = (1.0/PI);
      
      if(ang<0) ang=ang+2*PI;    
      if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX){
	i = min((int)(ang*K*.5*overPI), K-1); 
	R = -2*PI*i*overK;
      }
      else if (gpmesh().corner(V) == GpMesh::CREASE_VERTEX){
	i=min(int(ang*K*overPI), K-1);
	R = -1*PI*i*overK;
      }  
      else if (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX){
	i=min((int)(ang*K*2.0*overPI), K-1);
	R = -.5*PI*i*overK;
      }  
      else if (gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX){
	i=min((int)(ang*K*(2.0/3.0)*overPI), K-1);
	R = -1.5*PI*i*overK;
      } 
      
      //------------------
      double cR, sR;
      sR = sin(R);
      cR = cos(R);
//      sincos(R, &sR, &cR);
      double st2xy[12];
      if(flags & EVAL_VALUE) {
	st2xy[0] = cR*x - sR*y;  st2xy[1] = sR*x + cR*y;
      }
      if(flags & EVAL_1ST_DERIV) {     
	st2xy[2] = cR;   st2xy[3] = sR;  
	st2xy[4] = -sR;  st2xy[5] = cR;
      }
      if(flags & EVAL_2ND_DERIV) {
	double* tmp = st2xy+6; for(int g=0; g<6; g++) tmp[g]=0;
      }
      
      //------------------
      double ra2st[12];
      double s = st2xy[0];  double t = st2xy[1];
      double r2 = s*s+t*t;  double r = sqrt(r2);  
      double a = atan2(t,s);  
      double r3 = r2*r;  double r4 = r2*r2;
      
      if(flags & EVAL_VALUE) {
	ra2st[0] = r;  ra2st[1] = a;
      }
      if(flags & EVAL_1ST_DERIV) {
	ra2st[2] = s/r;  ra2st[3] = -t/r2;  
	ra2st[4] = t/r;  ra2st[5] = s/r2; 
      }
      if(flags & EVAL_2ND_DERIV) {
	ra2st[6] = 1/r-s*s/r3;  ra2st[7] = 2*t*s/r4;  ra2st[8] = -s*t/r3;
	ra2st[9] = (-s*s+t*t)/r4;  ra2st[10] = 1/r-t*t/r3;  ra2st[11] = -2*t*s/r4;
	
      }      
      //------------------
      double cd2ra[12];
      double rpm2 = exp((p-2)*log(r)); //pow(r,p-2); //
      double rpm1 = rpm2*r;  
      double rp = rpm1*r; 
      double aq = a*q;
      double caq, saq;
      saq = sin(aq);
      caq = cos(aq);
//sincos(aq, &saq, &caq);
      if(flags & EVAL_VALUE) {
	cd2ra[0] = rp*caq;  cd2ra[1] = rp*saq;
      }
      if(flags & EVAL_1ST_DERIV) {
	cd2ra[2] = rpm1*p*caq; cd2ra[3] = rpm1*p*saq;  
	cd2ra[4] = -q*cd2ra[1];  cd2ra[5] = cd2ra[0]*q;
      }
      if(flags & EVAL_2ND_DERIV) {
	cd2ra[6] = rpm2*p*caq*(p-1);  cd2ra[7] = rpm2*p*saq*(p-1);  
	cd2ra[8] = -rpm1*p*saq*q;     cd2ra[9] = rpm1*p*caq*q;  
	cd2ra[10] = -q*rp*cd2ra[5];       cd2ra[11] = cd2ra[4]*q;
      }
      double ra2xy[12];
      compose(flags, 2, ra2st, st2xy, ra2xy);
      compose(flags, 2, cd2ra, ra2xy, cd);   
      f = i; 
    }
    return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::Vfcd2Vxy(int flags, int V, int f, double* cd, double* xy){
    assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
    
    int K = valence(V);//valence
    NumVector D;
    double p = 0; double q = 0;  //constants
    if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX) {//not boundary    
      D = _subdivMatrices->chooseMatrix(CCSubMatLib::INTERIOR, K).D();
      if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 4.0/K;	 q = 4.0/K;
      } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
	p = log(1.0/D(1))/log(2.0);	 q = 4.0/K;
      } else if(_chartType == ISODISTANCE_CHART) {
	p = 1.0;	 q = 4.0/K;
      }
    } 
    
    else if(gpmesh().corner(V) == GpMesh::CREASE_VERTEX){
      D = _subdivMatrices->chooseMatrix(CCSubMatLib::BOUNDARY, K).D();
      if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 2.0/K;	 q = 2.0/K;
      } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
	p = log(1.0/D(1))/log(2.0);	 q = 2.0/K;
      } else if(_chartType == ISODISTANCE_CHART) {
	p = 1.0;	 q = 2.0/K;
      }
    }
    else if(gpmesh().corner(V) == GpMesh::CONVEX_VERTEX){
      D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONVEX, K).D();
      if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 1.0/K;	 q = 1.0/K;
      } else if(_chartType == CHARACTERISTIC_MAP_CHART) { //CHARACTERISTIC MAP
	p = log(1.0/D(1))/log(2.0);	 q = 1.0/K;
      } else if(_chartType == ISODISTANCE_CHART) {
	p = 1.0;	 q = 1.0/K;
      }
    }
    else if(gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX){
      D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).D();
      if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 3.0/K;	 q = 3.0/K;
      } else {//CHARACTERISTIC = ISODISTANCE
	p = 1.0;	 q = 3.0/K;
      }
    }
    
    
    //---------------------
    if( abs(cd[0])<=SCL_EPS && abs(cd[1])<=SCL_EPS ) {
      //-------------------	 //CANNOT TRANSFROM AT C==0 && D==0	 
      double* val = xy;	 
      if(flags & EVAL_VALUE) {
	val[0] = 0; val[1] = 0;
      }
      if(flags & EVAL_1ST_DERIV) {assert(0);}
      if(flags & EVAL_2ND_DERIV) {assert(0);}
    } 
    else {
      
      double c = cd[0]; double d = cd[1];
      //-------------------
      double ra2cd[12];
      double r2 = c*c+d*d;  double r = sqrt(r2);  
      double a = atan2(d,c);  
      double r3 = r2*r;  double r4 = r2*r2;
      
      if(flags & EVAL_VALUE) {
	ra2cd[0] = r;  ra2cd[1] = a;
      }
      if(flags & EVAL_1ST_DERIV) {
	ra2cd[2] = c/r;  ra2cd[3] = -d/r2;  
	ra2cd[4] = d/r;  ra2cd[5] = c/r2; 
      }
      if(flags & EVAL_2ND_DERIV) {
	ra2cd[6] = 1/r-c*c/r3;  ra2cd[7] = 2*d*c/r4;  
	ra2cd[8] = -c*d/r3;     ra2cd[9] = (-c*c+d*d)/r4;  
	ra2cd[10] = 1/r-d*d/r3; ra2cd[11] = -2*d*c/r4;
      }
      //-------------------
      double st2ra[12];      
      double rpm2= exp((p-2)*log(r)); // pow(r,p-2);
      double  rpm1 = rpm2*r;//pow(r,p-1);    
      double  rp = rpm1*r; //pow(r,p);     
      double aq = a*q;
      double caq, saq;
      saq = sin(aq); caq= cos(aq);
      //sincos(aq, &saq, &caq);
      
      if(flags & EVAL_VALUE)
	st2ra[0] = rp*caq;  st2ra[1] = rp*saq;
      
      if(flags & EVAL_1ST_DERIV) {
	st2ra[2] = rpm1*p*caq;  st2ra[3] = rpm1*p*saq;  
	st2ra[4] = -q*st2ra[1];   st2ra[5] = st2ra[0]*q;
      }
      if(flags & EVAL_2ND_DERIV) {
	st2ra[6] = rpm2*p*caq*(p-1);  st2ra[7]  = rpm2*p*saq*(p-1);  
	st2ra[8] = -rpm1*p*saq*q;     st2ra[9]  = rpm1*p*caq*q;  
	st2ra[10] = -q*rp*st2ra[5];   st2ra[11] = st2ra[4]*q;
      }
      
      //-------------------
      double xy2st[12];
      double R = 0;
      double foverk = (double(f)/double(K));
      if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX) //not boundary    
	R = 2*PI*foverk;
      else if(gpmesh().corner(V) == GpMesh::CREASE_VERTEX)
	R = PI*foverk;
      else if(gpmesh().corner(V) == GpMesh::CONVEX_VERTEX)
	R = (PI/2.0)*foverk;
      else if(gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX)
	R = (3*PI/2.0)*foverk;
      
      double cR, sR;
      sR = sin(R); cR = cos(R);
      //sincos(R, &sR, &cR);
      double s = st2ra[0];  double t = st2ra[1];
      if(flags & EVAL_VALUE) {
	xy2st[0] = cR*s - sR*t;  xy2st[1] = sR*s + cR*t;
      }
      if(flags & EVAL_1ST_DERIV) {
	
	xy2st[2] = cR;  xy2st[3] = sR;  
	xy2st[4] = -sR; xy2st[5] = cR;
      }
      if(flags & EVAL_2ND_DERIV) {
	double* tmp = xy2st+6; 
	for(int g=0; g<6; g++) tmp[g]=0;
      }
      double st2cd[12];
      compose(flags, 2, st2ra, ra2cd, st2cd);
      compose(flags, 2, xy2st, st2cd, xy);  
    }
    return 0;    
}
// ---------------------------------------------------------------------- 
int BdSurf::Vfcd2Fcd(int flags, int V, int f, double* cd, int& F, double* st){    
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) ); //LEXING
  
  intpair Fv = gpmesh().Vf2Fv(V,f);
  double* val = st;
  double* ext = st+2; //first order
  double* fnl = st+6; //second order
  F = Fv.first;
  switch(Fv.second) {
  case 0:
    if(flags & EVAL_VALUE) {
      val[0]=cd[0];  val[1]=cd[1];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = 1;    ext[1] = 0;	
      ext[2] = 0;    ext[3] = 1;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  case 1:
    if(flags & EVAL_VALUE) {
      val[0]=1-cd[1];	val[1]=cd[0];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = 0;	ext[1] = 1;
      ext[2] = -1;	ext[3] = 0;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  case 2:
    if(flags & EVAL_VALUE) {
      val[0]=1-cd[0];	val[1]=1-cd[1];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = -1;	ext[1] = 0;
      ext[2] = 0;	ext[3] = -1;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  case 3:
    if(flags & EVAL_VALUE) {
      val[0]=cd[1];	val[1]=1-cd[0];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = 0;	ext[1] = -1;
      ext[2] = 1;	ext[3] = 0;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  default:
    assert(0);
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::Fcd2Vfcd(int flags, int F, double* cd, int v, int& V, int& f, double* st){
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) ); 
  
  intpair Vf = gpmesh().Fv2Vf(F,v);
  //return val
  double* val = st;
  double* ext = st+2;
  double* fnl = st+6;
  V = Vf.first;
  f = Vf.second;
  switch(v) {
  case 0:
    if(flags & EVAL_VALUE) {
      val[0] = cd[0];	val[1] = cd[1];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = 1;	ext[1] = 0;
      ext[2] = 0;	ext[3] = 1;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  case 1:
    if(flags & EVAL_VALUE) {
      val[0] = cd[1];	val[1] = 1-cd[0];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = 0;	ext[1] = -1;
      ext[2] = 1;	ext[3] = 0;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  case 2:
    if(flags & EVAL_VALUE) {
      val[0] = 1-cd[0];	val[1] = 1-cd[1];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = -1;	ext[1] = 0;
      ext[2] = 0;	ext[3] = -1;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  case 3:
    if(flags & EVAL_VALUE) {
      val[0] = 1-cd[1];	val[1] = cd[0];
    }
    if(flags & EVAL_1ST_DERIV) {
      ext[0] = 0;	ext[1] = 1;
      ext[2] = -1;	ext[3] = 0;
    }
    if(flags & EVAL_2ND_DERIV) {
      for(int g=0; g<6; g++) fnl[g] = 0;
    }
    break;
  default:
    assert(0);
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::blendFuncEval(int flags, double* cd, double lb, double ub, double* pou){
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
  
  double resu[3];  blendFunc1D(flags, cd[0], lb, ub, resu);
  double resv[3];  blendFunc1D(flags, cd[1], lb, ub, resv);
  double* val = pou;
  double* ext = pou+1;
  double* fnl = pou+3;
  if(flags & EVAL_VALUE) {
    val[0] = resu[0]*resv[0];
  }
  if(flags & EVAL_1ST_DERIV) {
    ext[0] = resu[1]*resv[0];
    ext[1] = resu[0]*resv[1]; 
  }
  if(flags & EVAL_2ND_DERIV) {
    fnl[0] = resu[2]*resv[0];
    fnl[1] = resu[1]*resv[1];
    fnl[2] = resu[0]*resv[2];
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::evalall(int _lvl, int _gen){
  
  int TTL = pow2(_lvl);
  int NS = TTL-1;// * 1/2; //to render half the chart
  //for(int V=0; V<numVertices(); V++) {
  for(int V=0; V<1; V++) {
    int K = valence(V);
    if(K==0) continue;
    if (K==4 &&_gen == 2) continue;
    for(int f=0; f<K; f++) {
      double step = 1.0/double(TTL);      
      for(int j=0; j<=NS; j++)
	for(int i=0; i<=NS; i++) {
	  double cd[2];  
	  cd[0] = i*step; cd[1] = j*step;
	  double xy[2];
	  Vfcd2Vxy(EVAL_VALUE, V, f, cd, xy); 
	  Point3 ret[6];
	  eval( EVAL_VALUE|EVAL_1ST_DERIV, V, xy, ret);
	  
	}
    }
  }
  
  return 0;
}



// ---------------------------------------------------------------------- 
int BdSurf::blendFunc1D(int flags, double u, double lb, double ub, double* res){
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
  
  double* val = res;
  double* ext = res+1;
  double* fnl = res+2;
  double t = (u-lb)/(ub-lb);
  double ERR = 1e-7;

  if(t<=0+ERR) {
      //--------------------------------------------------
      if(flags & EVAL_VALUE) val[0]=1;
      if(flags & EVAL_1ST_DERIV) ext[0]=0;
      if(flags & EVAL_2ND_DERIV) fnl[0]=0;
  } else if(t>=1-ERR) {
      //--------------------------------------------------
      if(flags & EVAL_VALUE) val[0]=0;
      if(flags & EVAL_1ST_DERIV) ext[0]=0;
      if(flags & EVAL_2ND_DERIV) fnl[0]=0;
  } 
  
  else
      
  {
      //--------------------------------------------------
      if(_blendFuncType==ANALYTIC_BLEND_FUNC ) {
          analyticBlendFunc1D(flags, t, lb, ub, res);
      } else if(_blendFuncType==BSPLINE_BLEND_FUNC) {
          splineBlendFunc1D(flags, _splineBlendFuncDeg, t, 2*(ub-lb), res);
          
      } else if (_blendFuncType == OPTIM_BLEND_FUNC){
          optimBlendFunc1D(flags, _splineBlendFuncDeg, 1-t, (ub-lb),res);
      } 

      res[1] = res[1]/(ub-lb);
      res[2] = res[2]/((ub-lb)*(ub-lb));
                       

  }
  return 0;

}
// ---------------------------------------------------------------------- 
int BdSurf::optimBlendFunc1D(int flags, int deg, double z, double scale,double* res){
   

    assert(deg<=10);
    //assert(z>=0 && z<=1); 
    double z2,z3,z4,z5,z6,z7,z8,z9,z10;
    z2 = z*z;  z3 = z2*z;  z4 = z3*z;  z5 = z4*z;  z6=z5*z;   z7= z6*z;    z8 = z7*z;  z9 = z8*z;  z10 = z9*z;

  
    switch(deg){
    case 2:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= .50000000000000000000) res[0] =2.*z2;
        else if(z <= 1.) res[0] = -1.+(4.-2.*z)*z;
        else res[0] = 1;
        break;
    case 3:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= .2500000000) res[0] = 5.333333333*z3;
        else if (z <= .7500000000) res[0] = .1666666667-2.*z+8.*z2-5.333333333*z3;
        else if (z <= 1.) res[0] = -4.333333333+16.*z-16.*z2+5.333333333*z3;
        else res[0] = 1;
        break;
    case 4:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= .14644660940672623780) res[0] = 16.*z4;
        else if (z <= .50000000000000000000) res[0] = -0.14718625761429707190e-1+.40202025355333863355*z-4.1177490060914376575*z2+18.745166004060958438*z3-16.*z4;
        else if (z <= .85355339059327376220) res[0] = 1.9852813742385702924-15.597979746446661366*z+43.882250993908562342*z2-45.254833995939041560*z3+16.*z4;
        else if ( z <= 1.) res[0] = -14.999999999999999999+64.*z-96.*z2+64.*z3-16.*z4;
        else res[0] = 1.0000000000000000010+0.53573911150189318880e-19*z-0.26786955575094659440e-19*z2;
        break;
    case 5:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= 0.95491502812526287949e-1) res[0] = 51.200000000000000000*z5;
        else if (z <= .34549150281252628795) res[0] = 0.81306187557833487499e-3-0.42572472504416375409e-1*z+.89164944001345942980*z2-9.3374741600201891447*z3+48.891649440013459430*z4-51.200000000000000000*z5;
        else if (z <= .65450849718747371205) res[0] = -.50325224750231333941+7.2523292150113563936*z-41.337474160020189144*z2+112.89164944001345943*z3-128.*z4+51.200000000000000000*z5;
        else if (z <= .90450849718747371205) res[0] = 11.795934690622108325-86.705098312484227236*z+245.77087639996635144*z2-325.77087639996635143*z3+207.10835055998654057*z4-51.200000000000000000*z5;
        else if ( z <= 1.) res[0] = -50.199999999999999992+255.99999999999999999*z-512.*z2+512.*z3-256.*z4+51.200000000000000000*z5;
        else res[0] = 1.0000000000000000005+0.50335829998180960801e-17*z-0.90945424859373685600e-17*z2+0.20480000000000000000e-17*z3-0.51200000000000000000e-18*z4;
        break;
    case 6:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= 0.66987298107780676618e-1) res[0] = 170.66666666666666667*z6;
        else if ( z <= .25000000000000000000) res[0] = -0.30841356309254049329e-4+0.27624362092913055274e-2*z-.10309552285743124926*z2+2.0520412231296636895*z3-22.974966311837064286*z4+137.18998652473482572*z5-170.66666666666666667*z6;
        else if (z <= .50000000000000000000) res[0] = 0.83302491977024079238e-1-1.9972375637907086944*z+19.896904477142568750*z2-104.61462544353700297*z3+297.02503368816293572*z4-374.81001347526517428*z5+170.66666666666666667*z6;
        else if ( z <= .75000000000000000000) res[0] = -5.2500308413563092643+62.002762436209291334*z-300.10309552285743126*z2+748.71870788979633037*z3-982.97496631183706432*z4+649.18998652473482570*z5-170.66666666666666667*z6;
        else if (z <= .93301270189221932338) res[0] = 55.499969158643690793-423.99723756379070876*z+1319.8969044771425687*z2-2131.2812921102036695*z3+1897.0250336881629356*z4-886.81001347526517424*z5+170.66666666666666667*z6;
        else if (z <= 1.) res[0] = -169.66666666666666681+1023.9999999999999999*z-2559.9999999999999998*z2+3413.3333333333333333*z3-2560.*z4+1024.*z5-170.66666666666666667*z6;
        else res[0]= 1.0000000000000001158-0.16720452122920761250e-15*z+0.18244964361477047929e-15*z2-0.17930021879695548175e-15*z3+0.62756570302591872875e-16*z4-0.40960000000000000000e-17*z5;
        break;
    case 7:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= 0.49515566048790436882e-1) res[0] = 585.14285714285714286*z7;
        else if (z <= .18825509907063323474) res[0] =0.85405165909716607756e-6-0.12073701445297728894e-3*z+0.73150944695255066160e-2*z2-.24622204871620700105*z3+4.9726190845438532088*z4-60.255222525103182580*z5+405.63151707169125893*z6-585.14285714285714286*z7;
        else if (z <= .38873953302184279786) res[0] = -0.98057613154314852871e-2+.36452440624340911689*z-5.8036054870927415800*z2+51.199224495620706333*z3-268.30256335837421735*z4+810.71778312098134436*z5-1136.5542545149362001*z6+585.14285714285714286*z7; 
        else if ( z <= .61126046697815720214) res[0] = 1.5602058488570656905-27.906542304612471398*z+212.37127236306429860*z2-884.19521195470236570*z3+2137.9215455689892500*z4-2903.1686182275957000*z5+2048.*z6-585.14285714285714286*z7;
        else if ( z <= .81174490092936676526) res[0] = -35.754160043984012320+399.40812200071269746*z-1884.8427013352945168*z2+4834.0817119732132800*z3-7216.9725345224895032*z4+6279.3922560313641436*z5-2959.4457454850638000*z6+585.14285714285714286*z7;
        else if ( z <= .95048443395120956312) res[0] = 236.03297034893468923-1944.3206905514169708*z+6776.9751057200789195*z2-12950.277629527747442*z3+14691.830737465603176*z4-9914.4661200949556290*z5+3690.3684829283087410*z6-585.14285714285714286*z7;
        else if(z <= 1.) res[0] = -584.14285714285714216+4095.9999999999999974*z-12287.999999999999996*z2+20479.999999999999998*z3-20480.*z4+12288.*z5-4096.*z6+585.14285714285714286*z7;
        else res[0]= 1.0000000000000000484-0.36175677281249441305e-15*z+0.63749634654633829355e-15*z2-0.26517790529631145793e-16*z3-0.42372151766370038882e-15*z4+0.36623100231065481024e-16*z5+0.16384000000000000000e-16*z6;
        break;
    case 8:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= 0.38060233744356621936e-1) res[0] = 2048.*z8;
        else if (z <= .14644660940672623780) res[0] = -0.18035639965439994271e-7+0.37909677773566962999e-5*z-0.34861549484613467622e-3*z2+0.18319146287314938525e-1*z3-.60164982204133232780*z4+12.646266464520426781*z5-166.13490276311756640*z6+1247.1577393350777876*z7-2048.*z8;
        else if (z <= .30865828381745511414) res[0] = 0.86653374158012347544e-3-0.47333694316756025668e-1*z+1.1309933743747886546*z2-15.432253515692972281*z3+131.27723700036884408*z4-707.77402790306240450*z5+2293.5338368755011943*z6-3551.6047577045275726*z7+2048.*z8;
        else if ( z <= .50000000000000000000) res[0] = -.33656327434686270681+8.6983848113431433435*z-98.040214315981076030*z2+627.16321785035092847*z3-2471.0969137047681297*z4+6037.2233679723528780*z5-8632.7846024904456260*z6+6562.5098864258416076*z7-2048.*z8;
        else if (z <= .69134171618254488586) res[0] = 15.663436725653137396-247.30161518865685704*z+1693.9597856840189244*z2-6540.8367821496490713*z3+15448.903086295231870*z4-22634.776632027647122*z5+20039.215397509554373*z6-9821.4901135741583924*z7+2048.*z8;
        else if(z <= .85355339059327376220) res[0] = -198.08456096638670317+2226.1264183148314156*z-10828.065017136323776*z2+29684.446633306717202*z3-50049.248130959109695*z4+53157.729081554865736*z5-34776.300532943808187*z6+12832.395242295472428*z7-2048.*z8;
        else if(z <= .96193976625564337806) res[0]= 955.91457248183607910-8589.8262441998840622*z+33522.803640873806616*z2-74236.102794031302533*z3+102138.87298221848013*z4-89481.850624077551436*z5+48780.030727417573053*z6-15136.842260664922211*z7+2048.*z8;
        else if(z <= 1.) res[0] = -2047.0000000000000043+16384.000000000000064*z-57344.000000000000130*z2+114688.00000000000009*z3-143360.00000000000002*z4+114688.*z5-57344.*z6+16384.*z7-2048.*z8;
        else res[0] = 1.0000000000000049597-0.54820642047180343368e-14*z+0.39224496026381535292e-14*z2-0.11250317116645488317e-13*z3+0.10613072687253747715e-13*z4-0.51840487459490730284e-14*z5+0.13760453616361108574e-14*z6-0.13107200000000000000e-15*z7;
        break;
    case 9:
        if      (z<= 0) res[0] = 0.0;
        else if (z <= 0.30153689607045807973e-1) res[0] = 7281.7777777777777778*z9;
        else if (z <= .11697777844051098240) res[0] = 0.30014530688699986480e-9-0.89584651072079832882e-7*z+0.11883739899099736676e-4*z2-0.91957988533361885860e-3*z3+0.45744645049278490318e-1*z4-1.5170496760233828764*z5+33.540388054068717260*z6-476.70624094417708557*z7+3952.3044041747081425*z8-7281.7777777777777778*z9;
        else if (z <= .25000000000000000000) res[0] = -0.59730418570177560056e-4+0.45954539075066661977e-2*z-.15713056001046843486*z2+3.1335706367386609790*z3-40.147661210465904882*z4+342.08159265890118616*z5-1924.6585503289824747*z6+6697.5462336543921967*z7-11380.206971579947342*z8+7281.7777777777777778*z9;
        else if ( z <= .41317591116653482557) res[0] = 0.55495825136985377258e-1-1.9954045460924933281*z+31.842869439989531566*z2-295.53309602992800583*z3+1751.8523387895340956*z4-6825.9184073410988144*z5+17190.008116337684193*z6-26070.453766345607804*z7+21387.793028420052658*z8-7281.7777777777777778*z9;
        else if (z <= .58682408883346517443) res[0]= -5.0552273567844730368+109.32887016157253927*z-1045.8992902455457948*z2+5790.8132711456916207*z3-20344.110431614285835*z4+46652.423802928302734*z5-69098.246296949255063*z6+63433.022751509310967*z7-32768.*z8+7281.7777777777777778*z9;
        else if (z <= .75000000000000000000) res[0] = 115.12660322810743222-1733.8750743359833380*z+11518.029011754604590*z2-44165.951051435753480*z3+107351.98778554328413*z4-170952.98921104980870*z5+178114.29678565342937*z6-117112.10953898518655*z7+44148.206971579947345*z8-7281.7777777777777778*z9;
        else if(z <= .88302222155948901760) res[0] = -978.37339677189257666+11388.124925664016721*z-58465.970988245395545*z2+173562.04894856424670*z3-328104.01221445671602*z4+409655.01078895019136*z5-337981.70321434657063*z6+177799.89046101481346*z7-54155.793028420052661*z8+7281.7777777777777778*z9;
        else if(z <= .96984631039295419203) res[0] = 3775.1114393095821432-37060.671129748843081*z+161002.10069787071451*z2-406369.18580927442303*z3+657031.84382374066480*z4-705986.05914739567384*z5+504308.21331499667623*z6-231002.27100754651194*z7+61583.695595825291858*z8-7281.7777777777777778*z9;
        else if(z <= 1.) res[0] = -7280.7777777777780210+65536.000000000000220*z-262143.99999999999959*z2+611669.33333333333273*z3-917503.99999999999968*z4+917503.99999999999990*z5-611669.33333333333333*z6+262144.*z7-65536.*z8+7281.7777777777777778*z9;
        else res[0] = .99999999999996338747-0.38084892777430476847e-13*z+0.38100254603956599662e-12*z2-0.53798512751231230287e-12*z3+0.29695038035979669252e-12*z4-0.45208610065940934494e-13*z5-0.86211581928152639460e-14*z6-0.34984443524634146274e-16*z7+0.39321600000000000000e-15*z8;
        break;
    case 10:
        if      (z<= 0) res[0] = 0.0;
        else if ( z <= 0.24471741852423213942e-1) res[0] = 26214.400000000000000*z10;
        else if (z <= 0.95491502812526287949e-1) res[0] = -0.40384875164793945666e-11+0.16502656577670215063e-8*z-0.30346002768152967394e-6*z2+0.33067803075241608518e-4*z3-0.23647133796461901013e-2*z4+.11595643958194340222*z5-3.9486509338409739117*z6+92.203161336119367790*z7-1412.9025105591736969*z8+12830.240592323261991*z9-26214.400000000000000*z10;
        else if(z <= .20610737385376343542) res[0] = 0.33053440291072819838e-5-0.34613889315882092304e-3*z+0.16311435276225912962e-1*z2-.45548358741953987547*z3+8.3455418873216695582*z4-104.78853899503713546*z5+911.52979780112876273*z6-5386.0909981065753686*z7+20100.637954127342092*z8-37234.808434250520462*z9+26214.400000000000000*z10;
        else if(z <= .34549150281252628795) res[0] =-0.72494865223342578718e-2+.35154771772884343771*z-7.6666856608031469320*z2+98.948976685545888113*z3-835.66989640124444868*z4+4809.2448942041588794*z5-18956.890164236811193*z6+49698.705936454644500*z7-80122.844904456821299*z8+70824.814388791403557*z9-26214.400000000000000*z10;
        else if(z <= .50000000000000000000) res[0] = 1.2631596940688140388-36.419521042260577435*z+471.27363001388107786*z2-3597.7386834608254167*z3+17888.968343881403582*z4-60227.277582968936504*z5+137912.67572410236643*z6-209757.02946215209234*z7+201493.05587002070501*z8-110312.23463778237889*z9+26214.400000000000000*z10;
        else if(z <= .65450849718747371205) res[0]= -49.936840305931193005+987.58047895773946460*z-8744.7263699861190370*z2+45554.261316539174767*z3-154143.03165611859655*z4+352649.52241703106350*z5-550215.32427589763360*z6+576674.97053784790771*z7-388330.94412997929500*z8+151831.76536221762111*z9-26214.400000000000000*z10;
        else if(z <= .79389262614623656458) res[0] = 706.41315638872070864-10568.419174901716989*z+70707.257318275144570*z2-278157.28316680560265*z3+711384.62043728070230*z4-1234240.7730875343176*z5+1470245.1972753673968*z6-1187318.7353027093976*z7+622347.51540533418924*z8-191319.18561120859644*z9+26214.400000000000000*z10; 
        else if(z <= .90450849718747371205) res[0] = -4507.7858074779671021+55110.474823368299615*z-301578.63281137226899*z2+972342.26927039137513*z3-2045126.9589535473005*z4+2932335.1420031392324*z5-2903326.8470495769608*z6+1960693.9090018934247*z7-864635.36204587265791*z8+224909.19156574947953*z9-26214.400000000000000*z10;
        else if(z <= .97552825814757678606) res[0]= 14709.693783341442707-157352.75410787171601*z+755442.42687295305590*z2-2143961.0430694622613*z3+3984148.4028061700548*z4-5066628.3355196876724*z5+4466203.5870620838613*z6-2695050.3555994998386*z7+1065588.7371796498158*z8-249313.75940767673801*z9+26214.400000000000000*z10;
        else if(z <= 1.) res[0]=-26213.399999999999546+262143.99999999999818*z-1179647.9999999999949*z2+3145727.9999999999930*z3-5505023.9999999999962*z4+6606028.7999999999996*z5-5505024.*z6+3145728.*z7-1179648.*z8+262144.*z9-26214.400000000000000*z10;
        else res[0]= 1.0000000000004344917-0.23771387168790120845e-11*z+0.60627974858379033630e-11*z2-0.70198079820713437333e-11*z3+0.37343903197696794288e-11*z4-0.77687433196916569664e-12*z5+0.45868652506776221160e-13*z6-0.21801870076643113956e-13*z7+0.54625363142958494196e-16*z8+0.15728640000000000000e-14*z9;
        break;

    }

    if(flags & EVAL_1ST_DERIV) {
      
        switch(deg){

        case 2:
            if(z<=0) res[1] = 0;
            else if(z <= .50000000000000000000) res[1] = 4.*z;
            else if(z <= 1.) res[1] = -4.*z+4.;
            else res[1] = 0;
            break;
        case 3:
            if(z<=0) res[1] = 0;
            else if(z < .25000000000000000000) res[1] =16.*z2;
            else if(z < .75000000000000000000) res[1] =-16.*z2+16.*z-2.; 
            else if(z < 1.) res[1]=16.*z2-32.*z+16.;
            else res[1] = 0;
            break;
        case 4:
            if(z<=0) res[1] = 0;
            else if (z < .14644660940672623781) res[1] = 64.*z3;
            else if (z < .50000000000000000000) res[1] = -64.*z3+56.235498012182875317*z2-8.2354980121828753150*z+.40202025355333863352;
            else if(z < .85355339059327376225) res[1] = 64.*z3-135.76450198781712468*z2+87.764501987817124684*z-15.597979746446661366;
            else if(z < 1.) res[1] = -64.*z3+192.*z2-192.*z+64;
            else res[1] = -0.53573911150189318881e-19*z+0.53573911150189318881e-19;
            break;
        case 5:
            if(z<=0) res[1] = 0;
            else if(z < 0.95491502812526287948e-1) res[1] = 256.*z4;
            else if(z < .34549150281252628792) res[1] = -256.*z4+195.56659776005383771*z3-28.012422480060567434*z2+1.7832988800269188596*z-0.42572472504416375407e-1;
            else if(z < .65450849718747371208) res[1] = 256.*z4-512.*z3+338.67494832004037829*z2-82.674948320040378282*z+7.2523292150113563936;
            else if(z < .90450849718747371204) res[1] = -256.*z4+828.43340223994616228*z3-977.31262919989905429*z2+491.54175279993270286*z-86.705098312484227238;
            else if(z < 1.) res[1] = 256.*z4-1024.*z3+1536.*z2-1024.*z+256.;
            else res[1] = -0.20480000000000000000e-17*z3+0.61440000000000000000e-17*z2-0.18189084971874737120e-16*z+0.50335829998180960800e-17;
            break;
        case 6:
            if(z<=0) res[1] = 0;
            else if(z < 0.66987298107780676620e-1) res[1] = 1024.*z5;
            else if(z < .25000000000000000000) res[1] = -1024.*z5+685.94993262367412860*z4-91.899865247348257152*z3+6.1561236693889910688*z2-.20619104571486249852*z+0.27624362092913055274e-2;
            else if(z < .50000000000000000000) res[1] = 1024.*z5-1874.0500673763258715*z4+1188.1001347526517428*z3-313.84387633061100891*z2+39.793808954285137496*z-1.9972375637907086945;
            else if(z < .75000000000000000000) res[1] = -1024.*z5+3245.9499326236741285*z4-3931.8998652473482572*z3+2246.1561236693889911*z2-600.20619104571486252*z+62.002762436209291333;
            else if(z < .93301270189221932341) res[1] = 1024.*z5-4434.0500673763258708*z4+7588.1001347526517424*z3-6393.8438763306110085*z2+2639.7938089542851374*z-423.99723756379070875;
            else if(z < 1.) res[1]=-1024.*z5+5120.*z4-10240.*z3+10240.*z2-5120.*z+1024.;
            else res[1] = -0.20480000000000000000e-16*z4+0.25102628121036749150e-15*z3-0.53790065639086644528e-15*z2+0.36489928722954095856e-15*z-0.16720452122920761250e-15;
            break;
        case 7:
            if(z<=0) res[1] = 0;
            else if(z < 0.49515566048790436879e-1) res[1] = 4096.*z6;
            else if(z < .18825509907063323475) res[1] = -4096.*z6+2433.7891024301475535*z5-301.27611262551591290*z4+19.890476338175412836*z3-.73866614614862100315*z2+0.14630188939051013232e-1*z-0.12073701445297728894e-3;
            else if(z < .38873953302184279784) res[1] = 4096.*z6-6819.3255270896172000*z5+4053.5889156049067218*z4-1073.2102534334968694*z3+153.59767348686211900*z2-11.607210974185483160*z+.36452440624340911689;
            else if(z < .61126046697815720216) res[1] = -4096.*z6+12288.*z5-14515.843091137978500*z4+8551.6861822759570000*z3-2652.5856358641070972*z2+424.74254472612859720*z-27.906542304612471398;
            else if(z < .81174490092936676525) res[1] = 4096.*z6-17756.674472910382800*z5+31396.961280156820718*z4-28867.890138089958013*z3+14502.245135919639841*z2-3769.6854026705890334*z+399.40812200071269745; 
            else if(z < .95048443395120956303) res[1] = -4096.*z6+22142.210897569852444*z5-49572.330600474778142*z4+58767.322949862412704*z3-38850.832888583242326*z2+13553.950211440157839*z-1944.3206905514169707;
            else if(z < 1.) res[1] = 4096.*z6-24576.*z5+61440.*z4-81920.*z3+61440.*z2-24575.999999999999992*z+4095.9999999999999974;
            else res[1] = 0.98304000000000000000e-16*z5+0.18311550115532740512e-15*z4-0.16948860706548015553e-14*z3-0.79553371588893437385e-16*z2+0.12749926930926765871e-14*z-0.36175677281249441306e-15;
            break;
        case 8:
            if(z<=0) res[1] = 0;
            else if(z < 0.38060233744356621938e-1) res[1] = 16384.*z7;
            else if(z < .14644660940672623781) res[1] = -16384.*z7+8730.1041753455445132*z6-996.80941657870539846*z5+63.231332322602133910*z4-2.4065992881653293110*z3+0.54957438861944815575e-1*z2-0.69723098969226935242e-3*z+0.37909677773566962996e-5; 
            else if(z < .30865828381745511415) res[1] = 16384.*z7-24861.233303931693010*z6+13761.203021253007166*z5-3538.8701395153120225*z4+525.10894800147537628*z3-46.296760547078916843*z2+2.2619867487495773094*z-0.47333694316756025672e-1;
            else if(z < .50000000000000000000) res[1] = -16384.*z7+45937.569204980891252*z6-51796.707614942673757*z5+30186.116839861764392*z4-9884.3876548190725192*z3+1881.4896535510527853*z2-196.08042863196215206*z+8.6983848113431433435;
            else if(z < .69134171618254488585) res[1] = 16384.*z7-68750.430795019108748*z6+120235.29238505732623*z5-113173.88316013823562*z4+61795.612345180927480*z3-19622.510346448947213*z2+3387.9195713680378490*z-247.30161518865685704;
            else if(z < .85355339059327376225) res[1] = -16384.*z7+89826.766696068306996*z6-208657.80319766284911*z5+265788.64540777432867*z4-200196.99252383643877*z3+89053.339899920151606*z2-21656.130034272647554*z+2226.1264183148314156;
            else if(z < .96193976625564337808) res[1] = 16384.*z7-105957.89582465445548*z6+292680.18436450543832*z5-447409.25312038775716*z4+408555.49192887392052*z3-222708.30838209390760*z2+67045.607281747613232*z-8589.8262441998840621;
            else if(z < 1.) res[1] = -16384.*z7+114688.*z6-344064.*z5+573440.*z4-573440.00000000000008*z3+344064.00000000000027*z2-114688.00000000000026*z+16384.000000000000064;
            else res[1] = -0.91750400000000000000e-15*z6+0.82562721698166651444e-14*z5-0.25920243729745365142e-13*z4+0.42452290749014990856e-13*z3-0.33750951349936464951e-13*z2+0.78448992052763070588e-14*z-0.54820642047180343371e-14;
            break;
        case 9:
            if(z<=0) res[1] = 0;
            else if(z < 0.30153689607045807973e-1) res[1] = 65536.*z8;
            else if(z < .11697777844051098241) res[1] = -65536.000000000000001*z8+31618.435233397665141*z7-3336.9436866092395990*z6+201.24232832441230355*z5-7.5852483801169143815*z4+.18297858019711396128*z3-0.27587396560008565757e-2*z2+0.23767479798199473350e-4*z-0.89584651072079832878e-7;
            else if(z < .25000000000000000000)res[1] = 65536.000000000000001*z8-91041.655772639578736*z7+46882.823635580745377*z6-11547.951301973894849*z5+1710.4079632945059308*z4-160.59064484186361953*z3+9.4007119102159829364*z2-.31426112002093686972*z+0.45954539075066661973e-2;
            else if(z < .41317591116653482558) res[1] = -65536.000000000000001*z8+171102.34422736042125*z7-182493.17636441925463*z6+103140.04869802610516*z5-34129.592036705494071*z4+7007.4093551581363820*z3-886.59928808978401749*z2+63.685738879979063132*z-1.9954045460924933281;
            else if(z < .58682408883346517442) res[1] = 65536.000000000000001*z8-262144.*z7+444031.15926056517677*z6-414589.47778169553038*z5+233262.11901464151367*z4-81376.441726457143344*z3+17372.439813437074863*z2-2091.7985804910915896*z+109.32887016157253927;
            else if(z < .75000000000000000000) res[1] = -65536.000000000000001*z8+353185.65577263957875*z7-819784.76677289630585*z6+1068685.7807139205762*z5-854764.94605524904355*z4+429407.95114217313656*z3-132497.85315430726044*z2+23036.058023509209180*z-1733.8750743359833380;
            else if(z < .88302222155948901759) res[1] = 65536.000000000000001*z8-433246.34422736042125*z7+1244599.2332271036942*z6-2027890.2192860794238*z5+2048275.0539447509570*z4-1312416.0488578268641*z3+520686.14684569274013*z2-116931.94197649079109*z+11388.124925664016721;
            else if(z < .96984631039295419207) res[1] = -65536.000000000000001*z8+492669.56476660233489*z7-1617015.8970528255834*z6+3025849.2798899800574*z5-3529930.2957369783692*z4+2628127.3752949626591*z3-1219107.5574278232691*z2+322004.20139574142904*z-37060.671129748843081;
            else if(z < 1.) res[1] = 65536.000000000000001*z8-524288.*z7+1835008.*z6-3670016.*z5+4587519.9999999999995*z4-3670015.9999999999987*z3+1835007.9999999999982*z2-524287.99999999999918*z+65536.000000000000220;
            else res[1] = 0.31457280000000000000e-14*z7-0.24489110467243902392e-15*z6-0.51726949156891583675e-13*z5-0.22604305032970467247e-12*z4+0.11878015214391867700e-11*z3-0.16139553825369369086e-11*z2+0.76200509207913199320e-12*z-0.38084892777430476847e-13;
            break;
        case 10:
            if(z<=0) res[1] = 0;
            else if(z < 0.24471741852423213942e-1) res[1] = 262144.*z9;
            else if(z < 0.95491502812526287948e-1) res[1] = -262144.*z9+115472.16533090935792*z8-11303.220084473389574*z7+645.42212935283557452*z6-23.691905603045843470*z5+.57978219790971701110*z4-0.94588535185847604044e-2*z3+0.99203409225724825554e-4*z2-0.60692005536305934786e-6*z+0.16502656577670215063e-8;
            else if(z < .20610737385376343542) res[1] = 262144.*z9-335113.27590825468416*z8+160805.10363301873673*z7-37702.636986746027579*z6+5469.1787868067725761*z5-523.94269497518567730*z4+33.382167549286678230*z3-1.3664507622586196264*z2+0.32622870552451825926e-1*z-0.34613889315882092304e-3;
            else if(z < .34549150281252628792) res[1] = -262144.*z9+637423.32949912263201*z8-640982.75923565457040*z7+347890.94155518251150*z6-113741.34098542086716*z5+24046.224471020794398*z4-3342.6795856049777948*z3+296.84693005663766433*z2-15.333371321606293864*z+.35154771772884343770;
            else if(z < .50000000000000000000) res[1] =262144.*z9-992810.11174004141001*z8+1611944.4469601656401*z7-1468299.2062350646463*z6+827476.05434461419858*z5-301136.38791484468250*z4+71555.873375525614328*z3-10793.216050382476250*z2+942.54726002776215574*z-36.419521042260577435;
            else if(z < .65450849718747371208) res[1] = -262144.*z9+1366485.8882599585900*z8-3106647.5530398343600*z7+4036724.7937649353540*z6-3301291.9456553858019*z5+1763247.6120851553176*z4-616572.12662447438620*z3+136662.78394961752431*z2-17489.452739972238074*z+987.58047895773946459;
            else if(z < .79389262614623656452) res[1] = 262144.*z9-1721872.6705008773680*z8+4978780.1232426735137*z7-8311231.1471189657839*z6+8821471.1836522043814*z5-6171203.8654376715880*z4+2845538.4817491228090*z3-834471.84950041680795*z2+141414.51463655028914*z-10568.419174901716989;
            else if(z < .90450849718747371204) res[1] = -262144.*z9+2024182.7240917453158*z8-6917082.8963669812632*z7+13724857.363013253973*z6-17419961.082297461765*z5+14661675.710015696162*z4-8180507.8358141892028*z3+2917026.8078111741257*z2-603157.26562274453794*z+55110.474823368299617;
            else if(z < .97552825814757678607) res[1] = 262144.*z9-2243823.8346690906421*z8+8524709.8974371985264*z7-18865352.489196498872*z6+26797221.522372503165*z5-25333141.677598438360*z4+15936593.611224680219*z3-6431883.1292083867842*z2+1510884.8537459061118*z-157352.75410787171601;
            else if(z < 1.) res[1] = -262144.*z9+2359296.*z8-9437184.*z7+22020096.*z6-33030144.*z5+33030144.*z4-22020095.999999999985*z3+9437183.9999999999790*z2-2359295.9999999999898*z+262143.99999999999818;
            else res[1] = 0.14155776000000000000e-13*z8+0.43700290514366795354e-15*z7-0.15261309053650179771e-12*z6+0.27521191504065732696e-12*z5-0.38843716598458284833e-11*z4+0.14937561279078717714e-10*z3-0.21059423946214031200e-10*z2+0.12125594971675806726e-10*z-0.23771387168790120845e-11;
            break;
            
        }
        
         
//        res[1] = res[1]/scale;
        
        
        
    }
    
    if(flags & EVAL_2ND_DERIV) {
        assert(0);     
    }
    
    return 0; 
}
// ---------------------------------------------------------------------- 
int BdSurf::analyticBlendFunc1D(int flags, double t, double lb, double ub, double* res){
    
    double* val = res;
  double* ext = res+1;
  double* fnl = res+2;
  
  double s = 1-t;
  double t2 = t*t;    double t3 = t2*t;  double t4 = t2*t2;
  double s2 = s*s;  double s3 = s2*s;  double s4 = s2*s2;
  double a = 2*exp(-1/t)/(t-1);
  double b = 2*exp(-1/s)/(s-1);
  double ea = exp(a);
  double eb = exp(b);
  double f = ea;
  double g = ea+eb;

  double da =  a*(1/(t*t) - 1/(t-1));
  double db = -b*(1/(s*s) - 1/(s-1));
  double dda = a*(-4*t3+7*t2-4*t+1+2*t4)/((t-1)*(t-1))/t4;
  double ddb = b*(-4*s3+7*s2-4*s+1+2*s4)/((s-1)*(s-1))/s4;
  
  double df = ea*da;
  double dg = ea*da + eb*db;
  double ddf = ea*da*da + ea*dda;
  double ddg = ea*da*da + ea*dda + eb*db*db + eb*ddb;
  
  if(flags & EVAL_VALUE) {
    val[0] = f/g;
  }
  if(flags & EVAL_1ST_DERIV) {
      ext[0] = (df/g - f/(g*g)*dg) ;
  }
  if(flags & EVAL_2ND_DERIV) {
      fnl[0] = (ddf/g - 2*df/(g*g)*dg + 2*f/(g*g*g)*dg*dg - f/(g*g)*ddg);
  }
  
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::splineBlendFunc1D(int flags,int deg, double t, double scale, double * res){ 
    //spline blending function of arbitrary degree
  if (deg==3){
    double x= 3.0*t;
    int tx= (int)ceil(x); 
    double t = x-tx+1;
    double t2 = t*t; double t3 = t2*t;
    switch(tx){
    case 1:{ 
      res[0] = (1.0/6.0)*((-t3+3.0*t2-3.0*t+1.0)+(4.0-6.0*t2+3.0*t3)+
			  (1.0+3.0*t+3.0*t2-3.0*t3));
      res[1] = -3.0*t2/2.;
    } break;
    case 2:{
      res[0] = (1.0/6.0)*((-t3+3.0*t2-3.0*t+1.0)+(4.0-6.0*t2+3.0*t3));
      res[1] = (1.0-2.0*t2+2.0*t)*(-3.0)/2.0;
    } break;
    case 3:{
      res[0] = (1.0/6.0)*(-t3+3.0*t2-3.0*t+1.0);
      res[1] = (1.0+t2-2.0*t)*(-3.0)/2.0;
    }break;
    default: assert(0); break;
    }
    
    return(0);
  }
  else{
    
    double *poucp = new double[2*deg];
    for (int i=0; i<deg; i++){ 
      poucp[i] = 1.0;
      poucp[i+deg] = 0.0;
    }
    
    double x = deg*t; 
    int tx = (int)ceil(x);
    int e = (int)floor((deg-1)/2.0);
    double ex = ceil((deg-1)/2.0);
    double *Nj = new double[2*deg];
    
    if(flags & EVAL_VALUE){
      
      res[0] = 0.0;
      SplineBasis(x, tx,deg,2*deg, Nj);
      for(int i=tx-(e+1); i<=tx+e; i++)
	res[0]+=Nj[i+e]*poucp[i+e];      
      res+=1;
    }
    
    if(flags &EVAL_1ST_DERIV){
      res[0] = 0.0;
      SplineBasis(x, tx, deg-1,2*deg,Nj);
      for(int i=tx-(e+1); i<tx+ex; i++){      
	res[0] += Nj[i+e]*(poucp[i+1+e] - poucp[i+e]);
      }
      res[0] *= (double(2*deg)/2.);
      res+=1;
    }
    if(flags & EVAL_2ND_DERIV){
      assert(0);
    }
    
    delete[] Nj;
    delete[] poucp;
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::chartEval(int flags, int V, double* xy, Point3* ret,int flag){ 

  assert( (flags==EVAL_VALUE) || 
	  (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || 
	  (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
  
  if(_basisType == SPLINE_BASIS &&
     _splineDeg == 3 && gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX 
     /*** DZ ***/
&& !LoadW
     /*** DZ ***/
){
    evalCubicBSpline(flags, V, xy, ret);  //faster version only for C2, and closed surfaces
    return(0);
  }
  else if(_basisType == POLYNOMIAL_BASIS) {
      evalPoly(flags, V, xy, ret, flag);   
  }
  else if (_basisType == SPLINE_BASIS){
      evalBSpline(flags, V,xy, ret);
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::evalPoly(int flags, int V, double* xy, Point3* ret, int flag){     
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV)) );
  
  Vector<Point3>* ivp = _controls[V];
  int K = valence(V);  
  assert( ivp!=NULL );
  Vector<Point3>& vp = *ivp;
  double x = xy[0];	 double y = xy[1];
  /*
  if(V==18){
      for(int i=0; i<vp.length(); i++)
          cout<<vp(i)(2)<<endl;
      exit(1);
  }
  */


  int DEG;
  
  if (K>=12 && _max_degree_valence != 0) DEG =_max_degree_valence+2;
  else if (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX && K>=6 ) DEG = K; 
  else if (gpmesh().boun(V) == GpMesh::BOUNDARY_VERTEX && K == 2) DEG = K+4;
  else if (K==1) DEG = K+3;
  else    DEG = K+2;
  
  int cn = (DEG+1)*DEG/2;
  assert(cn==vp.length());

  double xs[50], ys[50];
  xs[0] = 1;		ys[0] = 1;
  for(int d=1; d<DEG; d++) {
      xs[d] = xs[d-1]*x;
      ys[d] = ys[d-1]*y;
  }
  Vector<Point3> t; 
  t.resize(cn);
  if(flag == 1){
//      cout<<"here"<<endl;
      t(0) = Point3(1.0); 
      for(int i=1; i<cn; i++)
          t(i) = Point3(0.0);
  }

  double tmp[1000];
  Point3* cur = ret;
  if(flags & EVAL_VALUE) {
    {
      int j = 0;
      for(int d=0; d<DEG; d++) {
	for(int t=0; t<=d; t++)
	  tmp[j++] = xs[d-t] * ys[t];
      }
      assert(j==cn);
      *cur = Point3( 0.0); //clear
      for(int i=0; i<cn; i++){
          if(flag)
              (*cur) += tmp[i] * t(i);
          else
              (*cur) += tmp[i] * vp(i);
      }
      cur = cur+1;
    }
  }
  if(flags & EVAL_1ST_DERIV) {
    {
      int j = 0;
      for(int d=0; d<DEG; d++) {
	for(int t=0; t<d; t++)
	  tmp[j++] = (d-t) * xs[d-t-1] * ys[t];
	tmp[j++] = 0;
      }
         assert(j==cn);
      *cur = Point3( 0.0); //clear
      for(int i=0; i<cn; i++){
          if(flag)
              (*cur) += tmp[i] * t(i);
          else
              (*cur) += tmp[i] * vp(i);
      }
      cur = cur+1;
    }
    {
      int j = 0;
      for(int d=0; d<DEG; d++) {
	tmp[j++] = 0;
	for(int t=1; t<=d; t++) {
	  tmp[j++] = t * xs[d-t] * ys[t-1];
	}
      }
       assert(j==cn);
      *cur = Point3( 0.0); //clear
      for(int i=0; i<cn; i++){
          if(flag)
              (*cur) += tmp[i] * t(i);
          else
              (*cur) += tmp[i] * vp(i);
      }
      cur = cur+1;
    }
  }
  if(flags & EVAL_2ND_DERIV) {
      assert(0);
    {
      int j=0;
      for(int d=0; d<DEG; d++) {
	for(int t=0; t<=d; t++)
	  if(d-t-2>=0) 
	    tmp[j++] = (d-t)*(d-t-1)*xs[d-t-2] * ys[t];
	  else tmp[j++] = 0;
      }
      assert(j==cn);
      *cur = Point3(0.0);
      for(int i=0; i<cn; i++)
	(*cur) += tmp[i] * vp(i);
      cur = cur+1;
    }
    {
      int j = 0;
      for(int d=0; d<DEG; d++) {
	for(int t=0; t<=d; t++)
	  if(d-t-1>=0 && t-1>=0) 
	    tmp[j++] = (d-t)*t*xs[d-t-1]*ys[t-1];
	  else tmp[j++] = 0;
      }
      assert(j==cn);
      *cur = Point3(0.0);
      for(int i=0; i<cn; i++)
	(*cur) += tmp[i] * vp(i);
      cur = cur+1;
    }
    {
      int j=0;
      for(int d=0; d<DEG; d++) {
	for(int t=0; t<=d; t++)
	  if(t-2>=0) 
	    tmp[j++] = t*(t-1)*xs[d-t]*ys[t-2];
	  else tmp[j++] = 0;
      }
      assert(j==cn);
      *cur = Point3(0.0);
      for(int i=0; i<cn; i++)
	(*cur) += tmp[i] * vp(i);
      cur = cur+1;
    }
  }
  return 0;
}


// ---------------------------------------------------------------------- 
inline int BdSurf::evalCubicBSpline(int flags, int V, double* xy, Point3* ret){
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV))) ;
  Vector<Point3>* ivp = _controls[V];
  Vector<Point3>& vp = *ivp;
  
  double min = 1.8; 
  const double x = xy[0]; const double y = xy[1];
  int nrgrid = 3;
  
  const int np2 = nrgrid+2;  
  
  const double newx = (x+min)*(nrgrid-1)*(1.0/(2*min));
  const double newy = (y+min)*(nrgrid-1)*(1.0/(2*min));
  
  int tx = (int)newx; int ty = (int)newy;   
  if (tx != newx) tx++;
  if (ty != newy) ty++;
  
  double Nj[4];  double Njd[4];
  double t = newx-tx+1;
  double t2 = t*t; double t3 = t2*t;
  Nj[0] = (-t3+3.0*t2-3.0*t+1.0)*(1.0/6.0);
  Nj[1] = (4.0-6.0*t2+3.0*t3)*(1.0/6.0);
  Nj[2] = (1.0+3.0*t+3.0*t2-3.0*t3)*(1.0/6.0);
  Nj[3] = t3*(1.0/6.0);
  
  Njd[0] = (1.0+t2-2.0*t)*(1.0/2.0);
  Njd[1] = (1.0-2.0*t2+2.0*t)*(1.0/2.0);
  Njd[2] = t2*(1.0/2.0);
  
  double Nk[4]; double Nkd[4];
  t = newy-ty+1; 
  t2 = t*t; t3= t2*t;
  Nk[0] = (-t3+3.0*t2-3.0*t+1.0)*(1.0/6.0);
  Nk[1] = (4.0-6.0*t2+3.0*t3)*(1.0/6.0);
  Nk[2] = (1.0+3.0*t+3.0*t2-3.0*t3)*(1.0/6.0);
  Nk[3] = t3*(1.0/6.0);
  
  Nkd[0] = (1.0+t2-2.0*t)*(1.0/2.0);
  Nkd[1] = (1.0-2.0*t2+2.0*t)*(1.0/2.0);
  Nkd[2] = t2*(1.0/2.0);
  
  
  Vec3T<double>* vptr = &(vp.data()[(tx-1)*np2 + ty-1]);
  const Vec3T<double> v11 = vptr[0];
  const Vec3T<double> v12 = vptr[1];
  const Vec3T<double> v13 = vptr[2];
  const Vec3T<double> v14 = vptr[3];
  
  const Vec3T<double> v21 = vptr[np2];
  const Vec3T<double> v22 = vptr[np2+1];
  const Vec3T<double> v23 = vptr[np2+2];
  const Vec3T<double> v24 = vptr[np2+3];

  const Vec3T<double> v31 = vptr[2*np2];
  const Vec3T<double> v32 = vptr[2*np2+1];
  const Vec3T<double> v33 = vptr[2*np2+2];
  const Vec3T<double> v34 = vptr[2*np2+3];
  
  const Vec3T<double> v41  =vptr[3*np2];
  const Vec3T<double> v42 = vptr[3*np2+1];
  const Vec3T<double> v43 = vptr[3*np2+2];
  const Vec3T<double> v44 = vptr[3*np2+3];
  
  ret[0] = Point3(0.0);
  for(int j=0; j<=3; j++)
    for (int k=0; k<=3; k++)
      ret[0] +=  Nj[j]*Nk[k]*vptr[j*np2+k];
  
  
  ret[1] =  Njd[0]*((v21-v11)*Nk[0]+ (v22-v12)*Nk[1]+ (v23-v13)*Nk[2]+ (v24-v14)*Nk[3]);
  ret[1] += Njd[1]*((v31-v21)*Nk[0]+ (v32-v22)*Nk[1]+ (v33-v23)*Nk[2]+ (v34-v24)*Nk[3]);    
  ret[1] += Njd[2]*((v41-v31)*Nk[0]+ (v42-v32)*Nk[1]+ (v43-v33)*Nk[2]+ (v44-v34)*Nk[3]);
  
  ret[2] = Nj[0]*((v12-v11)*Nkd[0]+ (v13-v12)*Nkd[1]+ (v14-v13)*Nkd[2]);
  ret[2] += Nj[1]*((v22-v21)*Nkd[0]+ (v23-v22)*Nkd[1]+ (v24-v23)*Nkd[2]);
  ret[2] += Nj[2]*((v32-v31)*Nkd[0]+ (v33-v32)*Nkd[1]+ (v34-v33)*Nkd[2]);    
  ret[2] += Nj[3]*((v42-v41)*Nkd[0]+ (v43-v42)*Nkd[1]+ (v44-v43)*Nkd[2]);
  
  return 0;
} 



// ---------------------------------------------------------------------- 
int BdSurf::evalBSpline(int flags, int V, double* xy, Point3* ret){ 
  assert( (flags==EVAL_VALUE) || (flags==(EVAL_VALUE|EVAL_1ST_DERIV)));


  /**** DZ ****/
 if(::LoadW) { 
    Point3* cur = ret;  
    //    if(V == 0) cout << "xy: [" << xy[0] << "," << xy[1] << "],";
    evalSplineFromCP(
		     ::flags2chartId(gpmesh().boun(V),gpmesh().corner(V)), _splineDeg, valence(V), xy, *_controls[V], cur);			 
    cur++; 
    if (flags &EVAL_1ST_DERIV) 
    evalSplineDerFromCP(
			::flags2chartId(gpmesh().boun(V),gpmesh().corner(V)), _splineDeg, valence(V), xy, *_controls[V], cur);			 
    return 0; 
  }
  /*** end DZ ***/
  
  Vector<Point3>* ivp = _controls[V];
  assert( ivp!=NULL );
  Vector<Point3>& vp = *ivp;
  double x = xy[0];	 double y = xy[1];
  
  int cn = 0;
  int cm =0; int nrgrid =0;
  double xlimit = 0; double ylimit = 0; 
  double xmin = 0; double ymin = 0; 
  double totalx = 0; double totaly = 0;
  int degk =0 ; int ex = 0;
  initSplineVars(V, degk, cn, cm, xlimit, ylimit, totalx, totaly, xmin, ymin, ex, nrgrid);
  
  int gridpts= nrgrid+2*ex;
  if(degk%2 == 0 && (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX || gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX))
    gridpts--;
    
  Point3* cur = ret;
  
  double newx = (x-xmin)*(nrgrid-1)/totalx;
  double newy = (y-ymin)*(nrgrid-1)/totaly;
  
  //find relevant grid pts to compute basis functions for..
  int tx = (int)ceil(newx); int ty = (int)ceil(newy);   
  if(ty == 0) ty++;  
  if(tx == 0) tx++;
  double *Nj = new double[gridpts];
  double *Nk = new double[gridpts];	  
  SplineBasis(newx, tx, degk, gridpts, Nj);
  SplineBasis(newy, ty, degk, gridpts, Nk); 
  double * Njder = new double[gridpts-1];
  double * Nkder = new double[gridpts-1];   
  SplineBasis(newx, tx, degk-1,  gridpts-1,Njder);
  SplineBasis(newy, ty, degk-1,  gridpts-1,Nkder);      
  int maxx = tx+ex; int maxy = ty+ex;
  if(degk%2 == 0 && (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX || gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX)){
    maxx--; maxy--;
  }  
  if (flags &EVAL_VALUE){   
    *cur = Point3(0.0);
    for (int j=tx-(ex+1); j<=maxx; j++){
      for (int k=ty-(ex+1); k<=maxy; k++){
	(*cur) += Nj[j+ex] * Nk[k+ex] *vp((j+ex)*(gridpts) + (k+ex));
      }
    }
    cur = cur+1;
  }
  
  if (flags &EVAL_1ST_DERIV) {
    double scale, tj, tk;
    {	
      *cur = Point3(0.0);
      for (int j=tx-(ex+1); j<maxx; j++){
	tj= xmin+j+1;
	scale = 1.0/3.0;
	for (int k=ty-(ex+1); k<=maxy; k++){
	  (*cur) += (vp((j+1+ex)*(gridpts)+(k+ex)) - vp((j+ex)*(gridpts)+(k+ex)))*scale*Njder[j+ex]*Nk[k+ex]; 
	}
      }
      cur = cur+1;	
    }
    {
      *cur = Point3(0.0);
      for (int j=tx-(ex+1); j<=maxx; j++){
	for (int k=ty-(ex+1); k<maxy; k++){
	  tk=ymin+(k+1);
	  scale = 1.0/3.0;
	  (*cur)+= (vp((j+ex)*(gridpts)+(k+ex+1)) - vp((j+ex)*(gridpts)+(k+ex)))*scale*Nj[j+ex]*Nkder[k+ex];
	}	  
      }
      cur = cur+1;
    }
  }
  
  if(flags & EVAL_2ND_DERIV) {
    cerr <<"bdsurf, splines, second derivative not implemented"<<endl;
    assert(0);
  }
  
  delete[] Nj;     delete[] Nk;
  delete[] Njder;  delete[] Nkder;
  
  return 0;
}

// ---------------------------------------------------------------------- 
void BdSurf::gatherData(int V, vector<Point2>& tmpcd, vector<Point2>& tmpxy, 
                        vector<Point3>& tmpval){
   
  set<intpair> Vijset;
  int ctrlsize = pow2(_subdivCtrlLevel);  assert(_subdivCtrlLevel==2);
  int ss = (int)ceil( _evalUpperBnd * (ctrlsize) ); assert(ss==ctrlsize); 
  int K = valence(V);
  double ctrlstep = 1.0/ctrlsize;
  
  if (_basisType == SPLINE_BASIS)  ss = ss+1; 
  
  for(int f=0; f<K; f++) {
    for(int j=0; j<ss; j++)
      for(int i=0; i<ss; i++) {
	//----------remove redundancy
	if(i==0 && j==0 && f!=0) continue; 
	if(i==0 && j>0) continue;
	//----------get position in different parameterization
	
	double cd[2];             
	cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
	double xy[2];    
       
	Vfcd2Vxy(EVAL_VALUE, V, f, cd, xy);     //C position

        int F;
	double st[2]; 
	Vfcd2Fcd(EVAL_VALUE, V, f, cd, F, st); //F position
	int gh[2];             
	gh[0] = (int)round(st[0]*ctrlsize);            
	gh[1] = (int)round(st[1]*ctrlsize);
    	
	int MUL = pow2(20);
	intpair tVij( (int)round(xy[0]*MUL), (int)round(xy[1]*MUL) );
	if(Vijset.find(tVij)==Vijset.end()) { //not there yet
	  Vijset.insert(tVij);        //mark this point as processed
	  //1. pos
	  tmpcd.push_back( Point2(cd) );
	  tmpxy.push_back( Point2(xy) );   
	  tmpval.push_back(_ccLimitPos[F](gh[0],gh[1]));
	  //          if(V == 17 && _ccLimitPos[F](gh[0], gh[1])(0)<2+1e-5 && 
	  //             _ccLimitPos[F](gh[0], gh[1])(0)>2-1e-5)
	  //              cout<<xy[1]<<" "<<_ccLimitPos[F](gh[0], gh[1])(2)<<endl;
              //cout<<"V = "<<V<<" -- "<<xy[0]<<" "<<xy[1]<<" -- "<<_ccLimitPos[F](gh[0], gh[1])<<endl;
	} 
	else {               //REDUDANCY IS ALREADY REMOVED
	  assert(0);
	}
      }
  }
  
  //if on boundary need (ctrlsize-1) extra points.
  if (gpmesh().boun(V) == GpMesh::BOUNDARY_VERTEX){      
    
    int f = K-1;
    int i=0;
    for (int j=1; j<ss; j++){    
      double cd[2];          
      cd[0] = i*ctrlstep;  cd[1] = j*ctrlstep; //V position
      double xy[2];    
      Vfcd2Vxy(EVAL_VALUE, V, f, cd, xy);     //C position
      int F;
      double st[2];
      Vfcd2Fcd(EVAL_VALUE, V, f, cd, F, st); //F position
      
      int gh[2];            
      gh[0] = (int)round(st[0]*ctrlsize);            
      gh[1] = (int)round(st[1]*ctrlsize);
      
      int MUL = pow2(20);
      intpair tVij( (int)round(xy[0]*MUL), (int)round(xy[1]*MUL) );
      if(Vijset.find(tVij)==Vijset.end()) { //not there yet
	Vijset.insert(tVij);         //mark this point as processed
	//1. pos
	tmpcd.push_back( Point2(cd) );
	tmpxy.push_back( Point2(xy) );
	tmpval.push_back(_ccLimitPos[F](gh[0],gh[1]));
      } 
      else {              //REDUDANCY IS ALREADY REMOVED
	assert(0);
      }
    }
  }
}


// ---------------------------------------------------------------------- 
void BdSurf::preProcessPoly(int K, vector<Point3>& tmpval){
  if(_basisType == POLYNOMIAL_BASIS || _basisType==5) {
    cerr<<"PREPROCESSING"<<endl;
    //do projection
    BdULib::Entry& E = _matricesU->chooseEntry(K);
    NumMatrix& UF = E.U(); //full matrix
    NumMatrix& IUF = E.IU(); //full inverse matrix
    assert(UF.m()==tmpval.size() && UF.n()==6*K+1);
    assert(IUF.m()==6*K+1&& IUF.n()==tmpval.size());
    
    //LEXING: change based on the lib file
    int NUMCHK = tmpval.size();
    int NUMSIG = 6*K+1; 
    
    vector<Point3> tmp(NUMSIG, Point3(0,0,0));
    for(int i=0; i<NUMSIG; i++)
      for(int j=0; j<NUMCHK; j++)
	tmp[i] += IUF(i,j) * tmpval[j];
    
    for(int i=0; i<NUMCHK; i++) {
      tmpval[i] = Point3(0,0,0);
      for(int j=0; j<NUMSIG; j++)
	tmpval[i] += UF(i,j) * tmp[j];
    }
  }
}    

// ---------------------------------------------------------------------- 
void BdSurf::constructControls(int numv){
  DataLevel  ctrldu;
  DataLevel  ctrldv;
  DataLevel& ctrlpos = _ccSurf.pos(_subdivCtrlLevel); //control pos
  
  CCSurfOp::limit(gpmesh(), _subdivCtrlLevel, ctrlpos, _ccLimitPos, *_subdivMatrices);
  CCSurfOp::du(   gpmesh(), _subdivCtrlLevel, ctrlpos, ctrldu);
  CCSurfOp::dv(   gpmesh(), _subdivCtrlLevel, ctrlpos, ctrldv);
  
  _controls.resize(numv, NULL);
}

// ---------------------------------------------------------------------- 
void BdSurf::constructPoly(int V, vector<Point2>& tmpxy, 
                           vector<Point3>& tmpval){
  int cm = tmpxy.size();
  int K = valence(V);
  int DEG;
  if (K>=12 && _max_degree_valence != 0) DEG = _max_degree_valence+2; //cutoff to avoid high degree polynomials
  else if (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX && K>=6 ) DEG = K; //boundary special case(prob in cylinder
  else if (gpmesh().boun(V) == GpMesh::BOUNDARY_VERTEX && K==2 ) DEG = K+4; //boundary special case(prob in cylinder
  else if (K==1) DEG = K+3;
  else DEG = K+2;

  int cn = (DEG+1)*DEG/2; 
  
  if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX || !_indepBoundaryFit) //INTERIOR OR !INDEP
    constructPolyInterior( V, DEG, cm, cn, tmpxy, tmpval);
  else if (gpmesh().corner(V) == GpMesh::CREASE_VERTEX) //BOUNDARY    
    constructPolyBoundary( V, DEG, cm, cn, tmpxy, tmpval);
  else // CONVEX || CONCAVE 
    constructPolyCorner( V, DEG, cm, cn, tmpxy, tmpval); 
    
}

// ---------------------------------------------------------------------- 
void BdSurf::constructPolyInterior( int V, int DEG, int cm, int cn, 
                                    vector<Point2>& tmpxy, 
                                    vector<Point3>& tmpval){
  NumMatrix M(cm,cn); clear(M);
  for(int i=0; i<cm; i++) {
    double xs[50], ys[50];
    initPolyXsYs(DEG, tmpxy[i](0), tmpxy[i](1), &xs[0], &ys[0]);
    //----------------
    //M stores the evaluation of all monomials of degree<=K+1
    //at all 12k+1 ctrl nodes.
    int j = 0; 
    for(int d=0; d<DEG; d++) {
      for(int t=0; t<=d; t++)               
	M(i,j++) = xs[d-t] * ys[t];
    }
    assert(j==cn);
  }
  NumMatrix piM(cn,cm);
  //inverse of M (piM) is the coefficient matrix.
  pinv(M, 1e-10, piM);
  
  updatePolyControlPoint(V,/* cn, cm,*/ piM, tmpval);
  
}

// ---------------------------------------------------------------------- 
void BdSurf::constructPolyBoundary(int V, int DEG, int cm, int cn, 
                                   vector<Point2>& tmpxy, 
                                   vector<Point3>& tmpval){

  int ctrlsize = pow2(_subdivCtrlLevel);  assert(_subdivCtrlLevel==2);
  //first gather all boundary points in bounxy.
  vector<Point2> bounxy;
  for(int i=0; i<ctrlsize; i++)
    bounxy.push_back(tmpxy[i]);
  for (int i=cm-ctrlsize+1; i<cm; i++)
    bounxy.push_back(tmpxy[i]);
  
  
  //now, have all boundary points in bounxy.
  assert(bounxy.size() == 2*(ctrlsize-1)+1);
  
  int bpts = bounxy.size();    //nr of boundary nodes
  int ipts = cm-bounxy.size(); //nr of interior nodes
  int umon = DEG;              //nr of monomials with u only
  int uvmon =cn-DEG;           //nr of leftover monomials.
  
  //first solve least squares for the boundary nodes
  
  //B3 stores each 1-variable(u) monomial evaluated
  //at the boundary nodes
  NumMatrix B3(bpts,umon); clear(B3);
  for(int i=0; i<bpts; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, bounxy[i](0), bounxy[i](1), &xs[0], &ys[0]);    
    int j = 0;  
    for(int d=0; d<DEG; d++)  B3(i,j++) = xs[d];
    assert(j==umon);
  }    
  //M3 stores the pseudo inverse. the matrix of coefficients.
  NumMatrix M3(umon,bpts);    pinv(B3, 1e-10, M3);
  
  //B1 has "leftover" monomials evaluated at interior nodes
  NumMatrix B1(ipts, uvmon); clear(B1);
  for(int i=ctrlsize; i<=cm-ctrlsize; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, tmpxy[i](0), tmpxy[i](1), &xs[0], &ys[0]);    
    int j = 0; 
    for(int d=0; d<DEG; d++) {
      for(int t=1; t<=d; t++) 
	B1(i-ctrlsize,j++) = xs[d-t] * ys[t];
    }
    assert(j==uvmon);

  }
  
  //inverse of this gives the coefficient matrix.
  NumMatrix M1(uvmon,ipts);      pinv(B1, 1e-10, M1);
  
  //Compute M2
  //matrix for boundary node "leftover"  monomials
  //this is given by: M2 = -M1*B2*pinv(B1) = -M1*B2*M3
  
  NumMatrix B2(ipts, umon); clear(B2);
  for(int i=ctrlsize; i<=cm-ctrlsize; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, tmpxy[i](0), tmpxy[i](1), &xs[0], &ys[0]);
    int j = 0;  
    for(int d=0; d<DEG; d++) 
      B2(i-ctrlsize,j++) = xs[d];
    assert(j==umon);

  }
  multvalue<double>(M1,-1);  
  NumMatrix M2(uvmon, bpts);
  NumMatrix T(uvmon, umon);
  NumMatrixMult(M1, B2, T);
  NumMatrixMult(T, M3,M2);
  
  //write the values in the main coefficient matrix.
  NumMatrix U(cn,cm); clear(U);
  //store values for boundary points on monomials with x only (M3).
  set<int> degs;
  for (int i=0;i<umon; i++){
    int row = i*(i+1)/2;
    degs.insert(row);
    for (int j=0; j<ctrlsize; j++)
      U(row, j) = M3(i, j);
    for (int j=ctrlsize; j<bpts; j++)
      U(row, cm-bpts+j) = M3(i, j);
  }
  //store M1  (interiors that doesnt have x only)
  int count = 0;
  for (int i=0; i<cn; i++){
    if(degs.find(i) == degs.end()){ // if not already filled.
      for (int j=ctrlsize; j<=cm-ctrlsize; j++)
	U(i, j) = -1*M1(count, j-ctrlsize); 
      count++;
    }
  }
  //store M2
  count =0;
  for (int i=0;i<cn; i++){
    if(degs.find(i) == degs.end()){ // if not already filled.
      for (int j=0; j<ctrlsize; j++)
	U(i, j) = M2(count, j);
      for (int j=ctrlsize; j<bpts; j++)
	U(i, cm-bpts+j) = M2(count, j);
      count++;
    }
  }
  updatePolyControlPoint(V, /*cn, cm,*/ U, tmpval);
}
// ---------------------------------------------------------------------- 
void BdSurf::constructPolyCorner(int V, int DEG, int cm, int cn, 
                                 vector<Point2>& tmpxy, 
                                 vector<Point3>& tmpval){

  int ctrlsize = pow2(_subdivCtrlLevel);  assert(_subdivCtrlLevel==2);
  //first gather all boundary points of corner in bounxy.
  vector<Point2> bounxy;
  for(int i=0; i<ctrlsize; i++)        bounxy.push_back(tmpxy[i]);
  for (int i=cm-ctrlsize+1; i<cm; i++) bounxy.push_back(tmpxy[i]);
  
  //now, have all boundary points in bounxy.
  assert(bounxy.size() == 2*(ctrlsize-1)+1);
  
  int bpts = bounxy.size();       //nr of boundary nodes
  int ipts = cm-bounxy.size();    //nr of interior nodes
  int uvmon = 2*DEG-1;            //nr of monomials with u only
  int other  =cn-uvmon;           //nr of leftover monomials.
  
  //first solve least squares for the boundary nodes
  
  //B3 stores each 1-variable(u) monomial evaluated
  //at the boundary nodes
  
  NumMatrix B3(bpts,uvmon); clear(B3);
 
 
  for(int i=0; i<ctrlsize; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, bounxy[i](0), bounxy[i](1), &xs[0], &ys[0]);    
    int j = 0; 
    for(int d=0; d<DEG; d++)  
      B3(i,j++) = xs[d];
    assert(j==DEG);
  }
  
  for(int i=ctrlsize; i<bpts; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, bounxy[i](0), bounxy[i](1), &xs[0], &ys[0]);
    B3(i,0) = ys[0];
    int j = DEG;  
    for(int d=1; d<DEG; d++)  
      B3(i,j++) = ys[d];
    assert(j==uvmon);
  }
  NumMatrix M3(uvmon,bpts);clear(M3);     pinv(B3, 1e-10, M3);

  //B1 has "leftover" monomials evaluated at interior nodes
  NumMatrix B1(ipts, other); clear(B1);
  
  for(int i=ctrlsize; i<=cm-ctrlsize; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, tmpxy[i](0), tmpxy[i](1), &xs[0], &ys[0]);    
    int j = 0;   
    for(int d=1; d<DEG; d++) 
      for(int t=1; t<d; t++)
	B1(i-ctrlsize,j++) = xs[d-t] * ys[t];    
    assert(j==other);
  }
  //inverse of this gives the coefficient matrix.
  NumMatrix M1(other,ipts);  clear(M1);    pinv(B1, 1e-10, M1);
  //Compute M2
  //matrix for boundary node "leftover"  monomials
  //this is given by: M2 = -M1*B2*pinv(B1) = -M1*B2*M3  
  NumMatrix B2(ipts, uvmon); clear(B2);
   for(int i=ctrlsize; i<=cm-ctrlsize; i++){
    //construct prep
    double xs[50], ys[50];
    initPolyXsYs(DEG, tmpxy[i](0), tmpxy[i](1), &xs[0], &ys[0]);    
    int j = 0;   
    for(int d=0; d<DEG; d++) 
      B2(i-ctrlsize,j++) = xs[d];
    assert(j==DEG);
    for (int d=DEG; d<uvmon; d++)
      B2(i-ctrlsize, j++) = ys[d-DEG+1];
    assert(j==uvmon);
  }
  
   NumMatrix M2(other, bpts); clear(M2);
   NumMatrix T(other, uvmon); clear(T);
  NumMatrixMult(M1, B2, T);
  multvalue<double>(T,-1);
  NumMatrixMult(T, M3,M2);
  
  //write the values in the main coefficient matrix.
  NumMatrix U(cn,cm); clear(U);
  
  //store values for boundary points on monomials with x only (M3).
  set<int> degs;
  for (int i=0;i<DEG; i++){
    int row = i*(i+1)/2;
    degs.insert(row);
        for (int j=0; j<ctrlsize; j++)
      U(row, j) = M3(i, j);
      }
  for (int i=2;i<=DEG; i++){
    int row = i*(i+1)/2-1;
    degs.insert(row);
    U(row, 0) = M3(i-2+DEG, 0);
    for (int j=ctrlsize; j<bpts; j++)
      U(row, cm-bpts+j) = M3(i-2+DEG, j);
    }
  //store M1  (interiors that doesnt have x only)
  int count = 0;
  for (int i=0; i<cn; i++){
    if(degs.find(i) == degs.end()){ // if not already filled.
      for (int j=ctrlsize; j<=cm-ctrlsize; j++)
	U(i, j) = M1(count, j-ctrlsize); 
      count++;
    }
  }
  //store M2
  count =0;
  for (int i=0;i<cn; i++){
    if(degs.find(i) == degs.end()){ // if not already filled.
      for (int j=0; j<ctrlsize; j++)
	U(i, j) = M2(count, j);
      for (int j=ctrlsize; j<bpts; j++)
	U(i, cm-bpts+j) = M2(count, j);
      count++;
    }
  }
  updatePolyControlPoint(V, /*cn, cm,*/ U, tmpval);
  
}

// ---------------------------------------------------------------------- 
void BdSurf::initPolyXsYs(int n, double x, double y, double* xs, double* ys){
  xs[0] = 1; ys[0] = 1;
  for(int i=1; i<n; i++){
    xs[i] = xs[i-1] * x;
    ys[i] = ys[i-1] * y;
  }
}
// ---------------------------------------------------------------------- 
void BdSurf::updatePolyControlPoint(int V, /*int cn, int cm, */NumMatrix& M, vector<Point3>& Val){
 
  int cn = M.m();
  int cm = Val.size();
 

  assert(cm == M.n());
  
  Vector<Point3>* vpptr = new Vector<Point3>(cn);
  for(int i=0; i<cn; i++) {
    Point3& cur = (*vpptr)(i);
    cur = Point3(0.0);
    for(int j=0; j<cm; j++)    
      cur += M(i,j) * Val[j];
  }
  _controls[V] = vpptr;
}




// ---------------------------------------------------------------------- 
void BdSurf::constructSpline(int V, vector<Point2>& tmpxy, vector<Point3>& tmpval){
  
  int ctrlsize = pow2(_subdivCtrlLevel);  assert(_subdivCtrlLevel==2);
  int bpts = ctrlsize+1;
  cerr << "setup for vertex " << V << endl;


  
  if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX  || !_indepBoundaryFit) //interior or not indep
    constructSplineInterior(V, tmpxy, tmpval);
  else if( _indepBoundaryFit && gpmesh().boun(V) == GpMesh::BOUNDARY_VERTEX){ //if indep boundary
    if(gpmesh().corner(V) == GpMesh::CREASE_VERTEX) // smooth boun
      constructSplineIndepSmooth(V, bpts, tmpxy, tmpval);
    else if (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX ) //convex corner
      constructSplineIndepConvex(V, bpts, tmpxy, tmpval);
    else if (gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX ) //concave corner
      constructSplineIndepConcave(V, bpts, tmpxy, tmpval);  
  }

  /*** DZ ***/
    constructSplineCPfromW( ::flags2chartId(gpmesh().boun(V),gpmesh().corner(V)), _splineDeg, valence(V),  tmpval);
  if( ::LoadW) { _controls[V] = constructSplineCPfromW( ::flags2chartId(gpmesh().boun(V),gpmesh().corner(V)), _splineDeg, valence(V),  tmpval); return; }
  /*** DZ ***/
}
// ---------------------------------------------------------------------- 
inline int BdSurf::SplineBasis(double u, int k, int p, int n, double *N){ 
  assert(k>0);
  

  double x = u-k+1;
  double x2 = x*x;
  double x3 = x2*x;

  for(int i=0; i<n; i++) N[i] = 0.0;
  
  if (p==2){
    N[k-1] = (1.0+x2-2.0*x)/2.0;
    N[k]   = (1.0-2.0*x2+2.0*x)/2.0;
    N[k+1] = x2/(2.0);
    return(0);
  }
  else if (p==3){
    N[k-1] = (-x3+3.0*x2-3.0*x+1.0)*(1.0/6.0);
    N[k]   = (4.0-6.0*x2+3.0*x3)*(1.0/6.0);
    N[k+1] = (1.0+3.0*x+3.0*x2-3.0*x3)*(1/6.0);
    N[k+2] = x3*(1/6.0);
    return(0);
  }
  else if (p==4){
    double x4 = x2*x2;
    N[k-1] = (x4-4.0*x3+6.0*x2-4.0*x+1.0)/24.0;
    N[k] =   (11.0-12.0*x-6.0*x2+12.0*x3-4.0*x4)/24.0;
    N[k+1] = (11.0+12.0*x-6.0*x2-12.0*x3+6.0*x4)/24.0;
    N[k+2] = (1.0+4.0*x+6.0*x2+4.0*x3-4.0*x4)/24.0;
    N[k+3] = x4/24.0;
    return(0);
  }
  else {
    k--;
    N[k+p] = 1.0;
    
    for (int d=1; d<=p; d++){
        N[k-d+p] = ((k+1-u)/d)*N[k-d+1+p];
        for (int i=k-d+1; i<=k-1; i++)
            N[i+p] = ((u-i)/d)*N[i+p] + ((i+d+1-u)/d)*N[i+1+p];
        N[k+p] = ((u-k)/d)*N[k+p];
    }
  }
  return 0;
}


// ---------------------------------------------------------------------- 
void BdSurf::initSplineVars(int V, int &degk, int &cn, int &cm, double &xlimit, 
			    double &ylimit, double &totalx, double &totaly,
			    double &xmin, double &ymin, int &ex, int &nrgrid){
  
  //initialize variables related to spline fit
  int K = valence(V);  
  degk = _splineDeg;
  ex = (int)floor((degk)/2.0);  
  double scale = 0.0;
  int newdeg = 3;
  nrgrid = 3;

  if (gpmesh().boun(V) == GpMesh::INTERIOR_VERTEX){ // interior
    cn = 20*K+1; xlimit = ylimit = 1.8;
    xmin = ymin = -1*xlimit;
    totalx = totaly = 2*xlimit;
    scale =  2.0/3.0;
    newdeg = 5;
  }
  else if (gpmesh().boun(V) == GpMesh::BOUNDARY_VERTEX){
    scale = .5;
    newdeg = 4; 
    if (gpmesh().corner(V) == GpMesh::CREASE_VERTEX){ //boundary
      if(K==2) {
	scale = 2.0/3.0;
        newdeg = 4;
      }
      cn = 20*K+5; xlimit = ylimit = 1.8;
      xmin = -1*xlimit;   ymin = 0.0;
      totalx = 2*xlimit;  totaly = ylimit-ymin;
    }
    else if (gpmesh().corner(V) == GpMesh::CONVEX_VERTEX){ //convex 
      newdeg = 3;
      cn = 20*K+5; xlimit = ylimit = 1.8;
      xmin = ymin = 0.0;
      totalx = xlimit-xmin;  totaly = ylimit-ymin;
    }
    else if (gpmesh().corner(V) == GpMesh::CONCAVE_VERTEX){ //concave
      newdeg = 3;
      cn = 20*K+5; xlimit = ylimit = 1.8;
      xmin = ymin = -1*xlimit;
      totalx = totaly = 2*xlimit;
    }
  }

  cm = (nrgrid+(2*ex))*(nrgrid+(2*ex));
  if (cm >= scale*cn){
    nrgrid = 2;
    degk = newdeg;
    ex =(int) floor((degk)/2.0) ;
    cm = (nrgrid+(2*ex))*(nrgrid+(2*ex));
  }

  assert(cm < cn);
}
// ---------------------------------------------------------------------- 

void BdSurf::constructSplineInterior(int V, vector<Point2>& tmpxy, 
				     vector<Point3>& tmpval){
  
  int cn = 0; int cm =0;
  double xlimit = 0; double ylimit = 0; //rightmost(topmost) x(y) coordinate among all knots
  double totalx = 0; double totaly = 0; //length of each side of rect surrounding grid
  double xmin = 0; double ymin = 0;     //min values of x and y for knots
  int ex =0;  int nrgrid = 0; int degk = 0;
  initSplineVars(V, degk, cn, cm, xlimit, ylimit, totalx, totaly, xmin, ymin, ex, nrgrid);
  
  
  NumMatrix F(cn, cm); clear(F);  
  double newx, newy;
  int tx, ty;  
  double *Nj = new double[nrgrid+2*ex];
  double *Nk = new double[nrgrid+2*ex];

  vector<Point3> temp; //right hand side of eqn
  
  for (int i=0; i<cn; i++){
    newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
    newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
    //find relevant grid points to compute basis functions for.
    tx = (int)ceil(newx);  ty = (int)ceil(newy);
    if(ty==0) ty++;
    if(tx==0) tx++;
    SplineBasis(newx, tx, degk, nrgrid+2*ex, Nj);
    SplineBasis(newy, ty, degk, nrgrid+2*ex, Nk);   
    for (int j=tx-(ex+1); j<=tx+ex; j++){
      for (int k=ty-(ex+1); k<=ty+ex; k++)
	F(i, (j+ex)*(nrgrid+2*ex)+(k+ex)) = Nj[j+ex]*Nk[k+ex];      
    }

    temp.push_back(tmpval[i]);
  }
  //take inverse and multiply with  ctrl pts to get spline coefficients
  NumMatrix piF(cm,cn);      pinv(F, 1e-16, piF);
  
  //multiply with matrix to get control points
  Vector<Point3>* vpptr = new Vector<Point3>(cm);
  for(int i=0; i<cm; i++) {
    Point3& cur = (*vpptr)(i);
    cur = Point3(0.0);
    for(int j=0; j<cn; j++)
      cur += piF(i,j) * temp[j];
  }
  _controls[V] = vpptr;
  delete[] Nj; delete[] Nk;
}

// ---------------------------------------------------------------------- 
void BdSurf::constructSplineIndepSmooth(int V, int bpts,vector<Point2>& tmpxy, 
					vector<Point3>& tmpval){
  int cn = 0; int cm =0;
  double xlimit = 0; double ylimit = 0; //rightmost(topmost) x(y) coordinate among all knots
  double totalx = 0; double totaly = 0; //length of each side of rect surrounding grid
  double xmin = 0; double ymin = 0;     //min values of x and y for knots
  int ex =0;  int nrgrid = 0; int degk = 0;
  initSplineVars(V, degk, cn, cm, xlimit, ylimit, totalx, totaly, xmin, ymin, ex, nrgrid);
  vector<Point3> bcoef;  bcoef.resize(nrgrid+2*ex);

  computeSpBounCoefSmooth(degk, nrgrid, bpts, cn, totalx, xmin, tmpxy, tmpval, bcoef);

  int nre = nrgrid+2*ex; //nr of equations
  
  NumMatrix B1(nre, cm-nre); clear(B1);
  NumMatrix B2(nre, nre); clear(B2);
  NumMatrix A1(cn, cm-nre); clear(A1);
  NumMatrix A2(cn, nre); clear(A2);

  createSpIndepMatrices(BdSurf::SMOOTH_MAT,degk, nrgrid, tmpxy, 
		      xmin, ymin, totalx, totaly,A1, A2, B1, B2); 
 
  vector<Point3> p1;  p1.resize(cm-nre);  
  vector<Point3> p2;  p2.resize(nre);
  
  computeSpIndepControls(B1, B2, A1, A2, bcoef, tmpval, p1, p2);
  
  Vector<Point3>* vpptr = new Vector<Point3>(cm);
  
  for(int j=-ex; j<nrgrid+ex; j++)
    (*vpptr)((j+ex)*(nrgrid+2*ex) + (0+ex)) = p2[j+ex]; //k =0;
  for(int j=-ex; j<nrgrid+ex; j++){
    for(int k=-ex; k<0; k++)
      (*vpptr)((j+ex)*(nrgrid+2*ex) + (k+ex)) = p1[(j+ex)*(nre-1)+(k+ex)];
    for(int k=1; k<nrgrid+ex; k++)
      (*vpptr)((j+ex)*(nrgrid+2*ex) + (k+ex)) = p1[(j+ex)*(nre-1)+(k-1+ex)];
  }	    
  _controls[V] = vpptr;  
}
// ---------------------------------------------------------------------- 
void BdSurf::constructSplineIndepConvex(int V, int bpts,vector<Point2>& tmpxy, 
					vector<Point3>& tmpval){
  
  //indep bound, convex corner

  int cn = 0; int cm =0;
  double xlimit = 0; double ylimit = 0; //rightmost(topmost) x(y) coordinate among all knots
  double totalx = 0; double totaly = 0; //length of each side of rect surrounding grid
  double xmin = 0; double ymin = 0;     //min values of x and y for knots
  int ex =0;  int nrgrid = 0; int degk = 0;
  initSplineVars(V, degk, cn, cm, xlimit, ylimit, totalx, totaly, xmin, ymin, ex, nrgrid);

  int gridpts = nrgrid+2*ex;
  if(degk%2 == 0) gridpts = nrgrid+2*ex-1;
  int nre = 2*gridpts-1; //nr of equations 
  vector<Point3> bcoef; bcoef.resize(nre);
  computeSpBounCoefConvex(degk, nrgrid, bpts, cn, totalx, xmin, totaly,ymin,
			tmpxy, tmpval, bcoef);  

  
  //now fit surface to boundary  
  NumMatrix B1(nre, cm-nre); clear(B1);
  NumMatrix B2(nre, nre);    clear(B2);
  NumMatrix A1(cn, cm-nre); clear(A1);
  NumMatrix A2(cn, nre);    clear(A2);

  createSpIndepMatrices(BdSurf::CONVEX_MAT,degk, nrgrid, tmpxy, 
		      xmin, ymin, totalx, totaly,A1, A2, B1, B2); 

  vector<Point3> p1;  p1.resize(cm-nre);  
  vector<Point3> p2;  p2.resize(nre);
  
  computeSpIndepControls(B1, B2, A1, A2, bcoef, tmpval, p1, p2);
  
  Vector<Point3>* vpptr = new Vector<Point3>(cm);
  
  
  for(int j=-ex; j<gridpts-ex; j++)
    (*vpptr)((j+ex)*(gridpts) + (0+ex)) = p2[j+ex];
  for(int k=-ex; k<0; k++)
    (*vpptr)((0+ex)*(gridpts) + (k+ex)) = p2[gridpts+ex+k];
  for(int k=1; k<gridpts-ex;k++)
    (*vpptr)((0+ex)*(gridpts) + (k+ex)) = p2[gridpts+ex+k-1];
  
  for(int j=-ex; j<0; j++){
    for(int k=-ex; k<0; k++)
      (*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j+ex)*(gridpts-1)+(k+ex)];
    for(int k=1; k<gridpts-ex; k++)
      (*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j+ex)*(gridpts-1)+(k-1+ex)];
  }        
  
  for(int j=1; j<gridpts-ex; j++){
    for(int k=-ex; k<0; k++)
      (*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j-1+ex)*(gridpts-1)+(k+ex)];
    for(int k=1; k<gridpts-ex; k++)
      (*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j-1+ex)*(gridpts-1)+(k-1+ex)];
  }
  
  _controls[V] = vpptr;  
}
									
// ---------------------------------------------------------------------- 
void BdSurf::constructSplineIndepConcave(int V, int bpts,vector<Point2>& tmpxy, 
					vector<Point3>& tmpval){
  int cn = 0; int cm =0;
  double xlimit = 0; double ylimit = 0; //rightmost(topmost) x(y) coordinate among all knots
  double totalx = 0; double totaly = 0; //length of each side of rect surrounding grid
  double xmin = 0; double ymin = 0;     //min values of x and y for knots
  int ex =0;  int nrgrid = 0; int degk = 0;
  initSplineVars(V, degk, cn, cm, xlimit, ylimit, totalx, totaly, xmin, ymin, ex, nrgrid);
 
  int gridpts = nrgrid+2*ex;
  if(degk%2 == 0)  gridpts--;
  int nre; //nr of equations 
  if(nrgrid == 2) nre = 2*gridpts-1; 
  else if(nrgrid == 3) nre = (gridpts-1) + (gridpts-2); 

  vector<Point3> bcoef; bcoef.resize(nre);
  computeSpBounCoefConcave(degk, nrgrid, bpts, cn, totalx, xmin, totaly,ymin,
			tmpxy, tmpval, bcoef);  

  NumMatrix B2(nre, nre);    clear(B2);
  NumMatrix B1(nre, cm-nre); clear(B1);
  NumMatrix A1(cn, cm-nre); clear(A1);
  NumMatrix A2(cn, nre);    clear(A2);    

  if(nrgrid == 2){
   createSpIndepMatrices(BdSurf::CONCAVE2_MAT,degk, nrgrid, tmpxy, 
		      xmin, ymin, totalx, totaly,A1, A2, B1, B2); 
   vector<Point3> p1;  p1.resize(cm-nre);  
   vector<Point3> p2;  p2.resize(nre);

   
   computeSpIndepControls(B1, B2, A1, A2, bcoef, tmpval, p1,p2);  
   
   Vector<Point3>* vpptr = new Vector<Point3>(cm);
   for(int j=-ex; j<gridpts-ex; j++)
     (*vpptr)((j+ex)*(gridpts) + (-ex+ex)) = p2[j+ex];
   for(int k=-ex+1; k<gridpts-ex; k++)
     (*vpptr)((gridpts-ex-1+ex)*(gridpts) + (k+ex)) = p2[gridpts+ex+k-1];
   for(int j=-ex; j<gridpts-ex-1; j++){
     for(int k=-ex+1; k<gridpts-ex; k++)
       (*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j+ex)*(gridpts-1)+(k-1+ex)];
   }  
 
   
   _controls[V] = vpptr;  
  }
  
  else if (nrgrid==3){
    createSpIndepMatrices(BdSurf::CONCAVE3_MAT,degk, nrgrid, tmpxy, 
		       xmin, ymin, totalx, totaly,A1, A2, B1, B2); 
    vector<Point3> p1;  p1.resize(cm-nre);  
    vector<Point3> p2;  p2.resize(nre);
    
    computeSpIndepControls(B1, B2, A1, A2, bcoef, tmpval, p1,p2);  
    
    Vector<Point3>* vpptr = new Vector<Point3>(cm);
    
    for(int j=-ex+1; j<gridpts-ex; j++)
      (*vpptr)((j+ex)*(gridpts) + (1+ex)) = p2[j+ex-1];
    for(int k=-ex; k<=0; k++)
      (*vpptr)((1+ex)*(gridpts) + (k+ex)) = p2[gridpts+ex+k-1];    
    for(int k=2; k<gridpts-ex-1;k++)
      (*vpptr)((1+ex)*(gridpts) + (k+ex)) = p2[gridpts+ex+k-2];
    
    for(int j=-ex; j<=0; j++){
      for(int k=-ex; k<=0; k++)
	(*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j+ex)*(gridpts-1)+(k+ex)];
      for(int k=2; k<gridpts-ex; k++)
	(*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j+ex)*(gridpts-1)+(k-1+ex)];
    }  
    
    for(int j=2; j<gridpts-ex; j++){
      for(int k=-ex; k<=0; k++)
	(*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j-1+ex)*(gridpts-1)+(k+ex)];
      for(int k=2; k<gridpts-ex; k++)
	(*vpptr)((j+ex)*(gridpts) + (k+ex)) = p1[(j-1+ex)*(gridpts-1)+(k-1+ex)];
      
    } 
    (*vpptr)((-ex+ex)*(gridpts) + (1+ex)) = p1[cm-nre-2];
    (*vpptr)((1+ex)*(gridpts) + (gridpts-ex-1+ex)) = p1[cm-nre-1];
    
    _controls[V] = vpptr;  
  }

}  // END INDEP HIGH ORDER Concave

// ---------------------------------------------------------------------- 
void BdSurf::computeSpBounCoefSmooth(int degk, int nrgrid, int bpts, int cn, double totalx, double xmin,
				   vector<Point2> &tmpxy, vector<Point3> &tmpval, vector<Point3> &bcoef){
  int r = 2*bpts-1; //number of rows
  int ex = (int)floor((degk)/2.0);  
  double newx; int tx;
  double *Nj=  new double[nrgrid+2*ex];
  
  if(degk+1 <= r ) {
    NumMatrix M(r, nrgrid+2*ex); clear(M);
    vector<Point3> temp;
    
    for(int i=0; i<bpts; i++){ //first 5 boundary points	    
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx==0) tx++;
      SplineBasis(newx, tx, degk, nrgrid+2*ex, Nj);
      for (int j=tx-(ex+1); j<=tx+ex; j++){
	M(i, j+ex) = Nj[j+ex];
      }	
      temp.push_back(tmpval[i]);
    }
    for (int i=cn-bpts+1; i<cn; i++){ //the other 4
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx=(int)ceil(newx); 
      if(tx==0) tx++;
      SplineBasis(newx, tx, degk,nrgrid+2*ex, Nj);	    
      for (int j=tx-(ex+1); j<=(tx+ex); j++){
	M(i-cn+r, j+ex) = Nj[j+ex];
      }
      temp.push_back(tmpval[i]);  
    }
    
    NumMatrix piM(nrgrid+2*ex,r);	pinv(M, 1e-16, piM); // take inverse
    //multiply with data pts to get boun spline coefficients stored in bcoef.    
    matvecmult3(piM, temp, bcoef);
  }
  
  else { //don't have enough data points
    vector<Point3> cubic; cubic.resize(4);
    computeSpCubicControlsSmooth(bpts, cn, nrgrid, totalx, xmin, tmpxy, tmpval, cubic);
    vector<Point3> highorder; highorder.resize(2+2*ex);
    elevateDegree(degk, cubic, highorder);
    assert(nrgrid==2);
  }
  
  delete[] Nj; 
}



// ---------------------------------------------------------------------- 
void BdSurf::computeSpBounCoefConvex(int degk, int nrgrid, int bpts, int cn, 
				   double totalx, double xmin, double totaly, double ymin,
				   vector<Point2> &tmpxy, vector<Point3> &tmpval, vector<Point3> &bcoef){
  double newx, newy;
  int tx, ty;
  int ex = (int)floor((degk)/2.0);    
  int gridpts = nrgrid+2*ex;
  if(degk%2 == 0) gridpts = nrgrid+2*ex-1;
  double *Nj = new double[gridpts];
  double *Nk = new double[gridpts];
  
  vector<Point3> bcoefx;    bcoefx.resize(gridpts);
  vector<Point3> bcoefy;    bcoefy.resize(gridpts);
  int r = bpts;
  if(degk+1 <= r  && degk%2 == 1) {
    NumMatrix Mx(r, gridpts-ex); clear(Mx);
    vector<Point3> tempx;   
    for(int i=0; i<bpts; i++){ //5 points on x-axis
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx == 0) tx++;
      SplineBasis(newx, tx, degk, gridpts,Nj);      
      Mx(i,0) = Nj[0+ex];
      for(int j=1; j<=ex; j++)
	Mx(i,0) += 2*Nj[-j+ex];
      for(int j=1; j<=ex; j++)
	Mx(i, j) = Nj[j+ex] - Nj[-j+ex];
      for (int j=ex+1; j<gridpts-ex; j++)
	Mx(i, j) = Nj[j+ex];
      tempx.push_back(tmpval[i]);  
    }
    NumMatrix piMx(gridpts-ex,r);  pinv(Mx, 1e-16, piMx); // take inverse
    
    //coefficients of the boundary spline stored in bcoef
    vector<Point3> tmp; tmp.resize(gridpts-ex);
    matvecmult3(piMx, tempx, tmp);
    for(int i=0; i<gridpts-ex; i++)
      bcoefx[i+ex] = tmp[i];
    
    //now compute coefficients on y-axis
    r = bpts; //number of rows
    NumMatrix My(r, gridpts-ex); clear(My);
    vector<Point3> tempy;
    
    int i=0;
    newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
    ty = (int)ceil(newy);
    if(ty == 0) ty++;
    SplineBasis(newy, ty, degk,gridpts, Nk);
    My(i,0) = Nk[0+ex];
    for(int j=1; j<=ex; j++)
      My(i,0) += 2*Nk[-j+ex];
    for(int j=1; j<=ex; j++)
      My(i, j) = Nk[j+ex] - Nk[-j+ex];
    for (int j=ex+1; j<gridpts-ex; j++)
      My(i, j) = Nk[j+ex];
    tempy.push_back(tmpval[i]);
    
    for(int i=cn-bpts+1; i<cn; i++){
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      ty = (int)ceil(newy);
      if(ty==0) ty++;
      SplineBasis(newy, ty, degk,gridpts, Nk);
      My(i-cn+bpts,0) = Nk[0+ex];
      for(int j=1; j<=ex; j++)
	My(i-cn+bpts,0) += 2*Nk[-j+ex];
      for(int j=1; j<=ex; j++)
	My(i-cn+bpts, j) = Nk[j+ex] - Nk[-j+ex];
      for (int j=ex+1; j<gridpts-ex; j++)
	My(i-cn+bpts, j) = Nk[j+ex];
      tempy.push_back(tmpval[i]);
    }
   
    NumMatrix piMy(gridpts-ex,r);  pinv(My, 1e-16, piMy); // take inverse
    matvecmult3(piMy, tempy, tmp);
    for(int i=0; i<gridpts-ex; i++){
      bcoefy[i+ex] = tmp[i];
    }
    
    Point3 mid = Point3(bcoefy[0+ex](0),bcoefx[0+ex](1),bcoefx[0+ex](2));  
    bcoefy[0+ex] = bcoefx[0+ex] = mid;
    for(int i=-ex; i<0; i++){
      bcoefx[i+ex] = 2.0*mid - bcoefx[-i+ex];
      bcoefy[i+ex] = 2.0*mid - bcoefy[-i+ex];
    }
  }   
  else{ //not enough data pts on boundary or even deg   
    assert(nrgrid==2 || nrgrid == 3);
    if(nrgrid == 2){
      vector<Point3> cubicx; cubicx.resize(4);
      vector<Point3> cubicy; cubicy.resize(4);
      
      computeSpCubicControlsConvex(bpts, cn, nrgrid, totalx, xmin, totaly, ymin, 
				 tmpxy, tmpval, cubicx, cubicy);   
      int g=2;
      if(degk%2 == 0) g = 1; 
      vector<Point3> highorderx; highorderx.resize(g+2*ex);
      elevateDegree(degk, cubicx, highorderx);      
      vector<Point3> highordery; highordery.resize(g+2*ex);
      elevateDegree(degk, cubicy, highordery);
      
      for(int i=0; i<gridpts; i++){
        bcoefx[i] = highorderx[i];
        bcoefy[i] = highordery[i];  
      }
    }
    
    else if (nrgrid == 3){
      vector<Point3> cubicx; cubicx.resize(5);
      vector<Point3> cubicy; cubicy.resize(5);
      
      computeSpCubicControlsConvex(bpts, cn, nrgrid, totalx, xmin, totaly, ymin, 
				 tmpxy, tmpval, cubicx, cubicy);      
      int g=2;
      if(degk%2 == 0) g = 1; 
      vector<Point3> highx1; highx1.resize(g+2*ex);
      vector<Point3> highx2; highx2.resize(g+2*ex);
      vector<Point3> t1; t1.resize(4);
      vector<Point3> t2; t2.resize(4);
      
      for(int i=0; i<4; i++) t1[i] = cubicx[i];
      for(int i=1; i<5; i++) t2[i-1] = cubicx[i];
      
      elevateDegree(degk, t1, highx1);
      elevateDegree(degk, t2, highx2);
      
      for(int i=0; i<5; i++)        bcoefx[i] = highx1[i];
      for(int i=5; i<gridpts; i++)  bcoefx[i] = highx2[i-1];
      
      vector<Point3> highy1; highy1.resize(g+2*ex);
      vector<Point3> highy2; highy2.resize(g+2*ex);
      
      for(int i=0; i<4; i++) t1[i] = cubicy[i];
      for(int i=1; i<5; i++) t2[i-1] = cubicy[i];
      
      elevateDegree(degk, t1, highy1);
      elevateDegree(degk, t2, highy2);
      
      for(int i=0; i<5; i++)        bcoefy[i] = highy1[i];
      for(int i=5; i<gridpts; i++)  bcoefy[i] = highy2[i-1];
    }
  }

  for(int i=0; i<gridpts; i++)
    bcoef[i] = bcoefx[i];
  for(int i=0; i<ex; i++)
    bcoef[i+gridpts] = bcoefy[i];     
  for(int i=ex+1; i<gridpts; i++)
    bcoef[i-1+gridpts] = bcoefy[i];
  
  delete[] Nj; delete[] Nk;
}

// ---------------------------------------------------------------------- 
void BdSurf::computeSpBounCoefConcave(int degk, int nrgrid, int bpts, int cn, 
				    double totalx, double xmin, double totaly,double ymin,
				    vector<Point2>&tmpxy,vector<Point3>&  tmpval,vector<Point3>& bcoef){
  int ex = (int)floor((degk)/2.0); 
  double newx, newy;
  int tx, ty;
  int gridpts = nrgrid+2*ex;
  if(degk%2 == 0)  gridpts--;

  double *Nj = new double[gridpts];
  double *Nk = new double[gridpts];
  
  vector<Point3> bcoefx; bcoefx.resize(gridpts);
  vector<Point3> bcoefy; bcoefy.resize(gridpts);
  int r = bpts;
  
  if(degk+1 <= r && degk%2 == 1 && nrgrid == 3){
    NumMatrix Mx(r, gridpts-ex-1); clear(Mx);
    vector<Point3> tempx;   
    for(int i=0; i<bpts; i++){ //5 points on x-axis
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx == 0) tx++;
      SplineBasis(newx, tx, degk, gridpts,Nj);      
      Mx(i,0) = Nj[1+ex];
      for(int j=0; j<=ex; j++)
	Mx(i,0) += 2*Nj[-j+ex];
      for(int j=2; j<gridpts-ex; j++)
	Mx(i, j-1) = Nj[j+ex] - Nj[2-j+ex];
      tempx.push_back(tmpval[i]);  
    }
    NumMatrix piMx(gridpts-ex-1,r);  pinv(Mx, 1e-16, piMx); // take inverse
    
    vector<Point3> tmp; tmp.resize(gridpts-ex-1);
    matvecmult3(piMx, tempx, tmp);
    
    for(int i=1; i<gridpts-ex; i++)
      bcoefx[i+ex] = tmp[i-1];
    
    //now compute coefficients on y-axis
    r = bpts; //number of rows
    NumMatrix My(r, gridpts-ex-1); clear(My);
    vector<Point3> tempy;
    
    int i=0;
    newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
    ty = (int)ceil(newy);
    if(ty == 0) ty++;
    SplineBasis(newy, ty, degk,gridpts, Nk);
    My(i,0) = Nk[1+ex];
    for(int j=0; j<=ex; j++)
      My(i,0) += 2*Nk[-j+ex];
    for(int j=2; j<gridpts-ex; j++)
      My(i, j-1) = Nk[j+ex] - Nk[2-j+ex];
    tempy.push_back(tmpval[i]);
    
    for(int i=cn-bpts+1; i<cn; i++){
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      ty = (int)ceil(newy);
      if(ty==0) ty++;
      SplineBasis(newy, ty, degk,gridpts, Nk);
      My(i-cn+bpts,0) = Nk[1+ex];
      for(int j=0; j<=ex; j++)
	My(i-cn+bpts,0) += 2*Nk[-j+ex];
      for(int j=2; j<gridpts-ex; j++)
	My(i-cn+bpts, j-1) = Nk[j+ex] - Nk[2-j+ex];
      tempy.push_back(tmpval[i]);
    }
    NumMatrix piMy(gridpts-ex-1,r);  pinv(My, 1e-16, piMy); // take inverse
    matvecmult3(piMy, tempy, tmp);    
    for(int i=1; i<gridpts-ex; i++){
      bcoefy[i+ex] = tmp[i-1];
    }
  
    Point3 mid = Point3(bcoefy[1+ex](0),bcoefx[1+ex](1),bcoefx[1+ex](2));      
    bcoefy[1+ex] = bcoefx[1+ex] = mid;
    for(int i=-ex; i<=0; i++){
      bcoefx[i+ex] = 2.0*mid - bcoefx[-i+2+ex];
      bcoefy[i+ex] = 2.0*mid - bcoefy[-i+2+ex];
    }  
  }
  else{ //not enough pts on boundary (or even)
    if(nrgrid == 2){
      vector<Point3> cubicx; cubicx.resize(4);
      vector<Point3> cubicy; cubicy.resize(4);
      computeSpCubicControlsConcave(bpts, cn, nrgrid, totalx, xmin, totaly, ymin,
				  tmpxy, tmpval, cubicx, cubicy);
      int g=2;
      if(degk%2 ==0) g=1;
      
      vector<Point3> highorderx; highorderx.resize(g+2*ex);
      elevateDegree(degk, cubicx, highorderx);
      vector<Point3> highordery; highordery.resize(g+2*ex);
      elevateDegree(degk, cubicy, highordery);
      
      for(int i=0; i<gridpts; i++){
	bcoefx[i] = highorderx[i];
	bcoefy[i] = highordery[i];
      }
    }
    else if (nrgrid == 3){
      vector<Point3> cubicx; cubicx.resize(5);
      vector<Point3> cubicy; cubicy.resize(5);
      
      computeSpCubicControlsConcave(bpts, cn, nrgrid, totalx, xmin, totaly, ymin,
				  tmpxy, tmpval, cubicx, cubicy);
      int g=2; 
      if(degk%2 == 0) g=1;
      vector<Point3> highx1; highx1.resize(g+2*ex);
      vector<Point3> highx2; highx2.resize(g+2*ex);
      vector<Point3> t1; t1.resize(4);
      vector<Point3> t2; t2.resize(4);
      for(int i=0; i<4; i++)	t1[i] = cubicx[i];
      for(int i=1; i<5; i++)	t2[i-1] = cubicx[i];
      elevateDegree(degk, t1, highx1);
      elevateDegree(degk, t2, highx2);
      
      for(int i=0; i<g+2*ex; i++)	bcoefx[i] = highx1[i];
      for(int i=g+2*ex; i<gridpts; i++)	bcoefx[i] = highx2[i-1];
      
      vector<Point3> highy1; highy1.resize(g+2*ex);
      vector<Point3> highy2; highy2.resize(g+2*ex);
      for(int i=0; i<4; i++)	t1[i] = cubicy[i];
      for(int i=1; i<5; i++)	t2[i-1] = cubicy[i];
      
      elevateDegree(degk, t1, highy1);
      elevateDegree(degk, t2, highy2);
      
      for(int i=0; i<g+2*ex; i++)	bcoefy[i] = highy1[i];
      for(int i=g+2*ex; i<gridpts; i++)	bcoefy[i] = highy2[i-1];
    }    
  }
  
  if(nrgrid == 2){
    for(int i=0; i<gridpts; i++)
      bcoef[i] = bcoefx[i];
    for(int i=1; i<gridpts; i++)
      bcoef[i+gridpts-1] = bcoefy[i];
  }     
  else if (nrgrid==3){
    
    for(int i=1; i<gridpts; i++)
      bcoef[i-1] = bcoefx[i];
    for(int i=0; i<=ex; i++)
      bcoef[i+gridpts-1] = bcoefy[i];     
    for(int i=ex+2; i<gridpts-1; i++)
      bcoef[i+gridpts-2] = bcoefy[i];   
  }  
  delete [] Nj; delete[] Nk;
}
// ---------------------------------------------------------------------- 
void BdSurf::computeSpCubicControlsSmooth(int bpts, int cn, int nrgrid, double totalx, double xmin,
				  vector<Point2>& tmpxy, vector<Point3>& tmpval, vector<Point3>& cubic){
  double newx; int tx;
  double Nj[4];

  NumMatrix M(2*bpts-1, 4); clear(M);
  vector<Point3> temp;  
  for(int i=0; i<bpts; i++){ //first 5 boundary points	    
    newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
    tx = (int)ceil(newx);
    if(tx==0) tx++;
    SplineBasis(newx, tx, 3, 4, Nj);
    for (int j=0; j<4; j++){
      M(i, j) = Nj[j];
    }
    temp.push_back(tmpval[i]);
  }
  
  for (int i=cn-bpts+1; i<cn; i++){ //the other three
    newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
    tx=(int)ceil(newx); 
    if(tx==0) tx++;
    SplineBasis(newx, tx, 3,4, Nj);	    
    for (int j=0; j<4; j++){
      M(i-cn+bpts-1, j) = Nj[j];
    }
    temp.push_back(tmpval[i]);
  }
 
  NumMatrix piM(4,2*bpts-1);	pinv(M, 1e-16, piM); // take inverse
  matvecmult3(piM, temp, cubic); 
}
// ---------------------------------------------------------------------- 
void BdSurf::computeSpCubicControlsConvex(int bpts, int cn, int nrgrid, double totalx, double xmin, double totaly, double ymin,
				  vector<Point2>& tmpxy, vector<Point3>& tmpval, vector<Point3>& cubicx, vector<Point3>& cubicy){
  
  double newx; int tx;
  double newy; int ty;  
  int r = bpts;
  if(nrgrid == 2){    
    double Nj[4];
    double Nk[4];
    
    vector<Point3> rhsx;
    NumMatrix Mx(r, 3); clear(Mx);
    for(int i=0; i<bpts; i++){ //first 4 boundary points	    
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx==0) tx++;
      SplineBasis(newx, tx, 3, 4, Nj);
      Mx(i, 0) = 2*Nj[0]+Nj[1];
      Mx(i, 1) = Nj[2]-Nj[0];
      Mx(i, 2) = Nj[3];
      rhsx.push_back(tmpval[i]);
    }    
    NumMatrix piMx(3,r);	pinv(Mx, 1e-16, piMx); // take inverse
    vector<Point3> tmpx; tmpx.resize(3);
    matvecmult3(piMx, rhsx, tmpx);
    
    //fit y-boun
    NumMatrix My(r, 3); clear(My);
    vector<Point3> rhsy;
    int i=0;
    newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
    ty = (int)ceil(newy);
    if(ty==0) ty++;
    SplineBasis(newy, ty, 3, 4, Nk);
    My(i, 0) = 2*Nk[0]+Nk[1];
    My(i, 1) = Nk[2]-Nk[0];
    My(i, 2) = Nk[3];
    rhsy.push_back(tmpval[i]);
    
    for(int i=cn-bpts+1; i<cn; i++){
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      ty = (int)ceil(newy);
      if(ty==0) ty++;
      SplineBasis(newy, tx, 3, 4, Nk);
      My(i-cn+bpts, 0) = 2*Nk[0]+Nk[1];
      My(i-cn+bpts, 1) = Nk[2]-Nk[0];
      My(i-cn+bpts, 2) = Nk[3];
      rhsy.push_back(tmpval[i]);
    }    
    NumMatrix piMy(3,r);	pinv(My, 1e-16, piMy); // take inverse
    vector<Point3> tmpy; tmpy.resize(3);
    matvecmult3(piMy, rhsy, tmpy);
    
    Point3 mid = Point3(tmpy[0](0),tmpx[0](1),tmpx[0](2));
    cubicx[0] = 2.0*mid - tmpx[1];
    cubicy[0] = 2.0*mid - tmpy[1];
    for(int i=1; i<4; i++){
      cubicx[i] = tmpx[i-1];
      cubicy[i] = tmpy[i-1];
    }    
  }
  else if (nrgrid==3){

    double Nj[5];
    double Nk[5];
    
    NumMatrix Mx(r, 4); clear(Mx);
    vector<Point3> rhsx;
    for(int i=0; i<bpts; i++){ //first 4 boundary points	    
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx==0) tx++;
      SplineBasis(newx, tx, 3, 5, Nj);
      Mx(i, 0) = 2*Nj[0] + Nj[1]; 
      Mx(i, 1) = Nj[2] - Nj[0];
      Mx(i, 2) = Nj[3];
      Mx(i, 3) = Nj[4];
      rhsx.push_back(tmpval[i]);
    }
    NumMatrix piMx(4,r);	pinv(Mx, 1e-16, piMx); // take inverse
    vector<Point3> tmpx; tmpx.resize(4);
    matvecmult3(piMx, rhsx, tmpx);
    
    //fit y-boun
    NumMatrix My(r, 4); clear(My);
    vector<Point3> rhsy;

    int i=0;
    newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
    ty = (int)ceil(newy);
    if(ty==0) ty++;
    SplineBasis(newy, ty, 3, 5, Nk);
    My(i, 0) =  2*Nk[0] + Nk[1]; 
    My(i, 1) = Nk[2] - Nk[0];
    My(i, 2) = Nk[3];
    My(i, 3) = Nk[4];
    rhsy.push_back(tmpval[i]);    
    for(int i=cn-bpts+1; i<cn; i++){
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      ty = (int)ceil(newy);
      if(ty==0) ty++;
      SplineBasis(newy, tx, 3, 5, Nk);
      My(i-cn+bpts, 0) = 2*Nk[0]+Nk[1];
      My(i-cn+bpts, 1) = Nk[2]-Nk[0];
      My(i-cn+bpts, 2) = Nk[3];
      My(i-cn+bpts, 3) = Nk[4];
      rhsy.push_back(tmpval[i]);
    }
    NumMatrix piMy(4,r);	pinv(My, 1e-16, piMy); // take inverse
    vector<Point3> tmpy; tmpy.resize(4);
    matvecmult3(piMy, rhsy, tmpy);
 
    Point3 mid = Point3(tmpy[0](0),tmpx[0](1),tmpx[0](2));
    cubicx[0] = 2.0*mid - tmpx[1];
    cubicy[0] = 2.0*mid - tmpy[1];
    cubicx[1] = cubicy[1] = mid;
    
    for(int i=1; i<3; i++){
      cubicx[1+i] = tmpx[i];
      cubicy[1+i] = tmpy[i];
    }
  }
}

// ---------------------------------------------------------------------- 
void BdSurf::computeSpCubicControlsConcave(int bpts, int cn, int nrgrid, double totalx, double xmin, 
					 double totaly, double ymin, vector<Point2>& tmpxy, vector<Point3>& tmpval, 
					 vector<Point3>& cubicx, vector<Point3>& cubicy){
  
  double newx; int tx;
  double newy; int ty;
  int r = bpts;
  
  if(nrgrid == 2){
    
    double Nj[4];
    double Nk[4];
    double N0[4];
   
    newx = (tmpxy[0](0) -xmin)*(nrgrid-1)/totalx;
    tx = (int)ceil(newx);
    if(tx==0) tx++;
    SplineBasis(newx, tx, 3, 4, N0);

    NumMatrix Mx(bpts-1, 3); clear(Mx);
    vector<Point3> rhsx;
    for(int i=1; i<bpts; i++){
      newx = (tmpxy[i](0) -xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx==0) tx++;
      SplineBasis(newx, tx, 3, 4, Nj);
      for(int j=0; j<3; j++)
	Mx(i-1, j) = Nj[j+1] - ((N0[j+1]/N0[0])*Nj[0]);
      rhsx.push_back(tmpval[i] - (tmpval[0]*(Nj[0]/N0[0])));
    }
    
    NumMatrix piMx(3,bpts-1);	pinv(Mx, 1e-16, piMx); // take inverse
    vector<Point3> tmpx; tmpx.resize(3);
    matvecmult3(piMx, rhsx, tmpx);
  
    NumMatrix My(bpts-1, 3); clear(My);
    vector<Point3> rhsy;
    for(int i=cn-bpts+1; i<cn; i++){
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      ty = (int)ceil(newy);
      if(ty==0) ty++;
      SplineBasis(newy, ty, 3, 4, Nk);
      for(int j=0; j<3; j++)
	My(i-cn+bpts-1,j) = Nk[j+1] - ((1.0/N0[0])*Nk[0]*N0[j+1]);
      rhsy.push_back(tmpval[i]- (tmpval[0]*(1.0/N0[0])*Nk[0]));
    }

    NumMatrix piMy(3,bpts-1);	pinv(My, 1e-16, piMy); // take inverse
    vector<Point3> tmpy; tmpy.resize(3);
    matvecmult3(piMy, rhsy, tmpy);

    cubicx[0] = (1.0/N0[0])*(tmpval[0] -tmpx[0]*N0[1] - tmpx[1]*N0[2] - tmpx[2]*N0[3]);
    cubicy[0] = (1.0/N0[0])*(tmpval[0] -tmpy[0]*N0[1] - tmpy[1]*N0[2] - tmpy[2]*N0[3]);
    for(int i=1; i<4; i++){
      cubicx[i] = tmpx[i-1];
      cubicy[i] = tmpy[i-1];
    }
  }
  
  else if(nrgrid == 3){
    
    double Nj[5];
    double Nk[5];
    //fit x-boun
    
    vector<Point3> rhsx;
    NumMatrix Mx(r, 3); clear(Mx);
    for(int i=0; i<bpts; i++){ //first 4 boundary points	    
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      tx = (int)ceil(newx);
      if(tx==0) tx++;
      SplineBasis(newx, tx, 3, 5, Nj);
      Mx(i, 0) = 2*Nj[0]+2*Nj[1]+Nj[2];
      Mx(i, 1) = Nj[3]-Nj[1];
      Mx(i, 2) = Nj[4]-Nj[0];
      rhsx.push_back(tmpval[i]);
    }
    
    NumMatrix piMx(3,r);	pinv(Mx, 1e-16, piMx); // take inverse
    vector<Point3> tmpx; tmpx.resize(3);
    matvecmult3(piMx, rhsx, tmpx);
    
    NumMatrix My(r, 3); clear(My);
    vector<Point3> rhsy;
    
    int i=0;
    newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
    ty = (int)ceil(newy);
    if(ty==0) ty++;
    SplineBasis(newy, ty, 3, 5, Nk);
    My(i, 0) = 2*Nk[0]+2*Nk[1]+Nk[2];
    My(i, 1) = Nk[3]-Nk[1];
    My(i, 2) = Nk[4]-Nk[0];
    rhsy.push_back(tmpval[i]);
    for(int i=cn-bpts+1; i<cn; i++){
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      ty = (int)ceil(newy);
      if(ty==0) ty++;
      SplineBasis(newy, tx, 3, 5, Nk);
      My(i-cn+bpts, 0) = 2*Nk[0]+2*Nk[1]+Nk[2];
      My(i-cn+bpts, 1) = Nk[3]-Nk[1];
      My(i-cn+bpts, 2) = Nk[4]-Nk[0];
      rhsy.push_back(tmpval[i]);
    }
    
    NumMatrix piMy(3,r);	pinv(My, 1e-16, piMy); // take inverse
    vector<Point3> tmpy; tmpy.resize(3);
    matvecmult3(piMy, rhsy, tmpy);
    
    Point3 mid = Point3(tmpy[0](0),tmpx[0](1),tmpx[0](2));
    cubicx[0] = 2.0*mid - tmpx[2];
    cubicy[0] = 2.0*mid - tmpy[2];
    cubicx[1] = 2.0*mid - tmpx[1];
    cubicy[1] = 2.0*mid - tmpy[1];
    cubicx[2] = cubicy[2] = mid;
    
    for(int i=1; i<3; i++){
      cubicx[2+i] = tmpx[i];
      cubicy[2+i] = tmpy[i];
    }
  }
}

// ---------------------------------------------------------------------- 
void BdSurf::createSpIndepMatrices(int flag, int degk, int nrgrid, vector<Point2> &tmpxy,
				 double xmin, double ymin, double totalx, double totaly,
				 NumMatrix &A1, NumMatrix &A2, NumMatrix &B1, NumMatrix &B2){
  
  int ex = (int)floor((degk)/2.0);  
  double * Nk = new double[nrgrid+2*ex];
  double * Nj = new double[nrgrid+2*ex];
  int cn = A1.m();
  int nre = B2.m();
  int cm = B1.n() + B2.n();
  double newx, newy;
  int tx, ty;
  int gridpts = nrgrid+2*ex;
  if(degk %2 == 0) gridpts--;
  
  if(flag == SMOOTH_MAT){
    
    SplineBasis(0, 1, degk, nrgrid+2*ex, Nk);
    
    for(int i=0; i<nre; i++)
      B2(i, i) = Nk[0+ex];
    for(int i=0; i<nre; i++){
      for(int j=-ex; j<0; j++)
	B1(i, i*(nre-1)+j+ex) =  Nk[j+ex];
      for(int j=1; j<nre-ex; j++)
	B1(i, i*(nre-1)+j+ex-1) =Nk[j+ex];
    }
    
    for(int i=0; i<cn; i++){
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      tx = (int)ceil(newx); ty = (int)ceil(newy);
      if(tx==0) tx++; 
      if(ty==0) ty++;
      SplineBasis(newx, tx, degk, nrgrid+2*ex, Nj);
      SplineBasis(newy, ty, degk, nrgrid+2*ex, Nk);
      
      for(int j=-ex; j<nrgrid+ex; j++){
	for(int k=-ex; k<0; k++)
	  A1(i, (j+ex)*(nre-1)+(k+ex)) = Nj[j+ex]*Nk[k+ex];
	for(int k=1; k<nrgrid+ex; k++)
	  A1(i, (j+ex)*(nre-1)+(k-1+ex)) = Nj[j+ex]*Nk[k+ex];
	
	A2(i, j+ex) = Nj[j+ex]*Nk[0+ex];//k=0
      }
    }
  }
  
  else if(flag == CONVEX_MAT){
    
    SplineBasis(0, 1, degk, gridpts, Nk);
    
    //create B2
    for(int i=0; i<nre; i++)
      B2(i, i) = Nk[0+ex];
    for(int i=-ex; i<0; i++)
      B2(ex, gridpts+ex+i) = Nk[i+ex];
    for(int i=1; i<gridpts-ex; i++)
      B2(ex, gridpts+ex+i-1)= Nk[i+ex];
    
    //create B1
    for(int i=0; i<ex; i++){
      for(int j=-ex; j<0; j++)
	B1(i, i*(gridpts-1)+j+ex) =  Nk[j+ex];  
      for(int j=1; j<gridpts-ex; j++)
	B1(i, i*(gridpts-1)+j+ex-1) = Nk[j+ex];
    }
    for(int i=ex+1; i<gridpts; i++){
      for(int j=-ex; j<0; j++)
	B1(i, (i-1)*(gridpts-1)+j+ex) =  Nk[j+ex];
      for(int j=1; j<gridpts-ex; j++)
	B1(i, (i-1)*(gridpts-1)+j+ex-1) = Nk[j+ex];
    }
    for(int i=gridpts; i<nre; i++){
      for(int j=-ex; j<0; j++)   
	B1(i,   (j+ex-1)*(gridpts-1)+i-1) = Nk[j+ex];
      for(int j=1; j<gridpts-ex; j++)
	B1(i, (j+ex-2)*(gridpts-1)+i-1) = Nk[j+ex];
    }
    
    for(int i=0; i<cn; i++){
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      tx = (int)ceil(newx); ty = (int)ceil(newy);
      if(tx==0) tx++; 
      if(ty==0) ty++;
      SplineBasis(newx, tx, degk, gridpts, Nj);
      SplineBasis(newy, ty, degk, gridpts, Nk);
      
      for(int j=-ex; j<0; j++){
	for(int k=-ex; k<0; k++)
	  A1(i, (j+ex)*(gridpts-1)+(k+ex)) = Nj[j+ex]*Nk[k+ex];
	for(int k=1; k<gridpts-ex; k++)
	  A1(i, (j+ex)*(gridpts-1)+(k-1+ex)) = Nj[j+ex]*Nk[k+ex];
      }
      
      for(int j=1; j<gridpts-ex; j++){
	for(int k=-ex; k<0; k++)
	  A1(i, (j-1+ex)*(gridpts-1)+(k+ex)) = Nj[j+ex]*Nk[k+ex];
	for(int k=1; k<gridpts-ex; k++)
	  A1(i, (j-1+ex)*(gridpts-1)+(k-1+ex)) = Nj[j+ex]*Nk[k+ex];
      }
      
      for(int j=-ex; j<gridpts-ex; j++)
	A2(i, j+ex) = Nj[j+ex]*Nk[0+ex];
      for(int k=-ex; k<0; k++)
	A2(i, gridpts+ex+k) = Nj[0+ex]*Nk[k+ex];
      for(int k=1; k<gridpts-ex; k++)
	A2(i, gridpts+ex+k-1) = Nj[0+ex]*Nk[k+ex];
    }
  }

  else if(flag == CONCAVE2_MAT){
 
    SplineBasis(.5, 1, degk, gridpts, Nk);
    
    //create B2
    for(int i=0; i<gridpts; i++)
      B2(i, i) = Nk[-ex+ex]; 
    for(int i=gridpts; i<nre; i++)
      B2(i, i) = Nk[gridpts-ex-1+ex];
    for(int i=-ex+1; i<gridpts-ex; i++)
      B2(gridpts-1, gridpts+ex+i-1) = Nk[i+ex];
    
    //create B1
    for(int i=0; i<gridpts-1; i++){
      for(int j=-ex+1; j<gridpts-ex; j++)
	B1(i, (i)*(gridpts-1)+j+ex-1) =  Nk[j+ex];
    }
    for(int i=gridpts; i<nre; i++){
      for(int j=-ex; j<gridpts-ex-1; j++)
	B1(i, (j+ex-1)*(gridpts-1)+i-1) = Nk[j+ex]; 
    }
    
    //create A1 & A2
    for(int i=0; i<cn; i++){  
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      tx = (int)ceil(newx); ty = (int)ceil(newy);
      if(tx==0) tx++; 
      if(ty==0) ty++;
      SplineBasis(newx, tx, degk, gridpts, Nj);
      SplineBasis(newy, ty, degk, gridpts, Nk); 
      for(int j=-ex; j<gridpts-ex-1; j++){
	for(int k=-ex+1; k<gridpts-ex; k++)
	  A1(i, (j+ex)*(gridpts-1)+(k+ex-1)) = Nj[j+ex]*Nk[k+ex];
      }
      
      for(int j=-ex; j<gridpts-ex; j++)
	A2(i, j+ex) = Nj[j+ex]*Nk[-ex+ex];
      for(int k=-ex+1; k<gridpts-ex; k++)
	A2(i, gridpts+ex+k-1) = Nj[gridpts-ex-1+ex]*Nk[k+ex];
    }
  }

  else if(flag == CONCAVE3_MAT){

    SplineBasis(1, 1, degk, gridpts, Nk);

    //create B2
    for(int i=0; i<nre; i++)
      B2(i, i) = Nk[1+ex];
    for(int i=-ex; i<=0; i++)
      B2(ex, gridpts-1+ex+i) = Nk[i+ex];
    for(int i=2; i<gridpts-ex-1; i++)
      B2(ex, gridpts-1+ex+i-1)= Nk[i+ex];
    
    //create B1
    for(int i=0; i<ex; i++){
      for(int j=-ex; j<=0; j++)
	B1(i, (i+1)*(gridpts-1)+j+ex) =  Nk[j+ex];
      for(int j=2; j<gridpts-ex; j++)
	B1(i, (i+1)*(gridpts-1)+j+ex-1) = Nk[j+ex];   
    }
    B1(ex, cm-nre-1) = Nk[gridpts-ex+ex-1];
    for(int i=ex+1; i<gridpts-1; i++){
      for(int j=-ex; j<=0; j++)
	B1(i, (i)*(gridpts-1)+j+ex) =  Nk[j+ex];    
      for(int j=2; j<gridpts-ex; j++)
	B1(i, (i)*(gridpts-1)+j+ex-1) = Nk[j+ex];
    }
    for(int i=gridpts-1; i<nre; i++){
      for(int j=-ex; j<=0; j++) 
	B1(i, (j+ex-1)*(gridpts-1)+i) = Nk[j+ex];
      for(int j=2; j<gridpts-ex; j++)
	B1(i, (j+ex-2)*(gridpts-1)+i) = Nk[j+ex];
    }

    //create A1 & A2

    for(int i=0; i<cn; i++){        
      newx = (tmpxy[i](0)-xmin)*(nrgrid-1)/totalx;
      newy = (tmpxy[i](1)-ymin)*(nrgrid-1)/totaly;
      tx = (int)ceil(newx); ty = (int)ceil(newy);
      if(tx==0) tx++; 
      if(ty==0) ty++;
      SplineBasis(newx, tx, degk, gridpts, Nj);
      SplineBasis(newy, ty, degk, gridpts, Nk);
      
      for(int j=-ex; j<=0; j++){
	for(int k=-ex; k<=0; k++)
	  A1(i, (j+ex)*(gridpts-1)+(k+ex)) = Nj[j+ex]*Nk[k+ex];
	for(int k=2; k<gridpts-ex; k++)
	  A1(i, (j+ex)*(gridpts-1)+(k-1+ex)) = Nj[j+ex]*Nk[k+ex];
      }
      for(int j=2; j<gridpts-ex; j++){
	for(int k=-ex; k<=0; k++)
	  A1(i, (j-1+ex)*(gridpts-1)+(k+ex)) = Nj[j+ex]*Nk[k+ex];
	for(int k=2; k<gridpts-ex; k++)
	  A1(i, (j-1+ex)*(gridpts-1)+(k-1+ex)) = Nj[j+ex]*Nk[k+ex];
      }
      A1(i, cm-nre-2) =  Nj[-ex+ex]*Nk[1+ex];
      A1(i, cm-nre-1) =  Nj[1 +ex ]*Nk[gridpts-ex-1+ex];
      
      for(int j=-ex+1; j<gridpts-ex; j++)
	A2(i, j+ex-1) = Nj[j+ex]*Nk[1+ex];
      for(int k=-ex; k<=0; k++)
	A2(i, gridpts+ex+k-1) = Nj[1+ex]*Nk[k+ex];
      for(int k=2; k<gridpts-ex-1; k++)
	A2(i, gridpts+ex+k-2) = Nj[1+ex]*Nk[k+ex];
    }
  }
  else 
    assert(0);

  delete[] Nj; delete[]Nk;

}


// ---------------------------------------------------------------------- 
void BdSurf::computeSpIndepControls(NumMatrix& B1, NumMatrix& B2, NumMatrix& A1, NumMatrix& A2, 
				  vector<Point3> &bcoef,vector<Point3> &tmpval, 
				  vector<Point3>& p1, vector<Point3>& p2){
  
  int cn = A1.m();
  int nre = B2.m();
  int cm = B1.n() + B2.n();

  NumMatrix pB2(nre,nre);	pinv(B2, 1e-16, pB2); // take inverse
  NumMatrix T(cn, nre); 
  NumMatrixMult(A2, pB2, T);
  NumMatrix T2(cn, cm-nre); 
  NumMatrixMult(T, B1, T2);
  multvalue<double>(T2, -1);
  NumMatrix N(cn, cm-nre);
  NumMatrixAdd(A1, T2, N);
  NumMatrix pN(cm-nre, cn); pinv(N, 1e-16, pN); // take inverse
  
  vector<Point3> Tx; Tx.resize(cn);
  matvecmult3(T, bcoef, Tx);

  vector<Point3> rhs;
  for(int i=0; i<cn; i++)
    rhs.push_back(tmpval[i] - Tx[i]);
  
  matvecmult3(pN, rhs, p1);
  vector<Point3> temp;  temp.resize(nre);
  matvecmult3(B1, p1, temp);
  
  for(int i=0; i<nre; i++)
    temp[i] = bcoef[i]-temp[i];
  matvecmult3(pB2, temp, p2);
  
}


// ---------------------------------------------------------------------- 
void BdSurf::elevateDegree(int deg, vector<Point3>& cubic, vector<Point3> & high){
  int n = 3;       //starting degree
  int r = deg-n; //nr of elevations
   
  if (r==0){
    for(int i=0; i<4; i++)
      high[i] = cubic[i];
  }
  else{
    switchToBezier(cubic);
    vector<Point3> temp; temp.resize(deg+1);
    for(int i=0; i<4; i++)   temp[i] = cubic[i];
    for(int i=0; i<r; i++){
      high[0] = temp[0];
      high[n+1] = temp[n];
      for(int j=1; j<=n;j++){
	double k = (double)j/((double)(n+1));
	high[j] = k*temp[j-1] + (1-k)*temp[j];
      }
      n++;
      for(int j =0; j<n+1; j++) temp[j] = high[j];
    }
    switchToBSpline(high);
  }
}
// ---------------------------------------------------------------------- 
void BdSurf::switchToBezier(vector<Point3>&cubic){    
  NumMatrix M(4,4);
  M(0,0) = 1.0; M(0,1) = 4.0; M(0,2) = 1.0; M(0,3) = 0.0;
  M(1,0) = 0.0; M(1,1) = 4.0; M(1,2) = 2.0; M(1,3) = 0.0;
  M(2,0) = 0.0; M(2,1) = 2.0; M(2,2) = 4.0; M(2,3) = 0.0;
  M(3,0) = 0.0; M(3,1) = 1.0; M(3,2) = 4.0; M(3,3) = 1.0;
  
  multvalue<double>(M,(1.0/6.0));
  vector<Point3> temp; temp.resize(4);
  matvecmult3(M, cubic,temp);
  for(int i=0; i<4; i++)    cubic[i] = temp[i];
}

// ---------------------------------------------------------------------- 
void BdSurf::switchToBSpline(vector<Point3>&high){
  
  int n = high.size()-1;
  int kk, nc, fc;
  //create conversion matrix
  NumMatrix R(n+1, n+1); clear(R);
  for(int i=0; i<n+1; i++) R(i, i) = 1.0;
  for(int k = n-1; k>=1; k--){
    kk = n+1-k;
    nc = (int)round((n-k)/2.0);
    fc =n+1-nc;
    for(int j=fc; j<=n; j++){
      R(k-1, j-1) = kk*(R(k, j-1) - R(k, j));
      for(int i=k+2; i<=n+1; i++)
	R(i-1, j-1) = kk*(R(i-1, j-1) - R(i-1, j));
    }
    R(n, n) = kk*R(n,n);
    for(int j=n; j>=fc; j--){      
      R(k-1, j-1) = R(k, j-1) - R(k-1, j-1);
      for(int i=k+2; i<=n+1; i++)
	R(i-1, j-1) = R(i-2, j-1) + R(i-1, j-1);
    }
    if(kk%2 == 1)
      for(int i=k-1; i<=n; i++)
	R(i, n-nc-1) = R(n-(i-k+1), fc-1);
  }
  for(int j=1; j<=(int)round(n/2.0); j++){
    for(int i=1; i<=n+1;i++)
      R(i-1, j-1) = R(n-i+1, n+1-j);
  }

  //matrix done, now multiply
  vector<Point3> temp; temp.resize(n+1);
  matvecmult3(R, high, temp);
  for(int i=0; i<n+1; i++)    high[i] = temp[i];
  
}
// ---------------------------------------------------------------------- 

//  
//----------------------------------------------------------------------
int BdSurf::Vfcd2Vxy_val(int flags, int K,int bountype, int contype,  
int f, double* cd, double* xy){



     NumVector D;
     double p = 0; double q = 0;  //constants
     if (bountype== GpMesh::INTERIOR_VERTEX) {//not boundary


       D = _subdivMatrices->chooseMatrix(CCSubMatLib::INTERIOR, K).D();

if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 4.0/K;	 q = 4.0/K;
       } else if(_chartType == CHARACTERISTIC_MAP_CHART) { // CHARACTERISTIC MAP
	p = log(1.0/D(1))/log(2.0);	 q = 4.0/K;
       } else if(_chartType == ISODISTANCE_CHART) {
	p = 1.0;	 q = 4.0/K;
       }
     }

     else if(contype == GpMesh::CREASE_VERTEX){
       D = _subdivMatrices->chooseMatrix(CCSubMatLib::BOUNDARY, K).D();
       if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 2.0/K;	 q = 2.0/K;
       } else if(_chartType == CHARACTERISTIC_MAP_CHART) { // CHARACTERISTIC MAP
	p = log(1.0/D(1))/log(2.0);	 q = 2.0/K;
       } else if(_chartType == ISODISTANCE_CHART) {
	p = 1.0;	 q = 2.0/K;
       }
     }
     else if(contype == GpMesh::CONVEX_VERTEX){
       D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONVEX,  
K).D();
       if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 1.0/K;	 q = 1.0/K;
       } else if(_chartType == CHARACTERISTIC_MAP_CHART) { // CHARACTERISTIC MAP
	p = log(1.0/D(1))/log(2.0);	 q = 1.0/K;
       } else if(_chartType == ISODISTANCE_CHART) {
	p = 1.0;	 q = 1.0/K;
       }
     }
     else if(contype == GpMesh::CONCAVE_VERTEX){
       D = _subdivMatrices->chooseMatrix(CCSubMatLib::CORNERCONCAVE,  K).D();
       if(_chartType == FULLY_COMPLEX_CHART) { //FULLY COMPLEX
	p = 3.0/K;	 q = 3.0/K;
       } else {//CHARACTERISTIC = ISODISTANCE
	p = 1.0;	 q = 3.0/K;
       }
     }


     //---------------------
     if( abs(cd[0])<=SCL_EPS && abs(cd[1])<=SCL_EPS ) {
       //-------------------	 //CANNOT TRANSFROM AT C==0 && D==0	
       double* val = xy;	
       if(flags & EVAL_VALUE) {
	val[0] = 0; val[1] = 0;
       }
       if(flags & EVAL_1ST_DERIV) {assert(0);}
       if(flags & EVAL_2ND_DERIV) {assert(0);}
     }
     else {

       double c = cd[0]; double d = cd[1];
       //-------------------
       double ra2cd[12];
       double r2 = c*c+d*d;  double r = sqrt(r2);
       double a = atan2(d,c);
       double r3 = r2*r;  double r4 = r2*r2;

       if(flags & EVAL_VALUE) {
	ra2cd[0] = r;  ra2cd[1] = a;
       }
       if(flags & EVAL_1ST_DERIV) {
	ra2cd[2] = c/r;  ra2cd[3] = -d/r2;
	ra2cd[4] = d/r;  ra2cd[5] = c/r2;
       }
       if(flags & EVAL_2ND_DERIV) {
	ra2cd[6] = 1/r-c*c/r3;  ra2cd[7] = 2*d*c/r4;
	ra2cd[8] = -c*d/r3;     ra2cd[9] = (-c*c+d*d)/r4;
	ra2cd[10] = 1/r-d*d/r3; ra2cd[11] = -2*d*c/r4;
       }
       //-------------------
       double st2ra[12];
       double rpm2= exp((p-2)*log(r)); // pow(r,p-2);
       double  rpm1 = rpm2*r;//pow(r,p-1);
       double  rp = rpm1*r; //pow(r,p);
       double aq = a*q;
       double caq, saq;
       saq = sin(aq); caq= cos(aq);
       //sincos(aq, &saq, &caq);

       if(flags & EVAL_VALUE)
	st2ra[0] = rp*caq;  st2ra[1] = rp*saq;

       if(flags & EVAL_1ST_DERIV) {
	st2ra[2] = rpm1*p*caq;  st2ra[3] = rpm1*p*saq;
	st2ra[4] = -q*st2ra[1];   st2ra[5] = st2ra[0]*q;
       }
       if(flags & EVAL_2ND_DERIV) {
	st2ra[6] = rpm2*p*caq*(p-1);  st2ra[7]  = rpm2*p*saq*(p-1);
	st2ra[8] = -rpm1*p*saq*q;     st2ra[9]  = rpm1*p*caq*q;
	st2ra[10] = -q*rp*st2ra[5];   st2ra[11] = st2ra[4]*q;
       }

       //-------------------
       double xy2st[12];
       double R = 0;
       double foverk = (double(f)/double(K));
       if (bountype == GpMesh::INTERIOR_VERTEX) //not boundary
	R = 2*PI*foverk;
       else if(contype == GpMesh::CREASE_VERTEX)
	R = PI*foverk;
       else if(contype == GpMesh::CONVEX_VERTEX)
	R = (PI/2.0)*foverk;
       else if(contype == GpMesh::CONCAVE_VERTEX)
	R = (3*PI/2.0)*foverk;

       double cR, sR;
       sR = sin(R); cR = cos(R);
       //sincos(R, &sR, &cR);
       double s = st2ra[0];  double t = st2ra[1];
       if(flags & EVAL_VALUE) {
	xy2st[0] = cR*s - sR*t;  xy2st[1] = sR*s + cR*t;
       }
       if(flags & EVAL_1ST_DERIV) {
	
	xy2st[2] = cR;  xy2st[3] = sR;
	xy2st[4] = -sR; xy2st[5] = cR;
       }
       if(flags & EVAL_2ND_DERIV) {
	double* tmp = xy2st+6;
	for(int g=0; g<6; g++) tmp[g]=0;
       }
       double st2cd[12];
       compose(flags, 2, st2ra, ra2cd, st2cd);
       compose(flags, 2, xy2st, st2cd, xy);
     }
     return 0;
}
