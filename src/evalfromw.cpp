#include <iostream> 
#include <fstream>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <string>


using namespace std; 

string MatricesDir;

const int MAX_VALENCE = 13;
const int START_VALENCE = 3;
int NUM_PATCHES = 3; 
const int MAX_D = 5;

enum ChartType { CHART_INTERIOR = 0,  CHART_BOUN = 1, CHART_CONVEX = 2,
		 CHART_CONCAVE = 3, MAX_CHART_TYPE=3};

const char* ChartName[] = {"interior","boun","convex","concave"};
const int  ExtrPts[] = { 1,5,5,5};

double* W[MAX_CHART_TYPE+1][MAX_D+1][MAX_VALENCE+1];
double* SamplePts[MAX_CHART_TYPE+1][MAX_VALENCE+1];
double* Samples[MAX_CHART_TYPE+1][MAX_VALENCE+1];
double* CPs[MAX_CHART_TYPE+1][MAX_D+1][MAX_VALENCE+1];
double MinX[MAX_CHART_TYPE+1][MAX_VALENCE+1];
double MaxX[MAX_CHART_TYPE+1][MAX_VALENCE+1];
double MinY[MAX_CHART_TYPE+1][MAX_VALENCE+1];
double MaxY[MAX_CHART_TYPE+1][MAX_VALENCE+1];

void loadMatrix(istream& in, int chart_type) { 
  int d, valence, rowdim, coldim;
  while(in) { 
    //    in >> chart_type;
    in >> valence; if(!in) break;
    in >> rowdim;  if(!in) break;
    in >> coldim;  if(!in) break;
    in >> d;  if(!in) break;

    //    assert(chart_type <= MAX_CHART_TYPE);
    assert(d <= MAX_D);
    assert(d <= MAX_VALENCE);    
    cerr << valence << " " << rowdim <<  " " << d << " " << NUM_PATCHES << endl;
    assert(rowdim == (d + NUM_PATCHES)*(d + NUM_PATCHES));
    assert(coldim == 20*valence + ExtrPts[chart_type]);

    W[chart_type][d][valence] = new double[rowdim*coldim];
    int cnt = 0; 

    while( in && cnt < rowdim*coldim) { 
      in >> W[chart_type][d][valence][cnt];
      cnt++; 
    }
    if (cnt < rowdim*coldim) {
      cerr << "error reading matrix for d " << d <<  " valence " << valence << 
	" cnt " << cnt << endl;       
    }
  }
}

void loadSamplePoints(istream& in, int chart_type) { 
  int valence;
  while(in) { 
    in >> valence; if(!in) break;        
    cerr << "loading valence " << valence << endl;
    assert(chart_type <= MAX_CHART_TYPE);
    int num_samples = 20*valence + ExtrPts[chart_type];
    SamplePts[chart_type][valence] = new double[2*num_samples];
    int cnt = 0; 
    double minx = DBL_MAX, maxx = -DBL_MAX ,miny = DBL_MAX,maxy = -DBL_MAX;
    while( in && cnt < num_samples) { 
      in >> SamplePts[chart_type][valence][2*cnt];
      in >> SamplePts[chart_type][valence][2*cnt+1];

      minx = min(minx, SamplePts[chart_type][valence][2*cnt]);
      maxx = max(maxx, SamplePts[chart_type][valence][2*cnt]);
      miny = min(miny, SamplePts[chart_type][valence][2*cnt+1]);
      maxy = max(maxy, SamplePts[chart_type][valence][2*cnt+1]);      
      cnt++; 
    }

    if (cnt < num_samples) {
      cerr << "error reading sample points  for valence" << "val" << endl;       
    }
    MinX[chart_type][valence] = minx;
    MaxX[chart_type][valence] = maxx;
    MinY[chart_type][valence] = miny;
    MaxY[chart_type][valence] = maxy;

  }
}



void loadAllSamplePoints() { 
  ifstream in_S;
  cerr << "loading " << endl;
  for (int i=0; i <= MAX_CHART_TYPE; i++) {
    in_S.open( (string("dat/sample_pts_")+string(ChartName[i])+string(".dat")).c_str());     
    cerr << "loading from " << 
      (string("dat/sample_pts_")+string(ChartName[i])+string(".dat")).c_str() << endl; 
    loadSamplePoints(in_S,i);
    in_S.close();
  }
}


void loadAllMatrices() { 
  ifstream in_W;
  loadAllSamplePoints();
  cerr << "loading " << endl;
  for (int i=0; i <= MAX_CHART_TYPE; i++) {
    in_W.open( (MatricesDir+ string("W_")+string(ChartName[i])+string(".dat")).c_str());     
    cerr << "loading from " << 
      (MatricesDir + string("W_")+string(ChartName[i])+string(".dat")).c_str() << endl; 
    loadMatrix(in_W,i);
    in_W.close();
  }
}

double testfun(double u, double v) { 
  //    return   0.2*((u-(NUM_PATCHES-1)/2.)*(u-(NUM_PATCHES-1)/2.) + (v-(NUM_PATCHES-1)/2.)*(v-(NUM_PATCHES-1)/2.)); 
  //  return (u-1)*(u-1)+(v-1)*(v-1);
  //  return sin(u*M_PI) + sin(v*M_PI);
  //  return (u-0.8)*(u-0.8)+(v-0.8)*(v-0.8);
    return 
       0.3*(u-v)*(u-v)+
       pow(fabs(v+u-double(NUM_PATCHES-1)),1.5);
}


void  createSamples() {
  for(int chart_type = 0; chart_type <= MAX_CHART_TYPE; chart_type++) 
    for(int valence = 3; valence <= MAX_VALENCE; valence++) {
      //      if( SamplePts[chart_type] != 0) {
	Samples[chart_type][valence] = new double[20*valence + ExtrPts[chart_type]];
	for(int i = 0; i < 20*valence + ExtrPts[chart_type]; i++) {
	  double u =   (NUM_PATCHES-1)*(SamplePts[chart_type][valence][2*i]  -MinX[chart_type][valence])/(MaxX[chart_type][valence]-MinX[chart_type][valence]);
	  double v =   (NUM_PATCHES-1)*(SamplePts[chart_type][valence][2*i+1]-MinY[chart_type][valence])/(MaxY[chart_type][valence]-MinY[chart_type][valence]);
	  if( chart_type == CHART_CONCAVE) v =  NUM_PATCHES-1-v;
	  Samples[chart_type][valence][i] = 
	    testfun(u,v);

	}
	//    }
    }
}

void createCPs() { 
  
  for(int chart_type = 0; chart_type <= MAX_CHART_TYPE; chart_type++)
    for(int d = 0; d <= MAX_D; d++)
      for(int valence = 3; valence <= MAX_VALENCE; valence++) {
	int ncp = d+NUM_PATCHES;
	CPs[chart_type][d][valence] = new double[ncp*ncp];
	for(int l = 0; l < ncp*ncp; l++) {
	  CPs[chart_type][d][valence][l] = 0.0;
	  for(int i = 0; i < 20*valence + ExtrPts[chart_type]; i++) {
	  	    CPs[chart_type][d][valence][l] += 
		      W[chart_type][d][valence][ (20*valence + ExtrPts[chart_type])*l+i]*
		      Samples[chart_type][valence][i];
	  }
	}
      }
}



double N00(double u) { 
  if( u <0.0) return 0.0; 
  else if (u < 1.0) return 1.0; 
  else return 0.0;
}

extern "C" { 
extern double N0(double u);
extern double N1(double u);
extern double N2(double u);
extern double N3(double u);
extern double N4(double u);
extern double N5(double u);
extern double N6(double u);

extern double dN0(double u);
extern double dN1(double u);
extern double dN2(double u);
extern double dN3(double u);
extern double dN4(double u);
extern double dN5(double u);
extern double dN6(double u);
}

double N(int d, int shift, double u) { 
  assert(d >= 0 && d <=6);
  switch(d) {
  case 0: return N0 (u+d+1-shift);
  case 1: return N1 (u+d+1-shift);
  case 2: return N2 (u+d+1-shift);
  case 3: return N3 (u+d+1-shift);
  case 4: return N4 (u+d+1-shift);
  case 5: return N5 (u+d+1-shift);
  case 6: return N6 (u+d+1-shift);
  }
}

double dN(int d, int shift, double u) { 
  assert(d >= 0 && d <=6);
  switch(d) {
  case 0: return dN0 (u+d+1-shift);
  case 1: return dN1 (u+d+1-shift);
  case 2: return dN2 (u+d+1-shift);
  case 3: return dN3 (u+d+1-shift);
  case 4: return dN4 (u+d+1-shift);
  case 5: return dN5 (u+d+1-shift);
  case 6: return dN6 (u+d+1-shift);
  }
}


double evalSpline(int chart_type, int d, int valence, double u, double v) { 
  int ncp = d+NUM_PATCHES;  
  double res = 0.0; 
  for(int j = 0; j < ncp; j++)
    for(int k = 0; k < ncp; k++) { 
      res  += CPs[chart_type][d][valence][j*ncp+k]*
	N(d,j,u)*N(d,k,v);
    }

  return res; 
}


void evalSplinesSamples() { 
  for(int chart_type = 0; chart_type <= MAX_CHART_TYPE; chart_type++)
    for(int d = 1; d <= MAX_D; d++)
      for(int valence = 3; valence <= MAX_VALENCE; valence++) {
	cout << chart_type << " " << d << " " << valence << endl;
	for(int i = 0; i < 20*valence + ExtrPts[chart_type]; i++) { 
	  double u =   (NUM_PATCHES-1)*(SamplePts[chart_type][valence][2*i]  -MinX[chart_type][valence])/(MaxX[chart_type][valence]-MinX[chart_type][valence]);
	  double v =   (NUM_PATCHES-1)*(SamplePts[chart_type][valence][2*i+1]-MinY[chart_type][valence])/(MaxY[chart_type][valence]-MinY[chart_type][valence]);
	  if (chart_type == CHART_CONCAVE) v = (NUM_PATCHES-1)-v;
	  cout << evalSpline(chart_type,d, valence,u,v ) << " ";		  
	}
	cout << endl << endl;
      }
}

void evalSplinesRegular(int num_grid) { 
  for(int chart_type = 0; chart_type <= MAX_CHART_TYPE; chart_type++)
    for(int d = 1; d <= MAX_D; d++)
      for(int valence = 3; valence <= MAX_VALENCE; valence++) {
	cout << chart_type << " " << d << " " <<  valence << endl;
	for(int j = 0; j < num_grid; j++) 
	  for(int k = 0; k < num_grid; k++) {
	    cout << evalSpline(chart_type,d,valence, 
			       (NUM_PATCHES-1)*j/double(num_grid-1),
			       (NUM_PATCHES-1)*k/double(num_grid-1)
			       ) << " " ;	
	  }
	cout << endl << endl;
      }
  
}


#include "vec3t.hpp"
#include "numvector.hpp"

Vector<Point3>* constructSplineCPfromW( int chart_type, int degree, int valence,  
					const vector<Point3>& samples) { 

  cerr << "constr CP " << chart_type << " " << degree << " " << valence  << endl;
  int d = degree-1;
  int ncp = d + NUM_PATCHES;
  Vector<Point3>* CPs = new Vector<Point3>(ncp*ncp);

  if( 0 ) { //chart_type == 3 && valence == 7 ) { 
  cout <<  endl << "samples = [";
  for(int i=0; i< samples.size(); i++)  { cout << "[" << samples[i](0) << "," 
  					     << samples[i](1) << ","
  					     << samples[i](2)  << "]";
    if( i < samples.size()-1) cout << ",";
  }
  cout << "];" << endl;
  
  }

  for(int l = 0; l < ncp*ncp; l++) {
    (*CPs)(l) = Point3(0.0);
    for(int i = 0; i < 20*valence + ExtrPts[chart_type]; i++) {
      (*CPs)(l) += 
	W[chart_type][d][valence][ (20*valence + ExtrPts[chart_type])*l+i]*
	samples[i];
    }
  }
  return CPs;
}



void evalSplineFromCP(int chart_type, int degree, int valence, double* xy, 
			     const Vector<Point3>& CPs, Point3* ret) { 

  int d = degree -1;
  int ncp = d+NUM_PATCHES;  
  Point3 res = Point3(0.0); 
   

  double u =   (NUM_PATCHES-1)*(xy[0] - MinX[chart_type][valence])/(MaxX[chart_type][valence]-MinX[chart_type][valence]);
  double v =   (NUM_PATCHES-1)*(xy[1] - MinY[chart_type][valence])/(MaxY[chart_type][valence]-MinY[chart_type][valence]);  
  if(chart_type == CHART_CONCAVE) v = NUM_PATCHES-1-v;


  for(int j = 0; j < ncp; j++)
    for(int k = 0; k < ncp; k++) { 
      res  += CPs(j*ncp+k)*N(d,j,u)*N(d,k,v);
    }

  *ret = res;
}

void evalSplineDerFromCP(int chart_type, int degree, int valence, double* xy, 
			 const Vector<Point3>& CPs, Point3* ret) { 

  int d = degree-1;
  int ncp = d+NUM_PATCHES;  
  Point3 res_u = Point3(0.0); 
  Point3 res_v = Point3(0.0); 
  double u =   (NUM_PATCHES-1)*(xy[0] - MinX[chart_type][valence])/(MaxX[chart_type][valence]-MinX[chart_type][valence]);
  double v =   (NUM_PATCHES-1)*(xy[1] - MinY[chart_type][valence])/(MaxY[chart_type][valence]-MinY[chart_type][valence]);  

  if(chart_type == CHART_CONCAVE) v = NUM_PATCHES-1-v;
  for(int j = 0; j < ncp; j++)
    for(int k = 0; k < ncp; k++) { 
      res_u  += (dN(d,j,u)*N(d,k,v)*(NUM_PATCHES-1)/(MaxX[chart_type][valence]-MinX[chart_type][valence]))*CPs(j*ncp+k);
      res_v  += (N(d,j,u)*dN(d,k,v)*(NUM_PATCHES-1)/(MaxY[chart_type][valence]-MinY[chart_type][valence]))*CPs(j*ncp+k);	       }

  if(chart_type == CHART_CONCAVE) res_v = -res_v;
  ret[0] = res_u; 
  ret[1] = res_v; 
}
