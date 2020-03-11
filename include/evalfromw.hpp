#ifndef __EVALFROMW_H__


extern void loadAllMatrices();

extern Vector<Point3>* constructSplineCPfromW( int chart_type, int degree, int valence,  
					       const vector<Point3>& samples);

extern void evalSplineFromCP(int chart_type, int degree, int valence, double* xy, 
		       const Vector<Point3>& CPs, Point3* ret);			 


extern void evalSplineDerFromCP(int chart_type, int degree, int valence, double* xy, 
		       const Vector<Point3>& CPs, Point3* ret);			 

extern string MatricesDir;
extern int NUM_PATCHES;

#endif
