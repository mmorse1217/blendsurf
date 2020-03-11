#ifndef _COMMON_HPP_
#define _COMMON_HPP_

// define input files
#ifdef __DECCXX
#define __USE_STD_IOSTREAM
#endif

#if defined(__GNUG__) || defined(_STANDARD_C_PLUS_PLUS) || defined(__DECCXX) || defined(sun) || defined(__linux)
#include <iostream>
#include <fstream>
#include <sstream>
#else
#include <iostream.h>
#include <fstream.h>
#include <sstream.h>
#endif

#include <cfloat>
#include <cassert>
#include <cmath>
#include <cstring>
#include <memory>
#include <complex>

#include <vector>
#include <set>
#include <map>
#include <deque>
#include <queue>
#include <utility>
#include <algorithm>

using namespace std;

// make some global definitions
#define assertMsg(expr,msg) assert(expr)
#define warningMsg(msg) cerr<<msg<<endl

enum DimType {
    DIM0=0,
    DIM2=2,
    DIM3=3
};

#define SCL_EPS DBL_EPSILON
#define SCL_MAX DBL_MAX
#define SCL_MIN DBL_MIN

#define MAXDIM 3
#define PI M_PI

// define some global functions
inline DimType min(DimType a, DimType b) { return DimType(a&b); }
// ---------------------------------------------------------------------- 
inline void error(const char* msg) { cerr<<"ERROR "<<msg<<endl; }
inline int pow2(int l) { assert(l>=0); return (1<<l); }
// ---------------------------------------------------------------------- 
template <class F> inline int round(F f)  { return int(floor(f+0.5)); }

#endif //_COMMON_HPP_

