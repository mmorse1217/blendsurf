//catmull-clark surface

#include "ccsurf.hpp"
#include "ccsurfop.hpp"
#include "vecmatop.hpp"

// ---------------------------------------------------------------------- 
int CCSurf::setup(double _iniflatS, istream& sin, int flipNormal, Point3 scale, void* _iniSubdivMatrices)
{
  
 
    //0. set _subdivMatrices
    _subdivMatrices = (CCSubMatLib*)_iniSubdivMatrices;
    _flatS = _iniflatS;  

    //1. gpmesh initialize
    _gpmesh.setup(sin, flipNormal, scale);

    _pos.resize(MAXLEVELNUM); for(int i=0; i<MAXLEVELNUM; i++) _pos[i].resize(numFaces());
  
    //2. add gpmesh points
    DataLevel& tp = _pos[0];
    tp.resize(numFaces());
  
    for(int F=0; F<numFaces(); F++) {
        tp[F].resize(2,2);
        //this numbering convention is good to keep in mind.
        tp[F](0,0) = _gpmesh.vpoint( _gpmesh.Fv2Vf(F,0).first );
        tp[F](1,0) = _gpmesh.vpoint( _gpmesh.Fv2Vf(F,1).first );
        tp[F](1,1) = _gpmesh.vpoint( _gpmesh.Fv2Vf(F,2).first );
        tp[F](0,1) = _gpmesh.vpoint( _gpmesh.Fv2Vf(F,3).first );
    }

    //LEXING: VERY IMPORTANT
    CCSurfOp::share(_gpmesh, 0, _pos[0]);

    _inlvl = 0;
    assert(validpos(0));
  
    //3. extra reading
    //details are for multiresolution surfaces.
    char tmp[100];
    sin>>tmp;
    if(sin.eof()==false) { //not end yet
        assert(strcmp(tmp,"detail")==0);
        int lvl;
        sin>>lvl; assert(lvl>=0 && lvl<MAXLEVELNUM);
        int size = pow2(lvl)+1;
        DataLevel& tp = _pos[lvl];
        for(int F=0; F<numFaces(); F++) {
            tp[F].resize(size,size);
            for(int j=0; j<size; j++)
                for(int i=0; i<size; i++) {
                    double x,y,z;
                    sin>>x>>y>>z;
                    tp[F](i,j) = Point3(x,y,z);
                }
        }
        //LEXING: VERY IMPORTANT
        CCSurfOp::share(_gpmesh, lvl, _pos[lvl]);
        _inlvl = lvl;
        assert(validpos(lvl));
    }
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurf::dump(ostream& fout, int lvl)
{
    //dump data a certain lvl
    assert( validpos(lvl) );
    int size = pow2(lvl)+1;
    vector<Point3> ptvec;
    vector< CCRect<int> > il(numFaces()); //index level
    DataLevel& pl = _pos[lvl];
  
    //0. clear index level
    for(int F=0; F<numFaces(); F++) {
        il[F].resize(size, size);
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++)
                il[F](i,j) = -1;
    }
    int cnt = 0;
    GpMesh& gm = _gpmesh;
    //1. do vertices
    for(int V=0; V<numVertices(); V++) {
        int K = gm.valence(V);
        Point3 tmp;
        for(int f=0; f<K; f++) {
            intpair Fv = gm.Vf2Fv(V,f);
            int i,j;
            CCSurfOp::eno2index(lvl, Fv.second, 0, i, j);
            il[Fv.first](i,j) = cnt;
            tmp = pl[Fv.first](i,j);
        }
        // store
        ptvec.push_back(tmp);
        cnt ++;
    }
    //2. do edge
    for(int F=0; F<numFaces(); F++) {
        for(int e=0; e<4; e++) {
            intpair nFe = gm.Fe2Fe(F,e);
            int nF = nFe.first; 
            int ne = nFe.second;
	    if (nF == -1){
	      nF = F; ne = e;
	    }
	    
            else if(F>nF) continue; //skip if already done
            int ci, cj; CCSurfOp::eno2index(lvl, e, 0, ci, cj);
            int cdi,cdj;CCSurfOp::eno2index(lvl, e, 1, cdi,cdj);
            cdi=cdi-ci; cdj=cdj-cj; //difference
            int ni, nj; CCSurfOp::eno2index(lvl,ne,size-1, ni, nj);
            int ndi,ndj;CCSurfOp::eno2index(lvl,ne,size-2,ndi,ndj);
            ndi=ndi-ni; ndj=ndj-nj;
            for(int k=1; k<size-1; k++) {
                il[ F](ci+k*cdi,cj+k*cdj) = cnt;
                il[nF](ni+k*ndi,nj+k*ndj) = cnt;
                // store
                ptvec.push_back(pl[nF](ni+k*ndi,nj+k*ndj));
                cnt ++;
            }
        }
    }
    //3. do face
    for(int F=0; F<numFaces(); F++) {
        for(int i=1; i<size-1; i++)
            for(int j=1; j<size-1; j++) {
                il[F](i,j) = cnt;
                //store
                ptvec.push_back(pl[F](i,j));
                cnt ++;
            }
    }
  
    //4. write down 
    fout<<"points "<<ptvec.size()<<endl;
    for(int p=0; p<ptvec.size(); p++)
        fout<<ptvec[p]<<endl;
    fout<<"coords "<<numFaces()*(size-1)*(size-1)<<endl;
    for(int F=0; F<numFaces(); F++) {
        for(int i=0; i<size-1; i++)
            for(int j=0; j<size-1; j++) {
                fout<<il[F](i,j)<<" "<<il[F](i+1,j)<<" "<<il[F](i+1,j+1)<<" "<<il[F](i,j+1)<<endl;
            }
    }
    
    fout << endl<<"gids "<<numFaces()*(size-1)*(size-1)<<endl;
    // TODO fix for multiply connected domains
    for(int F=0; F<numFaces(); F++) {
        for(int i=0; i<size-1; i++)
            for(int j=0; j<size-1; j++) {
                fout << 0 << endl;
            }
    }
    fout << endl;
    fout << "groups " << 1 << endl;
    fout << "intpt " << 0 << " " << 0 << " " << 0 << endl;
    fout << "orient " << 1 << endl;
  
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurf::rotate(const Quaternion& r)
{
  
    for(int i=0; i<MAXLEVELNUM; i++) {
        if(validpos(i)) {		
            CCSurfOp::rotate(_gpmesh, i, _pos[i], r);
        }
    
    }
    return 0;
}
// ---------------------------------------------------------------------- 
int CCSurf::shift(const Point3& s)
{
  
    for(int i=0; i<MAXLEVELNUM; i++) {
        if(validpos(i)) {		CCSurfOp::shift(_gpmesh, i, _pos[i], s);	 }
    }
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurf::eval(int flags, int F, double* cd, Point3* ret)
{
  
    double c=cd[0];  double d=cd[1];
    int lastlvl = lastposlvl(); 
    int sz = pow2(lastlvl);
    double ss = 1.0/sz;
    int v = 0;
    bool adj = false;
    if(       c<ss && d<ss) {
        adj = true; v = 0;
    } else if(c>1-ss && d<ss) {
        adj = true; v = 1;
    } else if(c>1-ss && d>1-ss) {
        adj = true; v = 2;
    } else if(c<ss && d>1-ss) {
        adj = true; v = 3;
    } else{
        adj = false;
    }
    bool bad;
    if(adj==false) {
        bad = false;
    } else {
        intpair Vf = _gpmesh.Fv2Vf(F,v);
        int V = Vf.first;
        if(_gpmesh.valence(V) == 4) {
            bad = false;
        } else {
            bad = true;
        }
    }
    //LEXING: FOR THE TIME BEING, this is good enough
    assert(bad==false);
    CCRect<Point3> PD = _pos[lastlvl][F];  //Matrix<Point3> mp(sz+3, sz+3, false, PD.data()); //matrix
    int    mn[2];  
    mn[0] = sz+3;  mn[1] = sz+3;
    double ef[2];  
    ef[0] = double(sz+3)/double(sz);  ef[1] = double(sz+3)/double(sz);
    double ts[2];  
    ts[0] = cd[0]*sz;  ts[1] = cd[1]*sz;
    int    ij[2];  
    ij[0] = min(max(int(floor(ts[0])),0),sz-1);  
    ij[1] = min(max(int(floor(ts[1])),0),sz-1);
    double lf[2];  
    lf[0] = ts[0]-ij[0];	 lf[1] = ts[1]-ij[1];
    ij[0]++;  ij[1]++; 
    //CALL spline evaluation function
    spev2d(flags, DMFLAG_CLOSED, 3, (double*)(PD.data()), mn, ef, ij, lf, (double*)ret);
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurf::subdivide(int lvl, CCSubMatLib& sb)
{
  
    assert(lvl<MAXLEVELNUM-1);
    CCSurfOp::subdivide(_gpmesh, lvl,_flatS,  _pos[lvl], _pos[lvl+1], sb );
    return 0;
}
