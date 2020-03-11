//CCsurfop: catmull clark subdivion surface operations.

#include "ccsurfop.hpp"
#include "mat3t.hpp"

// ---------------------------------------------------------------------- 
int CCSurfOp::share(GpMesh& gm, int lvl, DataLevel& wk)
{//this is to keep track of neighbors of a face before/after subdivision.
    //each face has a rectangle associated with it that has an extra column/row
    //(at -1 and size)

  
    int numf = gm.numFaces();
    int size = pow2(lvl)+1;
    //--------------------
    //sharing edges : neighbors of the edge of a given face.
    for(int F=0; F<numf; F++) {
        for(int e=0; e<4; e++) {
            //for each (face,edge) double

            //return the indices of start point of edge
            int ci, cj; eno2index(lvl,e,0,ci,cj);
            //return the index of the next point on the edge
            int cdi,cdj;eno2index(lvl,e,1,cdi,cdj);
            cdi=cdi-ci; cdj=cdj-cj; //now (cdi, cdj) defines a vector      
            //perform necessary transformation to get to the correct position
            //in th rectangle.
            ci=ci+cdj; cj=cj-cdi;
            //have to end up at row/col -1 or size
            assert( (ci==-1)||(cj==-1)||(ci==size)||(cj==size));

            //retrieve the rectangle for current face
            CCRect<Point3>& CD = wk[F];      
            //find the neigboring face and edge
            intpair nFe = gm.Fe2Fe(F,e);
            int nF = nFe.first;
            int ne = nFe.second;

            if (nF != -1){ //interior face
                //similarly find the correct position in the neiboring face
                //rectangle and copy values to current.
                int ni, nj; eno2index(lvl,ne,size-1,ni,nj);
                int ndi,ndj;eno2index(lvl,ne,size-2,ndi,ndj);
                ndi=ndi-ni; ndj=ndj-nj;
                ni=ni+ndj; nj=nj-ndi;
                assert( (ni==1)||(nj==1)||(ni==size-2)||(nj==size-2));
                CCRect<Point3>& ND = wk[nF];
                for(int k=0; k<size; k++)
                    CD(ci+k*cdi,cj+k*cdj) = ND(ni+k*ndi,nj+k*ndj);
            }
            else{//face with an edge on boundary
                //if we're on the edge , want to share with the face itself
                //s.t. the position we want to fill should have value 2a-b
                //where a is the value on the boundary and b is the one right inside
                //of that. this way the average would give the value on boundary 
                //((2a-b+b)/2 = a)

                if (ci==-1)
                    for(int k=0; k<size; k++)
                        CD(ci+k*cdi,cj+k*cdj) = 2.0*CD(ci+k*cdi+1, cj+k*cdj) - 
                            CD(ci+k*cdi+2, cj+k*cdj);
	
                if (ci == size)
                    for(int k=0; k<size; k++)
                        CD(ci+k*cdi,cj+k*cdj) = 2.0*CD(ci+k*cdi-1, cj+k*cdj) - 
                            CD(ci+k*cdi-2, cj+k*cdj);
	
                if (cj == -1)
                    for(int k=0; k<size; k++)
                        CD(ci+k*cdi,cj+k*cdj) = 2.0*CD(ci+k*cdi, cj+k*cdj+1) - 
                            CD(ci+k*cdi, cj+k*cdj+2);
	
                if (cj == size)
                    for(int k=0; k<size; k++)  
                        CD(ci+k*cdi,cj+k*cdj) = 2.0*CD(ci+k*cdi, cj+k*cdj-1) - 
                            CD(ci+k*cdi, cj+k*cdj-2);
            }
        } 
    }
  

    //--------------------
    //sharing vertices
    //need to share the vertices that correspond to the corners in the rectangle
    for(int F=0; F<numf; F++) {
        for(int e=0; e<4; e++) {
            //initially proceed as in the edge case. each "edge" will be responsible 
            //for one vertex in ccw(?) direction
            int ci, cj; eno2index(lvl,e,0,ci,cj);
            int cdi,cdj;eno2index(lvl,e,1,cdi,cdj);
            cdi=cdi-ci; cdj=cdj-cj;
            ci=ci+cdj-cdi; cj=cj-cdi-cdj;
            CCRect<Point3>& CD = wk[F];
            intpair nFe = gm.Fe2Fe(F,e);
            int nF = nFe.first;
            int ne = nFe.second;
            if (nF!=-1){ //if there exists a face on the other side
                //beginning point of one edge is the end of the corresponding edge
                //on neighboring face. so retrieve indeces accordingly.
                int ni, nj; eno2index(lvl,ne,size-1,ni,nj);
                int ndi,ndj;eno2index(lvl,ne,size-2,ndi,ndj);
                ndi=ndi-ni; ndj=ndj-nj;
                //transform to point to correct vertex 
                ni=ni+ndj-ndi; nj=nj-ndi-ndj;
                CCRect<Point3>& ND = wk[nF];
                CD(ci,cj) = ND(ni,nj);
            }
            else{ //boundary face
                //since the current edge is on the boundary look at the face next
                //to current one. since the corner of a boundary face is equal to an 
                //edge share value of the neigboring face.
                intpair nFe = gm.Fe2Fe(F, gm.preveno(e));
                int nF = nFe.first; int ne =nFe.second;
                if (nF == -1){
		  nFe = gm.Fe2Fe(F, gm.nexteno(e));
		  nF = nFe.first; ne = nFe.second;
		}
		int ni, nj; eno2index(lvl, ne,0,ni,nj);
                int ndi,ndj;eno2index(lvl, ne,1,ndi,ndj);	
                ndi=ndi-ni; ndj=ndj-nj;
                ni=ni-ndj-ndi;   
                nj=nj-ndj+ndi;
                CCRect<Point3>& ND = wk[nF];
                CD(ci, cj) = ND(ni, nj);	
            }
        }
    }

    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::midsub(GpMesh&, int, DataLevel&, DataLevel&)
{
  
    //TODO
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::subdivide(GpMesh& gm, int lvl, double flats, DataLevel& fm, DataLevel& to, CCSubMatLib& submat)
{
  
  
    int numf = gm.numFaces();
    int numv = gm.numVertices();
    int size = pow2(lvl)+1;
  
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(2*size-1,2*size-1);
    }

    //subdivide all faces and edges according to regular cc rules
    for (int F=0; F<numf; F++){
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];

        {//subdivide face
            double cwt = 1.0/64.0;
            double ewt = 6.0/64.0;
            double swt = 36.0/64.0;
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    ND(2*i,2*j) = 
                        cwt*CD(i-1,j-1) + ewt*CD(i,j-1) + cwt*CD(i+1,j-1) +
                        ewt*CD(i-1,j  ) + swt*CD(i,j  ) + ewt*CD(i+1,j  ) +
                        cwt*CD(i-1,j+1) + ewt*CD(i,j+1) + cwt*CD(i+1,j+1);
        }
        {//subdivide edge   
            double cwt = 1.0/16.0;
            double ewt = 6.0/16.0;
            for(int i=0; i<size-1; i++)
                for(int j=0; j<size; j++)
                    ND(2*i+1,2*j) = 
                        cwt*CD(i,j-1) + cwt*CD(i+1,j-1) +
                        ewt*CD(i,j  ) + ewt*CD(i+1,j  ) +
                        cwt*CD(i,j+1) + cwt*CD(i+1,j+1);
            for(int i=0; i<size; i++)
                for(int j=0; j<size-1; j++)
                    ND(2*i,2*j+1) = 
                        cwt*CD(i-1,j  ) + ewt*CD(i,j  ) + cwt*CD(i+1,j  ) +
                        cwt*CD(i-1,j+1) + ewt*CD(i,j+1) + cwt*CD(i+1,j+1);
        }
        {//subdivide face
            double cwt = 1.0/4.0;
            for(int i=0; i<size-1; i++)
                for(int j=0; j<size-1; j++)
                    ND(2*i+1,2*j+1)= 
                        cwt*( CD(i,j  ) + CD(i+1,j  ) +
                              CD(i,j+1) + CD(i+1,j+1) );
        }
    }


   //subdivide vertex
    for(int V=0; V<numv; V++) {
      
        int isboun = gm.boun(V);
        int K = gm.valence(V);
	//if (K==1) cout<<V<<endl;

	//if (isboun == 1) assert(K>=2 || K==0);
        //else assert(K>=3 || K==0);
	
        if(K==0) continue; //unused vertices
        if(isboun == 0 && K!=4) {
            //if not on the boundary
            vector<Point3> onering;
            //get the ring around the vertex
            getonering(gm, lvl, fm, V, onering);
            Point3 center(0.0);
            //choose the relevant matrix from library
            NumMatrix& SM = submat.chooseMatrix(CCSubMatLib::INTERIOR, K).SM();
            for(int j=0; j<onering.size(); j++) {
                //find position of the center vertex
                center += SM(0,j) * onering[j];
            }
            putcenter(gm, lvl+1, to, V, center );
        }
        else if (isboun == 1 &&gm.corner(V) !=2){
            //if on the boundary the previous rules don't hold anymore so need to 
            //get the halfring and replace it with the full half ring rather than 
            //just the center vertex.
            vector<Point3> halfring;
            vector<Point3> newhalfring;
            gethalfring(gm, lvl, fm, V, halfring);
            newhalfring.resize(halfring.size());
	 
            for(int i=0; i<halfring.size(); i++)
                newhalfring[i] = (0.0);
            //choose relevant matrix from library
            NumMatrix SM;
            if (gm.corner(V) == 0)
                SM = submat.chooseMatrix(CCSubMatLib::BOUNDARY, K).SM();
            else if (gm.corner(V) == 1)
                SM = submat.chooseMatrix(CCSubMatLib::CORNERCONVEX, K).SM();
            else if (gm.corner(V) == 2)
                SM = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).SM();
            //cout<<SM.m()<<" "<<SM.n()<<endl;
            for(int j=0; j<halfring.size(); j++) {	
                for(int k=0; k<halfring.size(); k++)
                    newhalfring[j] += SM(j,k) * halfring[k];
            }

            /*
              if (V==1){
              for (int i=0; i<newhalfring.size();i++)
              cout<<newhalfring[i]<<endl;
              }
            */
            //put back...
            puthalfring(gm, lvl+1, to, V, newhalfring );   
        }
    }
    

    for(int V=0; V<numv; V++) {    

        int isboun = gm.boun(V);
        int K = gm.valence(V);
        if (isboun ==1 && gm.corner(V) == 2){
   
            vector<Point3> halfring;
            vector<Point3> newhalfring;
            gethalfring(gm, lvl, fm, V, halfring);
      
            newhalfring.resize(halfring.size());
            for(int i=0; i<halfring.size(); i++)
                newhalfring[i] = (0.0);
      
            NumMatrix &SM = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).SM();

            for(int j=0; j<halfring.size(); j++) {	
                for(int k=0; k<halfring.size(); k++)
                    newhalfring[j] += SM(j,k) * halfring[k];
            }

      
            NumMatrix &EV = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).V();
            NumMatrix &IV = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).IV();
            NumVector &D = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).D();
      
            Point3 a0(0.0);
            Point3 a2(0.0);
            Point3 a3(0.0);
      
            assert(D(2) == D(3) && D(2) == .5);
      
            for(int j=0; j<newhalfring.size(); j++) {
                a0 += IV(0,j) * newhalfring[j];
                a2 += IV(2,j) * newhalfring[j];
                a3 += IV(3,j) * newhalfring[j];
            }
      
      
            if (flats == -1)
                flats = 1.0-(1.0/(2.0+cos(3.0*PI/(2.0*K))-cos(3.0*PI/2.0)));
            //flats  = (2*D(1) - 1)/(2*D(1));
      
      
            for (int j=0; j<newhalfring.size();j++){
                newhalfring[j] = (1-flats)*newhalfring[j];
                newhalfring[j] += flats*(a0 + EV(j,2)*a2 + EV(j, 3)*a3);
            }
      
            puthalfring(gm, lvl+1, to, V, newhalfring );   
        }
    
    }  
    //LEXING: VERY IMPORTANT
    share(gm, lvl+1, to);
 
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::limit(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to, CCSubMatLib& submat)
{

  
  
    int numf = gm.numFaces();
    int numv = gm.numVertices();

    int size = pow2(lvl)+1;

    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }
  
    //compute limit regularly
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];

        double cwt = 1.0/36.0;
        double ewt = 4.0/36.0;
        double swt = 16.0/36.0;
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++)
                ND(i,j)=
                    cwt*CD(i-1,j-1) + ewt*CD(i,j-1) + cwt*CD(i+1,j-1) +
                    ewt*CD(i-1,j  ) + swt*CD(i,j  ) + ewt*CD(i+1,j  ) +
                    cwt*CD(i-1,j+1) + ewt*CD(i,j+1) + cwt*CD(i+1,j+1);
    }
  
    for(int V=0; V<numv; V++) {
        int isboun = gm.boun(V);
        int K = gm.valence(V);
        if (isboun == 0)  assert(K>=3 || K==0);
        //else assert(K>=2 || K==0);
        if(K==0) continue; //unused vertices
        if(isboun == 0 &&  K!=4) {
            //interior vertex
            vector<Point3> onering;
            getonering(gm, lvl, fm, V, onering);
            Point3 limit( 0.0);
            NumMatrix& IV = submat.chooseMatrix(CCSubMatLib::INTERIOR, K).IV();
     

            for(int j=0; j<onering.size(); j++) {
                limit += IV(0,j) * onering[j];
            }
            putcenter(gm, lvl, to, V, limit);
        }
        else if(isboun == 1){
            //boundary vertex
            vector<Point3> halfring;
            gethalfring(gm, lvl, fm, V, halfring);
            Point3 limit( 0.0);

            if (gm.corner(V) == 0){
                NumMatrix &IV = submat.chooseMatrix(CCSubMatLib::BOUNDARY, K).IV();
	
                for(int j=0; j<halfring.size(); j++) {
                    limit += IV(0,j) * halfring[j];
                }      
            }
            else if(gm.corner(V) == 1){
                NumMatrix &IV = submat.chooseMatrix(CCSubMatLib::CORNERCONVEX, K).IV();
	
                for(int j=0; j<halfring.size(); j++) {
                    limit += IV(0,j) * halfring[j];
	 
                }  
            }
            else if(gm.corner(V) == 2){
                NumMatrix &IV = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).IV();
       
                for(int j=0; j<halfring.size(); j++) {
                    limit += IV(0,j) * halfring[j];
                }  
            }


            //here don't need to replace the whole half ring.
            putcenter(gm, lvl, to, V, limit);
      
        }
    }
  
    //LEXING: VERY IMPORTANT
    share(gm, lvl, to);
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::normal(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to, CCSubMatLib& submat)
{

  
  
    int numf = gm.numFaces();
    int numv = gm.numVertices();
    int size = pow2(lvl)+1;
  
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }

    //normals are pretty much the same just replace one ring with half ring on the
    //boundary.

    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];
    
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++) {
                Point3 du = 
                    ( (CD(i+1,j  )-CD(i-1,j  ))*4.0 + 
                      (CD(i+1,j+1)-CD(i-1,j+1)) +
                      (CD(i+1,j-1)-CD(i-1,j-1)) )  / 12.0;
                Point3 dv = 
                    ( (CD(i,  j+1)-CD(i,  j-1)  )*4.0 + 
                      (CD(i+1,j+1)-CD(i+1,j-1)) +
                      (CD(i-1,j+1)-CD(i-1,j-1)) ) / 12.0;
                Point3 nor( cross(du,dv) ); nor=nor.dir();
                ND(i,j) = nor;
            }
    }
    //corners

    for(int V=0; V<numv; V++) {
        int isboun = gm.boun(V);
        int K = gm.valence(V);	
        if (isboun==0) assert(K>=3 || K==0);
        //else assert(K>=2 || K==0);
        if(K==0) continue;
        if(isboun == 0 && K!=4) {
            vector<Point3> onering;
            getonering(gm, lvl, fm, V, onering);
            NumMatrix& IV = submat.chooseMatrix(CCSubMatLib::INTERIOR, K).IV();

            Point3 du( 0.0);
            for(int j=0; j<onering.size(); j++) {
                du += IV(1,j) * onering[j];
            }
            Point3 dv( 0.0);
            for(int j=0; j<onering.size(); j++) {
                dv += IV(2,j) * onering[j];
            }
            Point3 nor( cross(du,dv) ); nor=nor.dir();
            putcenter(gm, lvl, to, V, nor);
        }
        else if (isboun == 1){
            vector<Point3> halfring;
            gethalfring(gm, lvl, fm, V, halfring);
      

            if (gm.corner(V) == 0){
                NumMatrix &IV = submat.chooseMatrix(CCSubMatLib::BOUNDARY, K).IV();
                Point3 du( 0.0);
                for(int j=0; j<halfring.size(); j++) {
                    du += IV(1,j) * halfring[j];
                }
                Point3 dv( 0.0);
                for(int j=0; j<halfring.size(); j++) {
                    dv += IV(2,j) * halfring[j];
                }
                Point3 nor( cross(du,dv) ); nor=nor.dir();
                putcenter(gm, lvl, to, V, nor);
            } 
      
            else if (gm.corner(V) == 1){
                NumMatrix &IV  = submat.chooseMatrix(CCSubMatLib::CORNERCONVEX, K).IV();
                Point3 du( 0.0);
                for(int j=0; j<halfring.size(); j++) {
                    du += IV(1,j) * halfring[j];
                }
                Point3 dv( 0.0);
                for(int j=0; j<halfring.size(); j++) {
                    dv += IV(2,j) * halfring[j];
                }
                Point3 nor( cross(du,dv) ); nor=nor.dir();
                putcenter(gm, lvl, to, V, nor);
            }
      
            else if (gm.corner(V) == 2){
                NumMatrix &IV = submat.chooseMatrix(CCSubMatLib::CORNERCONCAVE, K).IV();
	
                Point3 du( 0.0);
                for(int j=0; j<halfring.size(); j++) {
                    du += IV(2,j) * halfring[j];
                }
                Point3 dv( 0.0);
                for(int j=0; j<halfring.size(); j++) {
                    dv += IV(3,j) * halfring[j];
                }
                Point3 nor( cross(du,dv) ); nor=nor.dir();
                putcenter(gm, lvl, to, V, nor);
            }
      

      

        }
    }
  
    //LEXING: VERY IMPORTANT
    share(gm, lvl, to);
    return 0;
}

// ---------------------------------------------------------------------- 
//derivative functions are the same for interior and boundary points. 

int CCSurfOp::du(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
 
  
    int numf = gm.numFaces();
    // int numv = gm.numVertices();
    int size = pow2(lvl)+1;
    double h = 1.0/pow2(lvl);
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];
    
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++) {
                ND(i,j) = 
                    ( (CD(i+1,j  )-CD(i-1,j  ))*4.0 + 
                      (CD(i+1,j+1)-CD(i-1,j+1)) +
                      (CD(i+1,j-1)-CD(i-1,j-1)) )  / (12.0*h);
            }
    }
    //LEXING: IMPORTANT, DO NOT SHARE
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::dv(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{

  
    int numf = gm.numFaces();
    //int numv = gm.numVertices();
    int size = pow2(lvl)+1;
    double h = 1.0/pow2(lvl);
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];
    
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++) {
                ND(i,j) = 
                    ( (CD(i,  j+1)-CD(i,  j-1)  )*4.0 + 
                      (CD(i+1,j+1)-CD(i+1,j-1)) +
                      (CD(i-1,j+1)-CD(i-1,j-1)) ) / (12.0*h);
            }
    }
    //LEXING: IMPORTANT, DO NOT SHARE
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::duu(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{

  
    int numf = gm.numFaces();  
    // int numv = gm.numVertices();
    int size = pow2(lvl)+1;  
    double h = 1.0/pow2(lvl);
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++) {
                ND(i,j) = 
                    ( (CD(i+1,j+1)-2.*CD(i,j+1)+CD(i-1,j+1)) +
                      (CD(i+1,j  )-2.*CD(i,j  )+CD(i-1,j  ))*4.0 + 
                      (CD(i+1,j-1)-2.*CD(i,j-1)+CD(i-1,j-1)) )  / (6.0*h*h);
            }
    }
    //LEXING: IMPORTANT, DO NOT SHARE
    return 0;
}
// ---------------------------------------------------------------------- 
int CCSurfOp::dvv(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{

  
    int numf = gm.numFaces();  
    //  int numv = gm.numVertices();
    int size = pow2(lvl)+1;  
    double h = 1.0/pow2(lvl);
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++) {
                ND(i,j) = 
                    ( (CD(i+1,j+1)-2.*CD(i+1,j)+CD(i+1,j-1)) + 
                      (CD(i,  j+1)-2.*CD(i  ,j)+CD(i,  j-1)  )*4.0 + 
                      (CD(i-1,j+1)-2.*CD(i-1,j)+CD(i-1,j-1)) ) / (6.0*h*h);
            }
    }
    //LEXING: IMPORTANT, DO NOT SHARE
    return 0;
}
// ---------------------------------------------------------------------- 
int CCSurfOp::duv(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{

  
    int numf = gm.numFaces();  
    //int numv = gm.numVertices();
    int size = pow2(lvl)+1;  double h = 1.0/pow2(lvl);
    assert(fm.size()==numf);
    to.resize(numf);
    for(int F=0; F<numf; F++) {
        assert(fm[F].m()==size && fm[F].n()==size);
        to[F].resize(size,size); 
    }
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = fm[F];
        CCRect<Point3>& ND = to[F];
        for(int i=0; i<size; i++)
            for(int j=0; j<size; j++) {
                ND(i,j) = 
                    ( (CD(i+1,j+1)-CD(i+1,j-1)) -
                      (CD(i-1,j+1)-CD(i-1,j-1)) ) / (4.0*h*h);
            }
    }
    //LEXING: IMPORTANT, DO NOT SHARE
    return 0;
}
// ---------------------------------------------------------------------- 
int CCSurfOp::rotate(GpMesh& gm, int lvl, DataLevel& wk, Quaternion qr)
{
  
    NumMatrix tmp(3,3);  
    qr.toRotationMatrix(tmp);
    Matrix3 rot(tmp.data());
    int numf = gm.numFaces();
    int sz = pow2(lvl)+1;
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = wk[F];
        for(int i=-1; i<sz+1; i++)
            for(int j=-1; j<sz+1; j++)
                CD(i,j) = rot * CD(i,j);
    }
    return 0;
}
// ---------------------------------------------------------------------- 
int CCSurfOp::shift(GpMesh& gm, int lvl, DataLevel& wk, Point3 s)
{
  
    int numf = gm.numFaces();
    int sz = pow2(lvl)+1;
    for(int F=0; F<numf; F++) {
        CCRect<Point3>& CD = wk[F];
        for(int i=-1; i<sz+1; i++)
            for(int j=-1; j<sz+1; j++)
                CD(i,j) = s + CD(i,j);
    }
    return 0;
}


// ---------------------------------------------------------------------- 
int CCSurfOp::eno2index(int lvl, int e, int t, int& i, int& j)
{ //given a level and an edge and a position within the edge, return the indeces
  //in the rectangle data structure that points to that position on the edge
  
    int size = pow2(lvl)+1;
  
    switch(e) {
        case 0: i = t; j = 0; break;
        case 1: i = size-1; j = t; break;
        case 2: i = size-1-t; j = size-1; break;
        case 3: i = 0; j = size-1-t; break;
        default: assert(0);
    }
    return 0;
}

// ---------------------------------------------------------------------- 
int CCSurfOp::getcenter(GpMesh& gm, int lvl, DataLevel& wk, int V, Point3& val)
{ //get the center vertex in a data level

    
    intpair Fv = gm.Vf2Fv(V,0);
    int i, j;
    eno2index(lvl, Fv.second, 0, i, j);
    val = wk[Fv.first](i,j);  
    return 0;
}

// ---------------------------------------------------------------------- 

int CCSurfOp::putcenter(GpMesh& gm, int lvl, DataLevel& wk, int V, Point3& val)
{ //put center back
  
  
    int K = gm.valence(V);
    for(int f=0; f<K; f++) {
        intpair Fv = gm.Vf2Fv(V,f);
        int i, j;
        eno2index(lvl, Fv.second, 0, i, j);
        wk[Fv.first](i,j) = val;
    }
  
    return 0;
}
 
// ---------------------------------------------------------------------- 
int CCSurfOp::getonering(GpMesh& gm, int lvl, DataLevel& wk, int V, vector<Point3>& val)
{ //get the one ring around vertex V, store them in vector val. 
  //the ordering here should match that in orgcc.m.
  
  
    int K = gm.valence(V);
    val.resize(2*K+1);
  
    for(int f=0; f<K; f++) {
        intpair Fv = gm.Vf2Fv(V,f);
        int i,j;
        eno2index(lvl, Fv.second, 0, i, j);
        int di,dj;
        eno2index(lvl, Fv.second, 1, di, dj);
        val[0] = wk[Fv.first](i,j);
        val[1  +f] = wk[Fv.first](di,dj);
        val[1+K+f] = wk[Fv.first](di-(dj-j),dj+(di-i));
    }
  
    return 0;
}


// ---------------------------------------------------------------------- 
int CCSurfOp::putonering(GpMesh& gm, int lvl, DataLevel& wk, int V, vector<Point3>& val)
{ //do opposite of getonering to restore new values stored in val into the datalevel.
  
  
    int K = gm.valence(V);
    assert(val.size()==2*K+1);
    for(int f=0; f<K; f++) {
        intpair Fv = gm.Vf2Fv(V,f);
        int i,j;
        eno2index(lvl, Fv.second, 0, i, j);
        int di,dj;
        eno2index(lvl, Fv.second, 1, di, dj);
        wk[Fv.first](i,j) = val[0];
        wk[Fv.first](di,dj) = val[1  +f];
        wk[Fv.first](di-(dj-j),dj+(di-i)) = val[1+K+f];
    }
  
    return 0;
}


// ---------------------------------------------------------------------- 
int CCSurfOp::gethalfring(GpMesh& gm, int lvl, DataLevel& wk, int V, vector<Point3>& val)
{ //get the half ring around a boundary vertex store it in vector val.
  //the ordering here should match that in bouncc.m 

  
    int K = gm.valence(V);
    val.resize(2*K+2);

    for(int f=0; f<K; f++) {
        intpair Fv = gm.Vf2Fv(V,f);
        int i,j;    eno2index(lvl, Fv.second, 0, i, j);
        int di,dj;  eno2index(lvl, Fv.second, 1, di, dj);
        int ddi = di-(dj-j);   int ddj = dj+(di-i);
        int dddi = i-(dj-j);   int dddj = j+(di-i);
        val[0] = wk[Fv.first](i,j);
        val[f+1] = wk[Fv.first](di, dj);
        val[K+f+2] = wk[Fv.first](ddi, ddj);   
        if (f==K-1)
            val[K+1] = wk[Fv.first](dddi,dddj);
    }
  
    return 0;
}
// ---------------------------------------------------------------------- 

int CCSurfOp::puthalfring(GpMesh& gm, int lvl, DataLevel& wk, int V, vector<Point3>& val)
{ //do opposite to store the half ring in vector val into the data level
    

    int K = gm.valence(V);
    assert(val.size()==2*K+2);
  
    //note: in get half ring we only take 2*K+2 nodes around the ring. however,
    //after we're done with processing we want to re-store 3k+1 points back, to 
    //make sure the continuity among type "1" edges (in the interior not on the boundary)
    //is satisfied. this is not a problem with regular catmull clark, since the new added edge point
    //is always right in the middle of previos one.in HD(Henning-Denis) method this is not true.
    //so restore for edge midpoint, the two values coming from the two faces on both sides of it.
    for(int f=0; f<K; f++) {
        intpair Fv = gm.Vf2Fv(V,f);
        int i,j;    eno2index(lvl, Fv.second, 0, i, j);
        int di,dj;  eno2index(lvl, Fv.second, 1, di, dj);
        int ddi = di-(dj-j);   int ddj = dj+(di-i);

        int dddi = i-(dj-j);   int dddj = j+(di-i);
        wk[Fv.first](i,j) = val[0];
        wk[Fv.first](di, dj) = val[f+1];
        wk[Fv.first](ddi, ddj) = val[K+f+2]; 
    
        if (f==K-1)
            wk[Fv.first](dddi,dddj) = val[K+1];
        else wk[Fv.first](dddi, dddj) = val[f+2];
    }

    return 0;
  
}

