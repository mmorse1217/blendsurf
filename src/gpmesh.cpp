//general purpose quad mesh.

#include "gpmesh.hpp"

// ---------------------------------------------------------------------- 
int GpMesh::setup(istream& fin, int flipNormal, Point3 scale)
{
    
    assert(fin.good());
  
    //0. read in point and coordindex
    char tmp[100];
    fin>>tmp; cout << "tmp: " << tmp << endl;assert(strcmp(tmp,"points")==0);
    int nv;
    fin>>nv;
    cout << "num verts: " << nv << endl;
    vector<Point3> points;  points.resize(nv);
    for(int i=0; i<nv; i++) {
        double x,y,z;
        fin>>x>>y>>z;
        points[i] =  Point3(x,y,z);
        cout << i << ": " << points[i] << endl;
    }
    
    fin>>tmp;
    cout << "tmp: " << tmp << endl;
    assert(strcmp(tmp,"coords")==0);
    int nf;
    fin>>nf;
    vector< vector<int> > coords;  coords.resize(nf);
    for(int i=0; i<nf; i++) {
        coords[i].resize(4);
        for(int k=0; k<4; k++) {
            if(flipNormal==0)
                fin>>coords[i][k];
            else
                fin>>coords[i][3-k];
                cout << coords[i][k] << ", ";
        }
        cout << endl;
    } 
    int cnt;
    fin>>tmp; 
    cout << "tmp: " << tmp << endl;
    assert(strcmp(tmp,"gids") == 0);
    fin>>cnt;
    _face_group_ids.resize(cnt);
    //vector<int> face_group_ids;
    for(int i=0; i<cnt; i++) {
       int g; fin>>g;
       _face_group_ids[i] = g;
    }
    assert(_face_group_ids.size()==coords.size()); //both==number of faces

    // built _gids(); need a map from vertices -> group ids
    // face_group_ids maps faces -> group ids
    
    // Iterate over coordinates and assign groups to each vertex based on 
    // group ids of faces containing them. We assume here that options files are
    // written coherently, so all four faces containing a vertex will all be in
    // the same group.

    // One group id per input vertex
    // THIS IS NOT THE LIST OF gids IN THE WRL FILE
    // LOADS THOSE VALUES INTO face_group_ids AND DISTRIBUTES PER VERTEX
    _gids.resize(nv);
    int vertex;
    for(int i = 0; i < nf; i++){ // for each face...
        for(int j = 0; j < 4; j++){ // for each vertex in the ith face...
            vertex = coords[i][j];
            _gids[vertex] = _face_group_ids[i];
        }
    }


    fin>>tmp; assert(strcmp(tmp,"groups") == 0);
    fin>>cnt; 
    _interior_points.resize(cnt);
    _intpt_orientation.resize(cnt);
    double x,y,z; 
    int orientation;
    for(int i=0; i<cnt; i++) {
       fin>>tmp;
       assert(strcmp(tmp,"intpt") == 0);
       fin>>x>>y>>z; 
       _interior_points[i] = Point3(x,y,z); //interior point

       fin>>tmp;
       assert(strcmp(tmp,"orient") == 0);
       fin>>orientation; //orientation
       _intpt_orientation[i] = orientation;
    }
    
    
    _corner.resize(nv);

    fin>>tmp;
    int nr, temp;
    if(fin.eof()==false && strcmp(tmp, "convex") == 0) { //not end yet//there are convex corners
        fin>>nr;
	//        cout<<"nr convex is "<<nr<<endl;
        for(int i=0; i<nr; i++){
            fin>>temp;
	    _corner[temp] = CONVEX_VERTEX;
        }
    }
    fin>>tmp;
    if (fin.eof() == false && strcmp(tmp,"concave")==0){ //there are concave corners
        fin>>nr;
        //cout<<"nr concave is "<<nr<<endl;
        for(int i=0; i<nr; i++){
            fin>>temp;
	    _corner[temp] = CONCAVE_VERTEX;
        } 
    }
  
  
    //store points in the mesh with rescaling if necessary
    //1. _vpoint
    _vpoint.resize(nv);
    for(int V=0; V<nv; V++)
        _vpoint[V] = ewmul(points[V], scale); //rescaling
  
    //the following few data structures keep track of incidence relations.
  
    //2. Fe2Fe:Given a face and an edge, return the face on the other side and the edge number
    //that corresponds to the given edge in the neighboring face.
    //also takes care of _boun which stores a 1 for boundary vertices and 0 for interior 
    //vertices
 
    _Fe2Fe.resize(nf);
    for(int F=0; F<nf; F++) _Fe2Fe[F].resize(4);  //each face has 4 edges
  
    //given two vertices(in order) find the face and edge that corresponds to the pair.
    //store that in vv2fe
    map<intpair, intpair> vv2fe;
    for(int F=0; F<nf; F++)
        for(int e=0; e<4; e++) {
            int fm = coords[F][fromvno(e)]; //starting vertex
            int to = coords[F][gotovno(e)];
            // assert(vv2fe.find(intpair(fm,to))==vv2fe.end()); //cannot find 
            vv2fe[intpair(fm,to)] = intpair(F,e);
        }

    _boun.resize(nv);

    //no find in the map the opposite edge (edges go in opposite directions in 
    //neighboring faces)
    for(int F=0; F<nf; F++)
        for(int e=0; e<4; e++) {
            int fm = coords[F][fromvno(e)];
            int to = coords[F][gotovno(e)];      
            map<intpair, intpair>::iterator mi = vv2fe.find(intpair(to,fm));  
            if (mi==vv2fe.end()) { //if not found, then it's a boundary edge
                _boun[fm] = BOUNDARY_VERTEX;
                _boun[to] = BOUNDARY_VERTEX;
                //if boundary
                intpair a = intpair(-1, -1);//denote that with pair <-1, -1>
                _Fe2Fe[F][e] = a;
            }
            else // found it!- interior
                _Fe2Fe[F][e] = ((*mi).second);
        }
    vv2fe.clear();
  
  
    //3. Vf2Fv
    //Given a vertex and a face number, returns the global index of te face and the
    //vertex number of that face that coresponds to the center vertex V.
    //keeps track of all faces around a vertex in ccw order.

    _Vf2Fv.resize(nv);
  
    for(int F=0; F<nf; F++)
        for(int e=0; e<4; e++) {
            int V = coords[F][fromvno(e)]; //starting vertex
            if(_Vf2Fv[V].size()==0) { //have not touched yet
                vector<intpair> Fes;
                intpair cFe = intpair(F,e);
                int goback = 0;
                do {
                    if (cFe == intpair(-1,-1)) { //if we hit a boundary need to go back
                        //from the starting point to pick other faces
                        //to the boundary on other side.
                        goback = 1;
                        break;
                    }
                    //otherwise, keep storing the (face,edge) pairs in Fes
                    Fes.push_back(cFe);
                    cFe = _Fe2Fe[cFe.first][preveno(cFe.second)];
                } while (cFe != intpair(F,e));
	
                if (goback == 1){ //need to look in other direction
                    vector<intpair> Fes2;	  	  
                    intpair cFe = intpair(F,e);
                    cFe = _Fe2Fe[cFe.first][cFe.second];
                    //now keep track of a new set of "face, edge" pairs that lie on the other
                    //side of our starting point - until we find the boundary
                    if (cFe.first!=-1){
                        do{
                            cFe.second = nexteno(cFe.second);
                            Fes2.push_back(cFe);
                            cFe.second = preveno(cFe.second);
                            cFe = _Fe2Fe[cFe.first][nexteno(cFe.second)];
	      
                        }while(cFe != intpair(-1, -1));
                    }
	  
                    //now need to reverse the contents of Fes2 and append Fes1 to it, to follow
                    //the same direction around a vertex as in the interior case.
                    vector<intpair> Fes3;
                    for (int k=Fes2.size()-1; k>=0; k--)
                        Fes3.push_back(Fes2[k]);
	  
                    for (int k=0; k<Fes.size(); k++)
                        Fes3.push_back(Fes[k]);
	  
                    Fes = Fes3;
                }
	
		//                assert(Fes.size() >= 2);  //a vertex should at least have valance 2.
	
	
                int idx = 0;
                if (_boun[V] == INTERIOR_VERTEX) { //interior vertex
                    //need some kind of a start point to keep things uniform in the interior
                    //order them starting with the one with the least face numbers
	  
                    int mf = Fes[0].first;
                    int idx = 0;
                    for(int k=0; k<Fes.size(); k++) 
                        if(Fes[k].first < mf) {
                            mf = Fes[k].first;
                            idx = k;
                        }
                }
                //for boundary keep the original ordering (idx = 0)

                //store Fes in the original data structure.
                _Vf2Fv[V].insert(_Vf2Fv[V].end(), Fes.begin()+idx, Fes.end());
                _Vf2Fv[V].insert(_Vf2Fv[V].end(), Fes.begin(), Fes.begin()+idx);
                assert(_Vf2Fv[V].size()==Fes.size());
            }
        }
  
  
    //4. Fv2Vf
    //this is the exact opposite of Vf2fv
    //given a face and a vertex number find which face/vertex pair it corresponds to.
    _Fv2Vf.resize(nf);
    for(int F=0; F<nf; F++)
        _Fv2Vf[F].resize(4);
  
    for(int V=0; V<nv; V++)
        for(int f=0; f<_Vf2Fv[V].size(); f++) {
            intpair Fv = _Vf2Fv[V][f];
            _Fv2Vf[Fv.first][Fv.second] = intpair(V,f);
        }

    for(int i=0; i<nv; i++) { 
       cerr<< i << " " << valence(i)<<endl;
      if(  valence(i)  == 1 ) 	    _corner[i] = CONVEX_VERTEX;
    }
   cerr <<"gpmesh setup done"<<endl;
    /*
    ofstream out("debugout");
    dump(out);
    */
    return 0;
}



// ---------------------------------------------------------------------- 
int GpMesh::dump(ostream &fout)
{//output the contents of gpmesh.
  
  
    assert(fout.good());
  
    fout<<"vpoint"<<endl;
    for(int i=0; i<_vpoint.size(); i++)
        fout<<_vpoint[i]<<endl;
    fout<<endl;
  
    fout<<"Fe2Fe"<<endl;
    for(int F=0; F<numFaces(); F++) {
        for(int e=0; e<4; e++) {
            intpair tmp = _Fe2Fe[F][e];
            fout<<tmp.first<<" "<<tmp.second<<" | ";
        }
        fout<<endl;
    }
    fout<<endl;
  
    fout<<"Vf2Fv"<<endl;
    for(int V=0; V<numVertices(); V++) {
        for(int f=0; f<valence(V); f++) {
            intpair tmp = _Vf2Fv[V][f];
            fout<<tmp.first<<" "<<tmp.second<<" | ";
        }
        fout<<endl;
    }
  
    fout<<"Fv2Vf"<<endl;
    for(int F=0; F<numFaces(); F++) {
        for(int v=0; v<4; v++) {
            intpair tmp = _Fv2Vf[F][v];
            fout<<tmp.first<<" "<<tmp.second<<" | ";
        }
        fout<<endl;
    }
  
    return 0;
}


// ---------------------------------------------------------------------- 
Point3 GpMesh::ctr()
{//return center of mesh (average of all points in mesh)
    Point3 sum( 0.0);
    for(int V=0; V<_vpoint.size(); V++)
        sum += _vpoint[V];
    return sum/double(_vpoint.size());
}

// ---------------------------------------------------------------------- 
void GpMesh::bbox(Point3& bbmin, Point3& bbmax)
{//return the bounding box of a mesh.
    bbmin = Point3( SCL_MAX);
    bbmax = Point3(-SCL_MAX);
    for(int V=0; V<_vpoint.size(); V++) {
        bbmin = min(bbmin, _vpoint[V]);
        bbmax = max(bbmax, _vpoint[V]);
    }
}
