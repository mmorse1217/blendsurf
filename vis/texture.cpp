#include "texture.hpp"

Texture::Texture(BdSurf* surf, int lvl, char* texFile):_bdsurf(surf), _lvl(lvl)
{strcpy(_texFile,texFile);}
Texture::~Texture()
{}

void Texture::initData(vector< vector< Matrix< Point2 > > >& xy,
                       vector< vector< Matrix< Point3 > > >& pos,
                       vector< vector< Matrix< Point3 > > >& nor){
    _Vfxy = xy;
    _Vfpos = pos;
    _Vfnor = nor;
    _texturectrl = TEX_TYPE_IMAGE;
}

void Texture::updateTextureParam(){
    updateOrientation(_oldTextFaceNumber, _textureFaceNumber);
    getSurfOrientation();
    // cout << "center : " << _textureCenter << "\norient : "<<_textureOrientation<<cout<<"\nface : "<<_textureFaceNumber<<"\n";
    calculateTexture(_textureCenter, _textureOrientation, _textureFaceNumber);
}

GLubyte* Texture::read(int* width, int* height){
    _texImage =  readTexture(_texFile, width, height);
    _texSize[0] = *width;
    _texSize[1] = *height;
    return _texImage;
}

void Texture::storeOldPos(){
    _oldTextureNormal = _textureSurfNormal;
    _oldTextureSurfCenter = _textureSurfCenter;
    _oldTextureCenter = _textureCenter;
}

void Texture::storePos(int face, Point2& textCenter, Point3& textSurfCenter, Point3& textSurfNormal){
   _textureFaceNumber = face;
   _textureCenter = textCenter;
   _textureSurfCenter = textSurfCenter;
   double len = sqrt(dot(textSurfNormal,textSurfNormal));
   _textureSurfNormal = textSurfNormal / len;
}

void Texture::switchTextureMode(){
    _texturectrl = (_texturectrl+1) % TEX_TYPE_TTL;
    delete _texImage;
    _texImage = NULL;
}


void Texture::rotateTexture(int direction){
    double pi = acos(-1.);
    if (direction == ROTATE_CW){
        //    rotateTextureVector(pi/30);
        _newTextureRotation += -pi/30.;
    } else if (direction == ROTATE_CCW){
        //    rotateTextureVector(-pi/30);
        _newTextureRotation += pi/30.;
    }
}

void Texture::moveTexture(int direction){

    double d[2] = {0,0};
    double a = 0.05;
    switch(direction){
        case MOVE_LEFT:
            d[0] -= a;
            break;
        case MOVE_RIGHT:
            d[0] += a;
            break;
        case MOVE_UP:
            d[1] += a;
            break;
        case MOVE_DOWN:
            d[1] -= a;
            break;
    }

    moveTextureCenter(d);
    _moving = true;
}

void Texture::remove(){
    _textureFaces.clear();
    delete _texImage;
    _texImage = NULL;
}

void Texture::setTextureMode(int faceNum, Point2& textCenter, Point3& textSurfCenter, 
                        Point3& textSurfNormal){
    // default texture orientation
    _textureSurfOrient = Point3(0,0,0);
    _textureSurfNormal = Point3(0,0,0);
    _textureSurfNormal = textSurfNormal;
    _textureSurfCenter = textSurfCenter;
    _textureCenter = textCenter;
    _textureFaceNumber = faceNum;
    _oldTextFaceNumber = _textureFaceNumber;
    _oldTextureNormal = _textureSurfNormal;
    _oldTextureSurfCenter = _textureSurfCenter;
    _oldTextureCenter = _textureCenter;
    // make necessary initializations for texture
    _textureOrientation[0] = 1;  _textureOrientation[1] = 0;
    _centerVertex = 0;
    _dispMatrix[0] = -1;     _dispMatrix[1] = 0;
    _dispMatrix[2] = 0;      _dispMatrix[3] = -1;
    _orientMatrix[0] = 1;    _orientMatrix[1] = 0;
    _orientMatrix[2] = 0;    _orientMatrix[3] = 1;
    _textureScale = .4;
    _dispMapScale = .1;
    _oldTextureRotation = 0.;
    _newTextureRotation = 0.;
    _moving = false;
}

void Texture::getSurfOrientation(){

    int face = _textureFaceNumber;
    double cdCenter[2] = {_textureCenter[0], _textureCenter[1]};
    double cdVec[2] = {_textureOrientation[0], _textureOrientation[1]};

  
    int i = 0;
    bool done = false;
    bool inmove = _moving;

    /*
     * take the orientation in cd coordinates and transform to surface coordinates
     * do this just for once, when the texture is activated
     * other times we rely on the surface normal calculation below
     */
    while((!done && _textureSurfOrient == Point3(0,0,0)) || (!done && _moving)){
        intpair Vf = _bdsurf->gpmesh().Fv2Vf(face,i);
  
        int V = Vf.first;        //vertex number
        int f = Vf.second;       //face number

        double VcdCenter[6];
        _bdsurf->Fcd2Vfcd(BdSurf::EVAL_VALUE | BdSurf::EVAL_1ST_DERIV, face, cdCenter, i, V, f, VcdCenter);

        double VcdVec[2];
        VcdVec[0] = cdVec[0] * VcdCenter[2] + cdVec[1] * VcdCenter[4];
        VcdVec[1] = cdVec[0] * VcdCenter[3] + cdVec[1] * VcdCenter[5];

        if (VcdCenter[0] < (1-1./32.) && VcdCenter[1] < (1-1./32.)){

            double xyCenter[6];
            double xyVec[6];

            _bdsurf->Vfcd2Vxy(BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, V, f, VcdCenter, xyCenter);

            xyVec[0]  = xyCenter[2] * VcdVec[0]  + xyCenter[4] * VcdVec[1];
            xyVec[1]  = xyCenter[3] * VcdVec[0]  + xyCenter[5] * VcdVec[1];
   
            Point3 srfc[6];
            _bdsurf->eval(BdSurf::EVAL_VALUE | BdSurf::EVAL_1ST_DERIV, V,  xyCenter, srfc);

            NumMatrix J(3,2);
            J(0,0) = srfc[1][0];
            J(0,1) = srfc[2][0];
            J(1,0) = srfc[1][1];
            J(1,1) = srfc[2][1];
            J(2,0) = srfc[1][2];
            J(2,1) = srfc[2][2];
    
            // v at the surface
            Point3 srfVec;
            srfVec[0] = J(0,0) * xyVec[0] + J(0,1) * xyVec[1];
            srfVec[1] = J(1,0) * xyVec[0] + J(1,1) * xyVec[1];
            srfVec[2] = J(2,0) * xyVec[0] + J(2,1) * xyVec[1];

            Point3 srfNormal = cross(srfc[1], srfc[2]).dir();
            //      cout << _textureSurfNormal << endl << srfNormal << endl;
            if (_moving){
                _textureSurfNormal = srfNormal;
            }

            if (!_moving){
                double len = sqrt(dot(srfVec,srfVec));
                _textureSurfOrient = srfVec / len;
            }
            done = true;
            _moving = false;
        }
        i++;
    }

    //  cout << "1 --> "<< _textureSurfOrient << endl;
    //update the surface orientation according to the surface normal

    Point3& from = _oldTextureNormal;
    Point3& to   = _textureSurfNormal;

    double err = 1e-10;
    if( abs(from[0] - to[0])>err && abs(from[1] - to[1])>err && abs(from[2] - to[2])>err ){


        //     cout << "From " << from << endl;
        //     cout << "To   " << to << endl;

        double lenfrom = sqrt(dot(from,from));
        double lento = sqrt(dot(to,to));

        //     cout << "From " << lenfrom << endl;
        //     cout << "To   " << lento << endl;

        from /= lenfrom;
        to /= lento;

        Point3 v = cross(from, to).dir();       // rotation axis
        Point3 vs = cross(from,to);
        //cout << v << endl << vs << endl;


        double c = dot(from, to) ;              // cos alpha(rotation angle)
        double s = sqrt(1-c*c);                 // sin alpha
        s = sqrt(dot(vs,vs)) / (lenfrom*lento);

        //    double theta = atan(s/c);
        //cout << theta << endl;

        double u = 1 - c;

        NumMatrix R(3,3);

        double vx = v[0];    double vy = v[1];    double vz = v[2];
        double vsx = vs[0];    double vsy = vs[1];    double vsz = vs[2];

        R(0,0) = vx * vx * u + c;
        R(0,1) = vy * vx * u - vsz;
        R(0,2) = vz * vx * u + vsy;

        R(1,0) = vx * vy * u + vsz;
        R(1,1) = vy * vy * u + c;
        R(1,2) = vz * vy * u - vsx;

        R(2,0) = vx * vz * u - vsy;
        R(2,1) = vy * vz * u + vsx;
        R(2,2) = vz * vz * u + c;

        Point3& tso = _textureSurfOrient;
        Point3 t = _textureSurfOrient;
    
        tso[0] = R(0,0) * t[0] + R(0,1) * t[1] + R(0,2) * t[2];
        tso[1] = R(1,0) * t[0] + R(1,1) * t[1] + R(1,2) * t[2];
        tso[2] = R(2,0) * t[0] + R(2,1) * t[1] + R(2,2) * t[2];

    
        _oldTextureNormal = _textureSurfNormal;
    }
  
    Point3& tso = _textureSurfOrient;
    double len = sqrt(dot(tso, tso));
    tso /= len;
    //  cout << "2 --> "<< _textureSurfOrient << endl;

    // rotate the texture
    double textureRotation = _newTextureRotation - _oldTextureRotation;
    if(inmove){
        // textureRotation = _oldTextureRotation;
    }
    if(textureRotation != 0){
        Point3 t = _textureSurfOrient;
        rotateVector(t, _textureSurfNormal, textureRotation, _textureSurfOrient);
        textureRotation = 0.;
    }
    _oldTextureRotation = _newTextureRotation;
}

/*
 * a : vector to be rotated
 * v : rotation axis
 * ang : rotation angle
 * r : rotated vector
 */
void Texture::rotateVector(Point3& a, Point3& v, double ang, Point3& r){

    double s = sin(ang);
    double c = cos(ang);
    double u = 1 - c;
    Point3 vs = v*s;

    NumMatrix R(3,3);

    double vx = v[0];    double vy = v[1];    double vz = v[2];
    double vsx = vs[0];    double vsy = vs[1];    double vsz = vs[2];

    R(0,0) = vx * vx * u + c;
    R(0,1) = vy * vx * u - vsz;
    R(0,2) = vz * vx * u + vsy;

    R(1,0) = vx * vy * u + vsz;
    R(1,1) = vy * vy * u + c;
    R(1,2) = vz * vy * u - vsx;

    R(2,0) = vx * vz * u - vsy;
    R(2,1) = vy * vz * u + vsx;
    R(2,2) = vz * vz * u + c;
  
    r[0] = R(0,0) * a[0] + R(0,1) * a[1] + R(0,2) * a[2];
    r[1] = R(1,0) * a[0] + R(1,1) * a[1] + R(1,2) * a[2];
    r[2] = R(2,0) * a[0] + R(2,1) * a[1] + R(2,2) * a[2];
}

/*
 * calculates the rotation matrix for the texture orientation in cd coordinates
 * this matrix avoids flipping when switching from face to face.
 * this actually transformation from face coordinates to texture coordinates
 */
void Texture::updateOrientation(int of, int nf){
    int pe = -1;
    int ne = -1;
    double pi = acos(-1.);

    if(nf != of){
        //    cout << "face change --> " << of << " to " << nf << endl;

        bool found = false;
        for (int v=0; v<4; v++){
            int nV = _bdsurf->gpmesh().Fv2Vf(of,v).first;
            int val = _bdsurf->gpmesh().valence(nV);
            for(int f=0; f<val; f++){
                int nF = _bdsurf->gpmesh().Vf2Fv(nV, f).first;
                if(nF == nf){
                    pe = v;
                    ne = _bdsurf->gpmesh().Vf2Fv(nV, f).second;
                    found = true;
                    break;
                }
            }
            if(found) break;
        }

        //     for (int e=0; e<4; e++){
        //       int nF = _bdsurf->gpmesh().Fe2Fe(of, e).first;
        //       if(nF == nf){
        //   pe = e;
        //  ne = _bdsurf->gpmesh().Fe2Fe(of, e).second;
        //  break;
        //       }
        //     }

        // _moving toooooo fast
        if(found){

            _oldTextFaceNumber = _textureFaceNumber;
            assert( pe < 4 && pe>= 0);
            int k = pe - ne;
            double a = pi*k/2.;
      
            double dM[4];
            dM[0] = cos(a);
            dM[1] = sin(a);
            dM[2] = -sin(a);
            dM[3] = cos(a);
            //     double ndM[4];
            //     ndM[0] = _orientMatrix[0] * dM[0] + _orientMatrix[1] * dM[2];
            //     ndM[1] = _orientMatrix[0] * dM[1] + _orientMatrix[1] * dM[3];
            //     ndM[2] = _orientMatrix[2] * dM[0] + _orientMatrix[3] * dM[2];
            //     ndM[3] = _orientMatrix[2] * dM[1] + _orientMatrix[3] * dM[3];
      
            _orientMatrix[0] = dM[0];
            _orientMatrix[1] = dM[1];
            _orientMatrix[2] = dM[2];
            _orientMatrix[3] = dM[3];
      
            Point2 t;
            t[0] = _textureOrientation[0] * dM[0] + _textureOrientation[1] * dM[1];
            t[1] = _textureOrientation[0] * dM[2] + _textureOrientation[1] * dM[3];
      
            _textureOrientation = t;
      
            //_textureSurfOrient = Point3(0,0,0);
            //    cout << _oldTextureNormal <<endl << _textureSurfNormal << endl;
        }else{
            _textureSurfNormal = _oldTextureNormal;
            _textureSurfCenter = _oldTextureSurfCenter;
            _textureCenter = _oldTextureCenter;
        }
    }
}


/*
 * calculates the texture coordinates in the union of 4 charts sharing the face on which the texture 
 * center is located
 */
void Texture::calculateTexture(Point2 center, Point2 orient, int face){

    intpair Vf[4];

    vector<int>& faces = _textureFaces;  // faces on which texture is ON

    double cdCenter[2];       // x
    double xyCenter[4][6];    // x_i
    double cdVec[2];          // v
    //double cdVecP[2];         // v^p
    double xyVec[4][2];       // v_i
    double xyVecP[4][2];      // v_i^p
  
    cdCenter[0] = center[0];
    cdCenter[1] = center[1];

    cdVec[0] = orient[0];
    cdVec[1] = orient[1];
    //  cout << "Orientation (cd)--> " << cdVec[0] << "\t" << cdVec[1] << endl;
    // cout << "Center (cd)--> " << cdCenter[0] << "\t" << cdCenter[1] << endl;
  
    double consLen = _textureScale;
    //  cout << "Scale --> " << consLen << endl;

    //double pi = acos(-1.);
    int nsteps = pow2(_lvl);                 //number of intervals
    double step = 1 / double(nsteps);        //defined in a unit square, length of every interval
    bool chartUsed[4] = {false, false, false, false};

    Mat2T<double> invV[4];

    /*
     * calculate x_i, v_i and v_i^p for 4 charts
     */
    for (int i=0; i<4; i++){
        //cout <<endl << "Chart number --> " << i << endl;
        /*
         * switch from (Face / vertex) to (Vertex / face)
         */
        Vf[i] = _bdsurf->gpmesh().Fv2Vf(face,i);
        int V = Vf[i].first;        //vertex number
        int f = Vf[i].second;       //face number

        /*
         * calculate x_i
         *
         * 1st - calculate center point coordinates in Vfcd
         * 2nd - transform to xy plane
         */
        double VcdCenter[6];
        _bdsurf->Fcd2Vfcd(BdSurf::EVAL_VALUE | BdSurf::EVAL_1ST_DERIV, face, cdCenter, i, V, f, VcdCenter);

        double VcdVec[2];
        VcdVec[0] = cdVec[0] * VcdCenter[2] + cdVec[1] * VcdCenter[4];
        VcdVec[1] = cdVec[0] * VcdCenter[3] + cdVec[1] * VcdCenter[5];


        if (VcdCenter[0] < (1-1./32.) && VcdCenter[1] < (1-1./32.)){
            chartUsed[i] = true;


            _bdsurf->Vfcd2Vxy(BdSurf::EVAL_VALUE|BdSurf::EVAL_1ST_DERIV, V, f, VcdCenter, xyCenter[i]);

            //cout << i << "\t" << xyCenter[i][0] << "\t" << xyCenter[i][1] << endl;
            /*
             * calculate v_i
             */
            xyVec[i][0]  = xyCenter[i][2] * VcdVec[0]  + xyCenter[i][4] * VcdVec[1];
            xyVec[i][1]  = xyCenter[i][3] * VcdVec[0]  + xyCenter[i][5] * VcdVec[1];

            /*
             * calculate v_i^p
             */
            // map from xy to surface
            Point3 srfc[6];
            _bdsurf->eval(BdSurf::EVAL_VALUE | BdSurf::EVAL_1ST_DERIV, V,  xyCenter[i], srfc);

            NumMatrix J(3,2);
            J(0,0) = srfc[1][0];     J(0,1) = srfc[2][0];
            J(1,0) = srfc[1][1];     J(1,1) = srfc[2][1];
            J(2,0) = srfc[1][2];     J(2,1) = srfc[2][2];
    
            Point3 srfCenter;
            srfCenter = srfc[0];
    
            Point3 srfNormal = cross(srfc[1], srfc[2]).dir();   // normal direction at this point on surface
    
            // v at the surface
            Point3 srfVec;
            //       srfVec[0] = J(0,0) * xyVec[i][0] + J(0,1) * xyVec[i][1];
            //       srfVec[1] = J(1,0) * xyVec[i][0] + J(1,1) * xyVec[i][1];
            //       srfVec[2] = J(2,0) * xyVec[i][0] + J(2,1) * xyVec[i][1];
            //       double len0 = sqrt(srfVec[0]*srfVec[0] + srfVec[1]*srfVec[1] + srfVec[2]*srfVec[2]);
            //       srfVec[0] /= len0;
            //       srfVec[1] /= len0;
            //       srfVec[2] /= len0;

            //cout << "3 --> " << _textureSurfOrient << endl << "4 --> " << srfVec << endl;

            // we use the surface vector calculated by getSurfOrientation
            //       Point3 t;
            //       t = srfVec;
            //      _textureSurfOrient = t;
            //      getSurfOrientation();
            srfVec = _textureSurfOrient;
            //      _textureSurfOrient = t;
            // we need v to be unit length.
            double len = sqrt(srfVec[0]*srfVec[0] + srfVec[1]*srfVec[1] + srfVec[2]*srfVec[2]);
            srfVec[0] /= len;
            srfVec[1] /= len;
            srfVec[2] /= len;
            //      double l = sqrt(dot(srfNormal, srfNormal));
            //      cout << l << "\t" << len << "\t" << consLen << endl;

            // calculate v^p at the surface
            Point3 srfVecP;
            srfVecP = cross(srfNormal, srfVec).dir();
            //      cout << "v " << srfVec << endl << " vp " << srfVecP << endl;
            // cout << dot(srfVec,srfVecP) << "\n";
            // scale v and v^p
            srfVec[0] *= consLen;
            srfVec[1] *= consLen;
            srfVec[2] *= consLen;
    
            srfVecP[0] *= consLen;
            srfVecP[1] *= consLen;
            srfVecP[2] *= consLen;
    
            // pinv does not work
            NumMatrix piJ(2,3);               // pseudo inverso of Jacobian
            NumMatrix jTj(2,2);               // J transpose * J (J^T *J)
            NumMatrix ijTj(2,2);              // inverse of (J^T * J)
            NumMatrix jT(2,3);                // J transpose (J^T)
            tran(J, jT);                      // calculate J^T
            NumMatrixMult(jT, J, jTj);        // calculate jTj
            inv(jTj, ijTj);                   // calculate inv(J^T * J)
            NumMatrixMult(ijTj, jT, piJ);     // calculate inv(J^T * J) * J^T
            //     pinv(J, 1e-10, piJ);

            //calculate v^p from pinv(Jacobian)
            xyVecP[i][0] = piJ(0,0) * srfVecP[0] + piJ(0,1) * srfVecP[1] + piJ(0,2) * srfVecP[2] ;
            xyVecP[i][1] = piJ(1,0) * srfVecP[0] + piJ(1,1) * srfVecP[1] + piJ(1,2) * srfVecP[2] ;
    
            //calculate v from pinv(Jacobian)
            xyVec[i][0] = piJ(0,0) * srfVec[0] + piJ(0,1) * srfVec[1] + piJ(0,2) * srfVec[2] ;
            xyVec[i][1] = piJ(1,0) * srfVec[0] + piJ(1,1) * srfVec[1] + piJ(1,2) * srfVec[2] ;
    
            /*
             * define matrix for the coordinate axes to use in si, ti solution
             */
            Mat2T<double> vMat;
            vMat(0,0) = xyVec[i][0];
            vMat(1,0) = xyVec[i][1];
            vMat(0,1) = xyVecP[i][0];
            vMat(1,1) = xyVecP[i][1];
            /*
             * inverse the matrix for linear system solution
             */
            invV[i] = vMat.inverse();
        }
    }

    /*
     * calculate list of faces in the union of 4 charts defined by 4 vertices at the corners of face F
     */
    for (int i=0; i<4; i++){                 // loop in 4 vertices(charts) at corners of F

        if(chartUsed[i]){

            int V = Vf[i].first;                   // vertex number
            int K = _bdsurf->valence(V);           // valence of V
            intpair Fv;

            for (int f=0; f<K; f++){               // loop in valences of V

                Fv = _bdsurf->gpmesh().Vf2Fv(V,f);   // Get global face number for this valence
                int F = Fv.first;                    // global face number
                /*
                 * if F is not in the list, then add
                 */
                if(find(faces.begin(),faces.end(),F) == faces.end()){
                    faces.push_back(F);                // add this global face number to the list of faces to be processed 
                }
            }
        }
    }

    int nface = faces.size(); // number of faces with texture
    _Vftxy.resize(nface);
    _Vftpos.resize(nface);
    _Vftnor.resize(nface);
    _Vfdpos.resize(nface);
    _Vfdnor.resize(nface);

    for (int tf=0; tf<nface; tf++){           // loop in all faces in the union of 4 charts sharing F

        Matrix<Point2>& cmxy  = _Vftxy[tf];     // y_i
        Matrix<Point3>& cmpos = _Vftpos[tf];    // point coordinates on surface
        Matrix<Point3>& cmnor = _Vftnor[tf];    // normal direction on surface
        Matrix<Point3>& cmdpos = _Vfdpos[tf];    // point coordinates on surface
        Matrix<Point3>& cmdnor = _Vfdnor[tf];    // normal direction on surface

        vector<intpair> VfLocal;                // some or all of 4 vertices that are at the center of 4 charts sharing F
        vector<int> chartNumber;
        intpair VfAll[4];                       // all 4 vertices at corners of this face

        for (int i=0; i<4; i++){                // loop over 4 corners of face tf

            intpair V = _bdsurf->gpmesh().Fv2Vf(faces[tf],i);
            VfAll[i] = V;                         // store global vertex number for corner vertices

            for (int j=0; j<4; j++){              // loop over 4 vertices (charts) sharing F

                if(V.first == Vf[j].first && chartUsed[j]){
                    VfLocal.push_back(V);             // store global vertex number if this vertex is 1 of 4 vertices of face F
                    chartNumber.push_back(j);         // store chart number for this vertex
                }
            }
        }
    
        int ncharts = VfLocal.size();    // number of charts to be used for this face
        //    cout << tf << "\t" << ncharts << endl;

        /*
         * calculate s_i and t_i
         */
        cmxy.resize(nsteps+1, nsteps+1);
        cmpos.resize(nsteps+1, nsteps+1);
        cmnor.resize(nsteps+1, nsteps+1);
        cmdpos.resize(nsteps+1, nsteps+1);
        cmdnor.resize(nsteps+1, nsteps+1);
        /*
         * loop over every internal point to calculate y_i, s_i and t_i.
         */
        int nmax = nsteps/2;
        for (int ii=0; ii<=nsteps; ii++){      // loop in the unit square by 1/32 (1/2^5) length steps
            for(int jj=0; jj<=nsteps; jj++){

                double si[4] = {0,0,0,0};          // using 4 is safe and easy for initialization
                double ti[4] = {0,0,0,0};          // instead of using ncharts
                double wi[4] = {0,0,0,0};

                for (int i=0; i<ncharts; i++){     // loop over some or all of 4 charts that share this face

                    int ci = chartNumber[i];

                    int V = (VfLocal[i]).first;      // global vertex number (also a vertex of face F)
                    int f = (VfLocal[i]).second;     // local face number for vertex V
                    double cd[2];                    // coordinates in unit square
                    double yi[2];                    // coordinates in mapped plane
      
                    //cout << V << "\t" << f << "\t" << ci << endl;

                    /*
                     * calculate y_i
                     */
                    cd[0] = ii * step;    // cd is in Fcd coordinates, relatice to 0th vertex of any face F
                    cd[1] = jj * step;
                    double Vcd[2];
                    int ff = faces[tf];
                    int fv = _bdsurf->gpmesh().Vf2Fv(V,f).second;
                    _bdsurf->Fcd2Vfcd(BdSurf::EVAL_VALUE, ff, cd, fv, V, f, Vcd);
                    _bdsurf->Vfcd2Vxy(BdSurf::EVAL_VALUE, V, f, Vcd, yi);  // do mapping from unit square to xy plane
                    yi[0] -= xyCenter[ci][0];         // tnransform from corner to center x_i
                    yi[1] -= xyCenter[ci][1];

                    /*
                     * calculate surface position and normal for y_i
                     */
                    if(i == 0){   // just to do once
                        //        int vv = _bdsurf->gpmesh().Vf2Fv(V,f).second;
                        int v0 = 0;
                        int v1 = 1;
                        int v2 = 2;
                        int v3 = 3;

                        if(ii <= nmax && jj >= nmax){
                            cmpos(ii,jj) = _Vfpos[VfAll[v3].first][VfAll[v3].second](nmax - (jj-nmax),ii);
                            cmnor(ii,jj) = _Vfnor[VfAll[v3].first][VfAll[v3].second](nmax - (jj-nmax),ii);
                        }else if( ii >= nmax && jj >= nmax){
                            cmpos(ii,jj) = _Vfpos[VfAll[v2].first][VfAll[v2].second](nmax - (ii-nmax), nmax - (jj-nmax));
                            cmnor(ii,jj) = _Vfnor[VfAll[v2].first][VfAll[v2].second](nmax - (ii-nmax), nmax - (jj-nmax));
                        }else if(ii >= nmax && jj <= nmax){
                            cmpos(ii,jj) = _Vfpos[VfAll[v1].first][VfAll[v1].second](jj,nmax - (ii - nmax));
                            cmnor(ii,jj) = _Vfnor[VfAll[v1].first][VfAll[v1].second](jj,nmax - (ii - nmax));
                        }else{
                            cmpos(ii,jj) = _Vfpos[VfAll[v0].first][VfAll[v0].second](ii, jj);
                            cmnor(ii,jj) = _Vfnor[VfAll[v0].first][VfAll[v0].second](ii, jj);
                        }
                        cmnor(ii,jj) = cmnor(ii,jj).dir();
                    }

                    /*
                     * calculate s_i, t_i
                     * contribution from ith chart
                     */
                    _bdsurf->blendFuncEval(BdSurf::EVAL_VALUE, /*V, f,*/ Vcd, 1/32., 1-1./32., &wi[i]);
                    si[i] = invV[ci](0,0) * yi[0] + invV[ci](0,1) * yi[1];
                    ti[i] = invV[ci](1,0) * yi[0] + invV[ci](1,1) * yi[1];

                }
                /*
                 * average s_i and t_i.
                 */
                double sit = weight(wi, si);
                double tit = weight(wi, ti);
                cmxy(ii,jj) = Point2(sit+0.5,tit+0.5);             // store s_i and t_i
                //-------------------------------------------------------------------
                //     cout << si[0]<<","<<si[1]<<","<<si[2]<<","<<si[3]<<"," << ti[0]<<","<<ti[1]<<","<<ti[2]<<","<<ti[3] << endl;
                //cout << sit << "\t" << tit << endl;
            }
        }
    }
    if(_texturectrl == TEX_TYPE_DISP){
        applyDisplacementMap();
    }
}


void Texture::applyDisplacementMap(){
    int w = _texSize[0];
    int h = _texSize[1];
    int nfaces = _textureFaces.size();
    for (int f=0; f<nfaces; f++){
        Matrix<Point2>& cmxy  = _Vftxy[f];
        Matrix<Point3>& cmpos = _Vftpos[f];
        Matrix<Point3>& cmnor = _Vftnor[f];
        Matrix<Point3>& cmdpos = _Vfdpos[f];
        Matrix<Point3>& cmdnor = _Vfdnor[f];
        int mmax = cmpos.m();
        int nmax = cmpos.n();
        // update location in normal direction
        for(int j=0; j<nmax; j++) {
            for(int i=0; i<mmax; i++) {
                double x = (cmxy(i,j)[0]);
                double y = (cmxy(i,j)[1]);
                if(x<=1 && x>=0 && y<=1 && y>=0){
                    int icoor = (int)(x*w);
                    int jcoor = (int)(y*h);
                    int index = jcoor*w+icoor;
                    //    cerr<< index <<endl;
                    cmdpos(i,j) = cmpos(i,j) + cmnor(i,j)*(_dispMapScale*_dispMap[index]);
                }else{
                    cmdpos(i,j) = cmpos(i,j);
                }
            }
        }
        // update normals using new locations
        int im,ip,jm,jp;
        for(int j=0; j<nmax; j++) {
            for(int i=0; i<mmax; i++) {
                double x = (cmxy(i,j)[0]);
                double y = (cmxy(i,j)[1]);
                if(x<=1 && x>=0 && y<=1 && y>=0){
                    im = i-1;   ip = i+1;   jm = j-1;   jp = j+1;
                    if(i == 0) im=i;    if(j == 0) jm=j;
                    if(i == mmax-1) ip=i;    if(j == nmax-1) jp=j;
                    Point3 n1 = cross((cmdpos(ip,j) - cmdpos(i,j)),(cmdpos(i,jp)-cmdpos(i,j))).dir();
                    Point3 n2 = cross((cmdpos(i,jm) - cmdpos(i,j)),(cmdpos(ip,j)-cmdpos(i,j))).dir();
                    Point3 n3 = cross((cmdpos(im,j) - cmdpos(i,j)),(cmdpos(i,jm)-cmdpos(i,j))).dir();
                    Point3 n4 = cross((cmdpos(i,jp) - cmdpos(i,j)),(cmdpos(im,j)-cmdpos(i,j))).dir();
                    cmdnor(i,j) = ((n1+n2+n3+n4)/4.).dir();
                }else{
                    cmdnor(i,j) = cmnor(i,j);
                }
            }
        }
    }
    // update normals at the boundaries
    for (int f=0; f<nfaces; f++){
        Matrix<Point2>& cmxy  = _Vftxy[f];
        Matrix<Point3>& cmpos = _Vftpos[f];
        Matrix<Point3>& cmdpos = _Vfdpos[f];
        Matrix<Point3>& cmdnor = _Vfdnor[f];
        Matrix<Point3>& cmnor = _Vftnor[f];

        vector<Point3> _newPos;
        _newPos.resize(5);

        vector<Point3>& cmdposnew = _newPos;
        int mmax = cmpos.m();
        int nmax = cmpos.n();

        //1st edge
        int j = 0;
        for (int i=1; i<mmax-1; i++){
            double x = (cmxy(i,j)[0]);
            double y = (cmxy(i,j)[1]);
            if(x<=1 && x>=0 && y<=1 && y>=0){
                Point3 n1 = cross((cmdpos(i+1,j) - cmdpos(i,j)),(cmdpos(i,j+1)-cmdpos(i,j))).dir();
                Point3 n4 = cross((cmdpos(i,j+1) - cmdpos(i,j)),(cmdpos(i-1,j)-cmdpos(i,j))).dir();

                if (calcBoundryPos(i,j,f,0, cmdposnew)){
                    Point3 n2 = cross((cmdposnew[3] - cmdposnew[0]),(cmdposnew[2] - cmdposnew[0])).dir();
                    Point3 n3 = cross((cmdposnew[1] - cmdposnew[0]),(cmdposnew[3] - cmdposnew[0])).dir();
                    cmdnor(i,j) = ((n1+n4+n2+n3)/4.).dir();
                }else{
                    cmdnor(i,j) = ((n1+n4)/2.).dir();
                }
            }else{
                cmdnor(i,j) = cmnor(i,j);
            }
        }

        //4th edge
        int i = 0;
        for(int j=1; j<nmax-1; j++) {
            double x = (cmxy(i,j)[0]);
            double y = (cmxy(i,j)[1]);
            if(x<=1 && x>=0 && y<=1 && y>=0){
                Point3 n1 = cross((cmdpos(i+1,j) - cmdpos(i,j)),(cmdpos(i,j+1)-cmdpos(i,j))).dir();
                Point3 n2 = cross((cmdpos(i,j-1) - cmdpos(i,j)),(cmdpos(i+1,j)-cmdpos(i,j))).dir();

                if (calcBoundryPos(i,j,f,3, cmdposnew)){
                    Point3 n3 = cross((cmdposnew[1] - cmdposnew[0]),(cmdposnew[3] - cmdposnew[0])).dir();
                    Point3 n4 = cross((cmdposnew[4] - cmdposnew[0]),(cmdposnew[1] - cmdposnew[0])).dir();
                    cmdnor(i,j) = ((n1+n4+n2+n3)/4.).dir();
                }else{
                    cmdnor(i,j) = ((n1+n2)/2.).dir();
                }
            }else{
                cmdnor(i,j) = cmnor(i,j);
            }
        }

        //2nd edge
        i = mmax-1;
        for(int j=1; j<nmax-1; j++) {
            double x = (cmxy(i,j)[0]);
            double y = (cmxy(i,j)[1]);
            if(x<=1 && x>=0 && y<=1 && y>=0){
                Point3 n3 = cross((cmdpos(i-1,j) - cmdpos(i,j)),(cmdpos(i,j-1)-cmdpos(i,j))).dir();
                Point3 n4 = cross((cmdpos(i,j+1) - cmdpos(i,j)),(cmdpos(i-1,j)-cmdpos(i,j))).dir();

                if (calcBoundryPos(i,j,f,1, cmdposnew)){
                    Point3 n1 = cross((cmdposnew[2] - cmdposnew[0]),(cmdposnew[4] - cmdposnew[0])).dir();
                    Point3 n2 = cross((cmdposnew[3] - cmdposnew[0]),(cmdposnew[2] - cmdposnew[0])).dir();
                    cmdnor(i,j) = ((n1+n4+n2+n3)/4.).dir();
                }else{
                    cmdnor(i,j) = ((n3+n4)/2.).dir();
                }
            }else{
                cmdnor(i,j) = cmnor(i,j);
            }
        }

        //3rd edge
        j = nmax-1;
        for(int i=1; i<mmax-1; i++) {
            double x = (cmxy(i,j)[0]);
            double y = (cmxy(i,j)[1]);
            if(x<=1 && x>=0 && y<=1 && y>=0){
                Point3 n2 = cross((cmdpos(i,j-1) - cmdpos(i,j)),(cmdpos(i+1,j)-cmdpos(i,j))).dir();
                Point3 n3 = cross((cmdpos(i-1,j) - cmdpos(i,j)),(cmdpos(i,j-1)-cmdpos(i,j))).dir();
    
                if (calcBoundryPos(i,j,f,2, cmdposnew)){
                    Point3 n1 = cross((cmdposnew[2] - cmdposnew[0]),(cmdposnew[4] - cmdposnew[0])).dir();
                    Point3 n4 = cross((cmdposnew[4] - cmdposnew[0]),(cmdposnew[1] - cmdposnew[0])).dir();
                    cmdnor(i,j) = ((n1+n4+n2+n3)/4.).dir();
                }else{
                    cmdnor(i,j) = ((n2+n3)/2.).dir();
                }
            }else{
                cmdnor(i,j) = cmnor(i,j);
            }
        }

        //update normals at the corners
        double ttl = pow2(_lvl);
        double dd =1./ttl;
        intpair ijPair[4];
        ijPair[0].first = 0;      ijPair[0].second = 0;
        ijPair[1].first = mmax-1; ijPair[1].second = 0;
        ijPair[2].first = mmax-1; ijPair[2].second = nmax-1;
        ijPair[3].first = 0;      ijPair[3].second = nmax-1;

        Point2 ddVal[3];
        ddVal[0] = Point2(0,0);    ddVal[1] = Point2(dd,0);    ddVal[2] = Point2(0,dd);

        for (int ver=0; ver<4; ver++){

            //      cout << "Vertex " <<ver << " of face " << _textureFaces[f] << endl;
            i=ijPair[ver].first; j=ijPair[ver].second;

            double x = (cmxy(i,j)[0]);
            double y = (cmxy(i,j)[1]);
            if(x<=1 && x>=0 && y<=1 && y>=0){
                int V = _bdsurf->gpmesh().Fv2Vf(_textureFaces[f],ver).first;
                int K = _bdsurf->gpmesh().valence(V);
                Point3 nn(0,0,0);
                for (int ff=0; ff<K; ff++){
                    int F = _bdsurf->gpmesh().Vf2Fv(V,ff).first;
   
                    int texFace = -1;
                    for (int k=0; k<nfaces; k++){
                        if(_textureFaces[k] == F){
                            texFace = k;
                            break;
                        }
                    }
                    assert(texFace > -1);
                    Matrix<Point3>& cmdposnew = _Vfdpos[texFace];

                    int ni, nj;
                    Point2 Vcd(0,0);
                    Point2 Fcd(0,0);

                    for(int k=0; k<3; k++){

                        Vcd = ddVal[k];

                        _bdsurf->Vfcd2Fcd(BdSurf::EVAL_VALUE, V, ff, Vcd, F, Fcd);
                        ni = (int)(Fcd[0]*ttl);
                        nj = (int)(Fcd[1]*ttl);
                        _newPos[k] = cmdposnew(ni,nj);
                        //     cout << Vcd << "\t" << ni << "\t" << nj << endl;
                    }
                    Point3 n = cross((_newPos[1] - _newPos[0]), (_newPos[2] - _newPos[0])).dir();
                    //cout << n << endl;
    
                    nn += n;
                }
                nn /= K;
                cmdnor(i,j) = nn.dir();

                //cout << "Normal at vertex " << ver <<" global num " << V << " face " << _textureFaces[f] << " is " << cmdnor(i,j) << endl;
            }
        }
    }
}

bool Texture::calcBoundryPos(int i, int j, int of, int oe, vector<Point3>& pos){

    double ttl = pow2(_lvl);
    Point2 nVCd;  Point2 oVCd;  Point2 oCd;  Point2 nCd;
    Point2 oCdim;  Point2 oCdip;  Point2 oCdjm;  Point2 oCdjp;

    Point2 oCdall[5];

    double dd =1./ttl;
    Point2 im(-dd,0);
    Point2 ip(dd ,0);
    Point2 jm(0,-dd);
    Point2 jp(0, dd);
  
    oCdall[0][0] = oCd[0] = (double)i/ttl;
    oCdall[0][1] = oCd[1] = (double)j/ttl;
  
    oCdall[1] = oCdim = oCd + im; //i-1,j
    oCdall[2] = oCdip = oCd + ip; //i+1,j
    oCdall[3] = oCdjm = oCd + jm; //i,j-1
    oCdall[4] = oCdjp = oCd + jp; //i,j+1
  
    int gf = _textureFaces[of];
    int Fnew = _bdsurf->gpmesh().Fe2Fe(gf,oe).first;
    int enew = _bdsurf->gpmesh().Fe2Fe(gf,oe).second;
    
    if (Fnew == -1 && enew == -1) { // boundary edge
        return false;
    }
    int Vnew = _bdsurf->gpmesh().Fv2Vf(Fnew, enew).first;
    int fnew = _bdsurf->gpmesh().Fv2Vf(Fnew, enew).second;

    int nfaces = _textureFaces.size();
    int texFace = -1;
    for (int k=0; k<nfaces; k++){
        if(_textureFaces[k] == Fnew){
            texFace = k;
            break;
        }
    }
    assert(texFace > -1);
    Matrix<Point3>& cmdposnew = _Vfdpos[texFace];

    int V, ff;
    int ni, nj;
    ni = nj = -9;

    for (int i=0; i<5; i++){
        _bdsurf->Fcd2Vfcd(BdSurf::EVAL_VALUE, gf, oCdall[i], oe, V, ff, oVCd);
        nVCd[0] = 1-oVCd[0];
        nVCd[1] = -oVCd[1];
        _bdsurf->Vfcd2Fcd(BdSurf::EVAL_VALUE, Vnew, fnew, nVCd, Fnew, nCd);
        ni = (int)(nCd[0]*ttl);
        nj = (int)(nCd[1]*ttl);
        if(ni >= 0 && ni <= ttl && nj >= 0 && nj <= ttl){
            pos[i] = cmdposnew(ni,nj);
        }else{
            pos[i] = Point3(0,0,0);
        }
    }
    return true;
}

/*
 * weight function
 */
double Texture::weight(double w[4], double x[4]){
    double normcoeff = (w[0] + w[1] + w[2] + w[3]);
    if (normcoeff > 0.01){
        return (x[0]*w[0]+x[1]*w[1]+x[2]*w[2]+x[3]*w[3]) / normcoeff;
    }else{
        //    return (x[0]*w[0]+x[1]*w[1]+x[2]*w[2]+x[3]*w[3]) / 0.01;
        return 1e6;
    }
}

/*
 * rotates a given vector
 */
void Texture::rotateVector(double ang, double in[2], double* out){
    out[0] = in[0] * cos(ang) + in[1] * sin(ang);
    out[1] = in[1] * cos(ang) - in[0] * sin(ang);
}

/*
 * moves the texture center by d
 */
void Texture::moveTextureCenter(double d[2]){

    double pi = acos(-1.);
    bool faceChanged = false;
    _moving = true;

    double dx = d[0];    // movement in x direction
    double dy = d[1];    // movement in y direction

    double cd[6];        // current location of texture center in Fcd coordinates.
    double Vcd[6];       // current location of texture center in Vcd coordinates.

    int F = _textureFaceNumber;    // face number on which texture center is currently located.
    int v = _centerVertex;         // local vertex number according to which the displacement is defined.
    int V;
    int f;

    cd[0] = _textureCenter[0];
    cd[1] = _textureCenter[1];

    // convert to local vertex coordinates
    _bdsurf->Fcd2Vfcd(BdSurf::EVAL_VALUE, F, cd, v, V, f, Vcd);

    //  cout <<"Before " << cd[0] << "\t" << cd[1] <<"\t" << Vcd[0] <<  "\t" << Vcd[1] <<endl;

    // aply displacement
    // multiply by 180 degrees rotated disp matrix
    Vcd[0] += - _dispMatrix[0]*dx - _dispMatrix[1]*dy;
    Vcd[1] += - _dispMatrix[2]*dx - _dispMatrix[3]*dy;

    // get back to face coordinates
    _bdsurf->Vfcd2Fcd(BdSurf::EVAL_VALUE, V, f, Vcd, F, cd);

    //  cout << "After " << cd[0] << "\t" << cd[1] <<"\t" << Vcd[0] <<  "\t" << Vcd[1] <<endl;

    /*
     * decide if face is changed
     */
    int e = -1;
    int Fnew;
    int enew;
    int vnew;
    int Vnew;
    int fnew;
    if(cd[0] > (1.0+ SCL_EPS)){
        // going out from 1st edge
        e = 1;
        faceChanged = true;
    }else if(cd[0] < (0.0 - SCL_EPS)){
        // going out from 3rd edge
        e = 3;
        faceChanged = true;
    }else if(cd[1] < (0.0 - SCL_EPS)){
        // going out from 0th edge
        e = 0;
        faceChanged = true;
    }else if(cd[1] > (1.0 + SCL_EPS)){
        // going out from 2nd edge
        e = 2;
        faceChanged = true;
    }else{
        // still in same face
        faceChanged = false;
    }

    if(faceChanged){
        Fnew = _bdsurf->gpmesh().Fe2Fe(F,e).first;
        enew = _bdsurf->gpmesh().Fe2Fe(F,e).second;
        vnew = enew;
        Vnew = _bdsurf->gpmesh().Fv2Vf(Fnew, vnew).first;
        fnew = _bdsurf->gpmesh().Fv2Vf(Fnew, vnew).second;

        // change to local vertex coordinates of the edge on which face is changed
        int vt = e;
        _bdsurf->Fcd2Vfcd(BdSurf::EVAL_VALUE| BdSurf::EVAL_1ST_DERIV, F, cd, vt, V, f, Vcd);
    
        // update displacement matrix
        double dM[4];
        int k = vt - _centerVertex;
        double a = pi*k/2.;
        dM[0] = cos(a);
        dM[1] = sin(a);
        dM[2] = -sin(a);
        dM[3] = cos(a);
        double ndM[4];
        ndM[0] = -_dispMatrix[0] * dM[0] - _dispMatrix[1] * dM[2];
        ndM[1] = -_dispMatrix[0] * dM[1] - _dispMatrix[1] * dM[3];
        ndM[2] = -_dispMatrix[2] * dM[0] - _dispMatrix[3] * dM[2];
        ndM[3] = -_dispMatrix[2] * dM[1] - _dispMatrix[3] * dM[3];

        _dispMatrix[0] = ndM[0];
        _dispMatrix[1] = ndM[1];
        _dispMatrix[2] = ndM[2];
        _dispMatrix[3] = ndM[3];

        double newVCd[2];   // new Vcd coordinates in the new face
        double newCd[2];    // new Fcd coordinates in the new face
        newVCd[0] = 1-Vcd[0];
        newVCd[1] = -Vcd[1];

        //   cout << "Locals " << Vcd[0] << "\t" << Vcd[1] <<"\t" << newVCd[0] <<  "\t" << newVCd[1] <<endl;

        // convert these coordinates to Fcd
        _bdsurf->Vfcd2Fcd(BdSurf::EVAL_VALUE, Vnew, fnew, newVCd, Fnew, newCd);
        _textureCenter[0] = newCd[0];
        _textureCenter[1] = newCd[1];

        // cout << endl << "Face changed from " << _textureFaceNumber << " to " << Fnew << " through edge " << e << endl;
        //    cout << "Coordinates (cd--> "<< cd[0] << "\t" << cd[1] << endl;

        _textureFaceNumber = Fnew;
        _centerVertex = vnew;
    }else{
        _textureCenter[0] = cd[0];
        _textureCenter[1] = cd[1];
    }
    //  cout << "New texture center " << _textureCenter[0] << "\t" << _textureCenter[1] << endl;
}

/*
 * From Nate Robins OpenGL tutorials
 *
 * glmReadPPM: read a PPM raw (type P6) file.  The PPM file has a header
 * that should look something like:
 *
 *    P6
 *    #void  Texture::unpack(int v, int r[3]){
 r[0] = (v >> power) & (base-1);
 r[1] = (v >> (power/2)) & (base-1);
 r[2] = (v) & (base-1);
 }
 comment
 *    width height max_value
 *    rgbrgbrgb...
 *
 * where "P6" is the magic cookie which identifies the file type and
 * should be the only characters on the first line followed by a
 * carriage return.  Any line starting with a # mark will be treated
 * as a comment and discarded.   After the magic cookie, three integer
 * values are expected: width, height of the image and the maximum
 * value for a pixel (max_value must be < 256 for PPM raw files).  The
 * data section consists of width*height rgb triplets (one byte each)
 * in binary format (i.e., such as that written with fwrite() or
 * equivalent).
 *
 * The rgb data is returned as an array of unsigned chars (packed
 * rgb).  The malloc()'d memory should be free()'d by the caller.  If
 * an error occurs, an error message is sent to stderr and NULL is
 * returned.
 *
 * filename   - name of the .ppm file.
 * width      - will contain the width of the image on return.
 * height     - will contain the height of the image on return.
 *
 */
GLubyte* Texture::readTexture(char* filename, int* width, int* height)
{
    FILE* fp;
    int i, w, h, d;
    unsigned char* image;
    unsigned char* kimage;
    char head[70];          /* max line <= 70 in PPM (per spec). */
    
    fp = fopen(filename, "rb");
    if (!fp) {
        perror(filename);
        return NULL;
    }
    
    /* grab first two chars of the file and make sure that it has the
       correct magic cookie for a raw PPM file. */
    fgets(head, 70, fp);
    if (strncmp(head, "P6", 2)) {
        fprintf(stderr, "%s: Not a raw PPM file\n", filename);
        return NULL;
    }
    
    /* grab the three elements in the header (width, height, maxval). */
    i = 0;
    while(i < 3) {
        fgets(head, 70, fp);
        if (head[0] == '#')     /* skip comments. */
            continue;
        if (i == 0)
            i += sscanf(head, "%d %d %d", &w, &h, &d);
        else if (i == 1)
            i += sscanf(head, "%d %d", &h, &d);
        else if (i == 2)
            i += sscanf(head, "%d", &d);
    }
    
    /* grab all the image data in one fell swoop. */
    image = (unsigned char*)malloc(sizeof(unsigned char)*w*h*3);
    fread(image, sizeof(unsigned char), w*h*3, fp);
    fclose(fp);

    kimage = (unsigned char*)malloc(sizeof(unsigned char)*w*h*4);
    for (int i=0; i<w*h; i++){
        for (int j=0; j<3; j++){
            kimage[i*4+j] = image[i*3+j];
        }
        if(kimage[i*4] == 255 && kimage[i*4+1] == 255 && kimage[i*4+2] == 255){
            kimage[i*4+3] = 0; // convert white to transparent
        }else{
            if(_texturectrl== TEX_TYPE_DISP){
                kimage[i*4+3] = 0;
            }else{
                kimage[i*4+3] = 255;
            }
        }
    }
    *width = w;
    *height = h;
    return kimage;
}

void Texture::calcDisplacementMap(){
    int npix = _texSize[0]*_texSize[1];
    _dispMap = (double*)malloc(sizeof(double)*npix);
    for(int i=0; i<npix; i++){
        int r = _texImage[i*4];
        int g = _texImage[i*4+1];
        int b = _texImage[i*4+2];
        _dispMap[i] = 1-(double)(r*256*256 + g*256 + b) / (double)(256*256*255 + 256*255 +255);
    }
}

Matrix < Point2 >& Texture::Vfxy(int i){
   return _Vftxy[i];
}
Matrix < Point3 >& Texture::Vfpos(int i){
   if(_texturectrl == TEX_TYPE_DISP){
      return _Vfdpos[i];
   }else{
      return _Vftpos[i];
   }
}
Matrix < Point3 >& Texture::Vfnor(int i){
   if(_texturectrl == TEX_TYPE_DISP){
      return _Vfdnor[i];
   }else{
      return _Vftnor[i];
   }
}
vector < Matrix < Point2 > >& Texture::Vfxy(){
   return _Vftxy;
}
vector < Matrix < Point3 > >& Texture::Vfpos(){
   if(_texturectrl == TEX_TYPE_DISP){
      return _Vfdpos;
   }else{
      return _Vftpos;
   }
}
vector < Matrix < Point3 > >& Texture::Vfnor(){
   if(_texturectrl == TEX_TYPE_DISP){
      return _Vfdnor;
   }else{
      return _Vftnor;
   }
}
