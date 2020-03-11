#include "ccsurfobj.hpp"
#include "reflmap.h"

//---------------------------------------------------
void CCSurfObj::render()
{
    //------------------todo at first rendering
    if(_lmt.size()==0) {
        //1. make cubemap
        checkCubeMapExtension();  makeCubeMap(1024, 64); cubeMapOff();
        //2. make checkboard
        int chksz=8;
        Matrix<GLubyte> chkimage(chksz,chksz);
        for(int i=0; i<chksz; i++)
            for(int j=0; j<chksz; j++) {
                chkimage(i,j) = ((i%2==0)^(j%2==0)) * 255;
            }
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glGenTextures(1, &_chkname);
        glBindTexture(GL_TEXTURE_2D, _chkname);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, chksz, chksz,	0, GL_LUMINANCE, GL_UNSIGNED_BYTE, chkimage.data());
        //3. generate data
        CCSurfOp::limit(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _lmt, *(_ccSurf->subdivMatrices()) );
        CCSurfOp::du(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _du );
        CCSurfOp::dv(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _dv );
        CCSurfOp::normal(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _nor, *(_ccSurf->subdivMatrices()) );
        CCSurfOp::duu(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _duu );
        CCSurfOp::duv(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _duv );
        CCSurfOp::dvv(_ccSurf->gpmesh(), _lvl, _ccSurf->pos(_lvl), _dvv );
    }
  
    //compute active faces
    GpMesh& gpmesh = _ccSurf->gpmesh();
    vector<int> activeFs;  
    for(int f=0; f<gpmesh.valence(_activevert); f++)
        activeFs.push_back( gpmesh.Vf2Fv(_activevert,f).first );
  
    double zOffset = 3e-4;
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_SURF) {
        glDepthRange(zOffset,1.0);//	 	 if(_mapOn==1) cubeMapOn();
        //glCullFace(GL_BACK);
        //glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        setMaterial();	 
	glColor3f(0.5,0.5,1);
        int size = pow2(_lvl)+1;
        if(       _surfctrl==SURF_NONE) { //
            //---------------------------------------
            glEnable(GL_LIGHTING);
            for(int F=0; F<_ccSurf->numFaces(); F++) {
                CCRect<Point3>& cmlmt = _lmt[F];
                CCRect<Point3>& cmnor = _nor[F];
                for(int j=0; j<size-1; j++) {
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0; i<size; i++) {
                        glNormal3dv(cmnor(i,j+1));
                        glVertex3dv(cmlmt(i,j+1));
                        glNormal3dv(cmnor(i,j));
                        glVertex3dv(cmlmt(i,j));
                    }
                    glEnd();
                }
            }
            glFlush();
        } else if(_surfctrl==SURF_CUBEMAP) { //
            //---------------------------------------
            
		
            glEnable(GL_LIGHTING);
	    glMatrixMode(GL_TEXTURE);
            glScalef(0.5, 0.5, 1.0);
	    cubeMapOn();
            for(int F=0; F<_ccSurf->numFaces(); F++) {
                CCRect<Point3>& cmlmt = _lmt[F];
                CCRect<Point3>& cmnor = _nor[F];
                for(int j=0; j<size-1; j++) {
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0; i<size; i++) {
                        glNormal3dv(cmnor(i,j+1));
			glVertex3dv(cmlmt(i,j+1));				
			glNormal3dv(cmnor(i,j));
			glVertex3dv(cmlmt(i,j));
                    }
                    glEnd();
                }
            }
            cubeMapOff();
		
            glLoadIdentity();
            glMatrixMode(GL_MODELVIEW);
            glFlush();
        } else if(_surfctrl==SURF_CHECKBOARD) { //
            //---------------------------------------
            glEnable(GL_LIGHTING);
            for(int F=0; F<_ccSurf->numFaces(); F++) {		  
                CCRect<Point3>& cmlmt = _lmt[F];
                CCRect<Point3>& cmnor = _nor[F];
                if(find(activeFs.begin(), activeFs.end(), F)!=activeFs.end()) {
                    glEnable(GL_TEXTURE_2D);
                    for(int j=0; j<size-1; j++) {
                        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
                        glBindTexture(GL_TEXTURE_2D, _chkname);
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<size; i++) {
                            glTexCoord2d( double(i)/(size-1), double(j+1)/(size-1) );				glNormal3dv(cmnor(i,j+1));				glVertex3dv(cmlmt(i,j+1));
                            glTexCoord2d( double(i)/(size-1), double(j)/(size-1) );				glNormal3dv(cmnor(i,j));				glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                    glDisable(GL_TEXTURE_2D);
                } else {
                    for(int j=0; j<size-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<size; i++) {
                            glNormal3dv(cmnor(i,j+1));				glVertex3dv(cmlmt(i,j+1));
                            glNormal3dv(cmnor(i,j));				glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                }
            }
            glDisable(GL_LIGHTING);
            glDepthRange(0,1.0-zOffset);
            glLineWidth(2.0);
            glColor3f(0.1,0.4,0.0);
            for(int F=0; F<_ccSurf->numFaces(); F++) {		  
                CCRect<Point3>& cmlmt = _lmt[F];
                if(find(activeFs.begin(), activeFs.end(), F)!=activeFs.end()) {
                    //some lines
                    glBegin(GL_LINE_STRIP);		for(int i=0; i<size; i++)		  glVertex3dv(cmlmt(i,0));		glEnd();
                    glBegin(GL_LINE_STRIP);		for(int i=0; i<size; i++)		  glVertex3dv(cmlmt(i,size-1));		glEnd();
                    glBegin(GL_LINE_STRIP);		for(int j=0; j<size; j++)		  glVertex3dv(cmlmt(0,j));		glEnd();
                    glBegin(GL_LINE_STRIP);		for(int j=0; j<size; j++)		  glVertex3dv(cmlmt(size-1,j));		glEnd();
                }
            }
            glFlush();
        } else if(_surfctrl>=SURF_GX && _surfctrl<=SURF_GYY) {
            //---------------------------------------
            glDisable(GL_LIGHTING);
            vector< CCRect<Point3> >* dp = NULL;
            switch(_surfctrl) {
                case SURF_GX: dp = &_du; break;
                case SURF_GY: dp = &_dv; break;
                case SURF_GXX:dp = &_duu; break;
                case SURF_GXY:dp = &_duv; break;
                case SURF_GYY:dp = &_dvv; break;
            }
            vector< CCRect<Point3> >& data = *dp;
            for(int F=0; F<_ccSurf->numFaces(); F++) {
                CCRect<Point3>& cmlmt = _lmt[F];
                CCRect<Point3>& cmnor = _nor[F];
                if(find(activeFs.begin(), activeFs.end(), F)!=activeFs.end()) {
                    //get min and max
                    double mmin = SCL_MAX, mmax = -SCL_MAX;
                    CCRect<double> mags(size, size);
                    for(int i=0; i<size; i++)
                        for(int j=0; j<size; j++) {
                            double gg = data[F](i,j).length();
                            mmin = min(mmin, gg);				mmax = max(mmax, gg);
                            mags(i,j) = gg;
                        }
                    mmin -= 1e-2;		  mmax += 1e-2;
                    //mmin = 2*mmin - mmax; //make min = grey
                    for(int j=0; j<size-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        double tmp;
                        for(int i=0; i<size; i++) {
                            tmp = (mags(i,j+1)-mmin)/(mmax-mmin);				glColor3f(tmp,tmp,tmp);				glVertex3dv(cmlmt(i,j+1));
                            tmp = (mags(i,j)-mmin)/(mmax-mmin); 				glColor3f(tmp,tmp,tmp);				glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                } else {
                    for(int j=0; j<size-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<size; i++) {
                            glNormal3dv(cmnor(i,j+1));				glVertex3dv(cmlmt(i,j+1));
                            glNormal3dv(cmnor(i,j));				glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                }
                glFlush();
            }
        }
	/*
        glDepthRange(zOffset,1.0); 
        //glCullFace(GL_BACK);  //in case of boundary can see back face so don't cull.
        //glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        setMaterial();

        size = pow2(_lvl)+1;
        if(       _surfctrl==SURF_NONE) { //SHOW ALL
            //---------------------------------------
            glEnable(GL_LIGHTING);
            for(int F=0; F<_ccSurf->numFaces(); F++) { 
                //replace lmt with pos in case I want to see the current position rather than the limit..
                CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];//
                CCRect<Point3>& cmnor = _nor[F];
                for(int j=0; j<size-1; j++) {
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0; i<size; i++) {
                        glNormal3dv(cmnor(i,j+1)); glVertex3dv(cmlmt(i,j+1));
                        glNormal3dv(cmnor(i,j));   glVertex3dv(cmlmt(i,j));
                    }
                    glEnd();
                }
            }
            glFlush();
        } else if(_surfctrl==SURF_CUBEMAP) { //SHOW ALL
            //---------------------------------------
            glEnable(GL_LIGHTING);
            cubeMapOn();
            for(int F=0; F<_ccSurf->numFaces(); F++) {
                CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];// 
                CCRect<Point3>& cmnor = _nor[F];
                for(int j=0; j<size-1; j++) {
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0; i<size; i++) {
                        glNormal3dv(cmnor(i,j+1));  glVertex3dv(cmlmt(i,j+1));
                        glNormal3dv(cmnor(i,j));	glVertex3dv(cmlmt(i,j));
                    }
                    glEnd();
                }
            }
            cubeMapOff();
            glLoadIdentity();
            glMatrixMode(GL_MODELVIEW);
      
            glFlush();
        } else if(_surfctrl==SURF_CHECKBOARD) { //SHOW ONE
            //---------------------------------------
            glEnable(GL_LIGHTING);
            for(int F=0; F<_ccSurf->numFaces(); F++) {		  
                CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];//
                CCRect<Point3>& cmnor = _nor[F];
                if(find(activeFs.begin(), activeFs.end(), F)!=activeFs.end()) {
                    glEnable(GL_TEXTURE_2D);
                    for(int j=0; j<size-1; j++) {
                        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
                        glBindTexture(GL_TEXTURE_2D, _chkname);
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<size; i++) {
                            glTexCoord2d( double(i)/(size-1), double(j+1)/(size-1) );
                            glNormal3dv(cmnor(i,j+1));  glVertex3dv(cmlmt(i,j+1));
                            glTexCoord2d( double(i)/(size-1), double(j)/(size-1) );
                            glNormal3dv(cmnor(i,j));	  glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                    glDisable(GL_TEXTURE_2D);
                } else {
                    for(int j=0; j<size-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<size; i++) {
                            glNormal3dv(cmnor(i,j+1));  glVertex3dv(cmlmt(i,j+1));
                            glNormal3dv(cmnor(i,j));	  glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                }
            }
            glDisable(GL_LIGHTING);
            glDepthRange(0,1.0-zOffset);
            glLineWidth(2.0);
            glColor3f(0.1,0.4,0.0);
            for(int F=0; F<_ccSurf->numFaces(); F++) {		  
                CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];//
                if(find(activeFs.begin(), activeFs.end(), F)!=activeFs.end()) {
                    //some lines
                    glBegin(GL_LINE_STRIP);
                    for(int i=0; i<size; i++)
                        glVertex3dv(cmlmt(i,0));
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                    for(int i=0; i<size; i++)
                        glVertex3dv(cmlmt(i,size-1));
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                    for(int j=0; j<size; j++)
                        glVertex3dv(cmlmt(0,j));
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                    for(int j=0; j<size; j++)
                        glVertex3dv(cmlmt(size-1,j));
                    glEnd();
                }
            }
            glFlush();
        } else if(_surfctrl>=SURF_GX && _surfctrl<=SURF_GYY) {
            //---------------------------------------
            glDisable(GL_LIGHTING);
            vector< CCRect<Point3> >* dp = NULL;
            switch(_surfctrl) {
                case SURF_GX: dp = &_du; break;
                case SURF_GY: dp = &_dv; break;
                case SURF_GXX:dp = &_duu; break;
                case SURF_GXY:dp = &_duv; break;
                case SURF_GYY:dp = &_dvv; break;
            }
            vector< CCRect<Point3> >& data = *dp;
            for(int F=0; F<_ccSurf->numFaces(); F++) {
                CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];//
                CCRect<Point3>& cmnor = _nor[F];
                if(find(activeFs.begin(), activeFs.end(), F)!=activeFs.end()) {
                    //get min and max
                    double mmin = SCL_MAX, mmax = -SCL_MAX;
                    CCRect<double> mags(size, size);
                    for(int i=0; i<size; i++)
                        for(int j=0; j<size; j++) {
                            double gg = data[F](i,j).length();
                            mmin = min(mmin, gg);
                            mmax = max(mmax, gg);
                            mags(i,j) = gg;
                        }
                    mmin -= 1e-2;   mmax += 1e-2;
	 
                    for(int j=0; j<size-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        double tmp;
                        for(int i=0; i<size; i++) {
                            tmp = (mags(i,j+1)-mmin)/(mmax-mmin);
                            glColor3f(tmp,tmp,tmp);  glVertex3dv(cmlmt(i,j+1));
                            tmp = (mags(i,j)-mmin)/(mmax-mmin);
                            glColor3f(tmp,tmp,tmp);  glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                } else {
                    for(int j=0; j<size-1; j++) {
                        glBegin(GL_QUAD_STRIP);
                        for(int i=0; i<size; i++) {
                            glNormal3dv(cmnor(i,j+1));  glVertex3dv(cmlmt(i,j+1));
                            glNormal3dv(cmnor(i,j));	  glVertex3dv(cmlmt(i,j));
                        }
                        glEnd();
                    }
                }
                glFlush();
            }
        }
    
	*/
    }	
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_FRAME) {
        glDepthRange(0, 1-zOffset);
        glDisable(GL_LIGHTING);
        glColor3f(0.0, 0.0, 0.0);
        int size = pow2(_lvl)+1;
        for(int F=0; F<_ccSurf->numFaces(); F++) {
            CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];//
            for(int j=0; j<size; j++) {
                assert(!glCheck());
                glBegin(GL_LINE_STRIP);
                for(int i=0; i<size; i++) {
                    glVertex3dv(cmlmt(i,j));
                }
                glEnd();
            }
            for(int i=0; i<size; i++) {
                assert(!glCheck());
                glBegin(GL_LINE_STRIP);
                for(int j=0; j<size; j++) {
                    glVertex3dv(cmlmt(i,j));
                }
                glEnd();
            }
        }
        glFlush();
    }
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_INFLBDRY) {
        glDepthRange(0,1.0-zOffset);
        glDisable(GL_LIGHTING);
        int size = pow2(_lvl)+1;
        for(int F=0; F<_ccSurf->numFaces(); F++) {
            CCRect<Point3>& cmlmt = _lmt[F];//_ccSurf->pos(_lvl)[F];//
            glBegin(GL_LINE_STRIP);
            for(int i=0; i<size; i++)
                glVertex3dv(cmlmt(i,0));
            glEnd();
            glBegin(GL_LINE_STRIP);
            for(int i=0; i<size; i++)
                glVertex3dv(cmlmt(i,size-1));
            glEnd();
            glBegin(GL_LINE_STRIP);
            for(int j=0; j<size; j++)
                glVertex3dv(cmlmt(0,j));
            glEnd();
            glBegin(GL_LINE_STRIP);
            for(int j=0; j<size; j++)	
                glVertex3dv(cmlmt(size-1,j));
            glEnd();
        }
        glFlush();
    }
    //-------------------------------------------------------------------------------------
    if(_renderctrl & RENDER_CPT_FACE) {
        glDisable(GL_LIGHTING);
        //face

        glDepthRange(zOffset,1.0);
        glColor3f(0.7, 0.7, 0.7);
        glBegin(GL_QUADS);
        for(int F=0; F<gpmesh.numFaces(); F++) {
            for(int v=0; v<4; v++) {
                int V = gpmesh.Fv2Vf(F,v).first;
                glVertex3dv(gpmesh.vpoint(V).array());
            }
        }
        glEnd();
    }
    if(_renderctrl & RENDER_CPT_FACE || _renderctrl & RENDER_CPT) {
        glDisable(GL_LIGHTING);
        //line
        glDepthRange(0,1.0-zOffset);
        glColor3f(0.0, 0.0, 0.0);
        glLineWidth(2.0);
        for(int F=0; F<gpmesh.numFaces(); F++) {
            glBegin(GL_LINE_LOOP);
            for(int v=0; v<4; v++) {
                int V = gpmesh.Fv2Vf(F,v).first;
                glVertex3dv(gpmesh.vpoint(V).array());
            }
            glEnd();
        }
	
        glFlush();
    }


    if(_renderctrl & RENDER_CCPTS) {
        glDisable(GL_LIGHTING);
        int ctrllvl = 2;

        glColor3f(1.0, 0.0, 0.0);
        glPointSize(3.0);
        for(int F=0; F<_ccSurf->numFaces(); F++) { 
      
            CCRect<Point3>& cmpos = _ccSurf->pos(ctrllvl)[F];

            for(int j=0; j<pow2(ctrllvl); j++) {
                glBegin(GL_POINTS);
                for(int i=0; i<pow2(ctrllvl)+1; i++) {
                    glVertex3dv(cmpos(i,j+1));
                    glVertex3dv(cmpos(i,j));
                }
                glEnd();
            }
        
        }
        glFlush();
        setMaterial();

    } 
}

//---------------------------------------------------
void CCSurfObj::key(unsigned char k)
{
    //LEXING: q,c,w is reserved
    switch(k) {
        case 'r':
            _renderctrl = _renderctrl ^ RENDER_SURF; 	 break;
        case 'f':
            _renderctrl = _renderctrl ^ RENDER_FRAME; 	 break;
        case 'i':
            _renderctrl = _renderctrl ^ RENDER_INFLBDRY;   break;
        case 'p':
	  if( _renderctrl & RENDER_CPT )
            _renderctrl = (_renderctrl | RENDER_CPT_FACE) & ~RENDER_CPT;  
	  else if( _renderctrl & RENDER_CPT_FACE)
	  _renderctrl = _renderctrl & ~RENDER_CPT_FACE; 
	  else _renderctrl =  _renderctrl | RENDER_CPT;
	  break;
        case 'k':
            _renderctrl = _renderctrl ^ RENDER_CCPTS;   break;
        case 's':
            _surfctrl = (_surfctrl+1) % SURF_TTL;	 break;
        case 'v':
            _activevert = (_activevert+1) % _ccSurf->numVertices();	 break;
        case 'h':
            printf("r: suRf\n");
            printf("f: wireFrame\n");
            printf("i: Influence boundary\n");
            printf("p: control Points\n");
            printf("s: Surface mode\n");
            break;
    }
    glutPostRedisplay();
}



//---------------------------------------------------
void CCSurfObj::setMaterial()
{
    if(0) {
        float matSpecular[] = {0.478814, 0.457627, 0.5, 1.0};
        float matDiffuse[] =  {0.15, 0.452647, 0.154303, 1,0};
        float matAmbient[] = {0.15, 0.452647, 0.154303, 1.0};
        float matShininess = 10.0;
        float matBackSpecular[] = {0.5, 0.1, 0.75, 1.0};
        float matBackAmbient[] =  {0.15, 0.154303, 0.754303, 1.0};
        float matBackDiffuse[] =  {0.15, 0.154303, 0.754303, 1.0};
        float matBackShininess = 10.0;
        glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
        glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbient);
        glMaterialf(GL_FRONT, GL_SHININESS, matShininess);
        glMaterialfv(GL_BACK, GL_DIFFUSE, matBackDiffuse);
        glMaterialfv(GL_BACK, GL_SPECULAR, matBackSpecular);
        glMaterialfv(GL_BACK, GL_AMBIENT, matBackAmbient);
        glMaterialf(GL_BACK, GL_SHININESS, matBackShininess);
    }
    if(1) { //HENNING's material
        float matSpecular[] = {0.478814, 0.457627, 0.5};
        float matAmbient[] =  {0.75, 0.652647, 0.154303};
        float matDiffuse[] =  {0.75, 0.652647, 0.154303};
        float matFrontShininess = 25.0;
        float matBackSpecular[] = {0.1596,   0.1525,   0.1667};
        float matBackAmbient[] =  {0.3750,   0.3263,   0.0772};
        float matBackDiffuse[] =  {0.3750,   0.3263,   0.0772};
        float matBackShininess = 100.0;
        glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
        glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbient);
        glMaterialf(GL_FRONT, GL_SHININESS, matFrontShininess);


	/*        glMaterialfv(GL_FRONT, GL_DIFFUSE, matBackDiffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, matBackSpecular);
        glMaterialfv(GL_FRONT, GL_AMBIENT, matBackAmbient);
        glMaterialf(GL_FRONT, GL_SHININESS, matBackShininess);
	*/

        glMaterialfv(GL_BACK, GL_DIFFUSE, matBackDiffuse);
        glMaterialfv(GL_BACK, GL_SPECULAR, matBackSpecular);
        glMaterialfv(GL_BACK, GL_AMBIENT, matBackAmbient);
        glMaterialf(GL_BACK, GL_SHININESS, matBackShininess);	 
    }  
    glCheck();
}
