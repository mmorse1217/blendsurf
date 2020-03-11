#include <GL/gl.h>
#include "psOpenGL.h"
#include <GL/glu.h>
#include <GL/glut.h>

#include "glcheck.hpp"
#include "ballviewer.hpp"
#include "geoobject.hpp"
#include "uiaction.hpp"

extern int WriteAndQuit;

BallViewer::BallViewer(char* name, int w, int h) : Viewer3D(name, w, h),  _uiAction(0), _drawBall(false) 
{
    setLight();   //setLightAndMaterial();
    _arcball.ShowResult();
    _imageCnt = 0;
}

BallViewer::~BallViewer() 
{ ; }
/*
  void BallViewer::addObject(GeoObject* geoObject) 
  {Viewer3D::addObject(geoObject);  //positionArcBall();}*/

void BallViewer::renderObject() 
{
    // render object
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glColor3f(1.0, 1.0, 0.0);
    _obj->setPerspectiveParam(_camera->fovy(), _camera->aspect(), _camera->znear(), _camera->zfar());
    _obj->render();
    glCheck();
}

void BallViewer::display() 
{
    //if(_needPosition && (getObject() != 0))
  if(_needPosition && 
       _obj!=NULL)
	{
	    //int vp[4] = {0, 0, getWidth(), getHeight()};
	    //spositionCamera(_camera, getObject(), vp);
	    positionCamera();
	    positionArcBall();  
	    _needPosition = false;
	}


    glClearColor( 1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
    glPushMatrix();
  
    glMultMatrixd(HMatrix::Translation(_arcball.Position()));
    glMultMatrixd(_arcball.Rotation());
    glMultMatrixd(HMatrix::Translation(-_arcball.Position()));
  
    renderObject();
    glPopMatrix();
  
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
  
    if(_drawBall)
        _arcball.DrawCompleteBall(getCamera()->viewWorldPosition());
	
    glEnable(GL_DEPTH_TEST);

    glutSwapBuffers();

    glCheck();
    //    writePPM();
    //     exit(0);
    if(WriteAndQuit) { 
      writePPM();
      while(_obj->advanceRender()) {
	glClearColor( 1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);	
	glPushMatrix();       
	glMultMatrixd(HMatrix::Translation(_arcball.Position()));
	glMultMatrixd(_arcball.Rotation());
	glMultMatrixd(HMatrix::Translation(-_arcball.Position()));	
	renderObject();
	glPopMatrix();
	glFlush();
	glutSwapBuffers();	
	writePPM();
      }
      exit(0);
    }
}

int specKey;
void BallViewer::mouse(int button, int state, int x, int y) 
{
    y = getHeight() - y; 
    if(state == GLUT_DOWN) {
        if(_uiAction)
            delete _uiAction;
        // map left+shift to right button
        if((button == GLUT_LEFT_BUTTON) && 
           (glutGetModifiers() & GLUT_ACTIVE_ALT))
            button = GLUT_RIGHT_BUTTON;
        switch(button) 		{
            case GLUT_LEFT_BUTTON: 
                _uiAction = new CameraRotateAction(x, y, &_arcball, getCamera()); 
                break;
            case GLUT_MIDDLE_BUTTON: 
                _uiAction = new CameraTransXYAction(x, y, &_arcball, getCamera()); 
                break;
            case GLUT_RIGHT_BUTTON: 
                _uiAction = new CameraTransZAction(x, y, &_arcball, getCamera()); 
                break;
        }
    } else {
        delete _uiAction;
        _uiAction = 0;
    }
    specKey = glutGetModifiers();
    _obj->mouse(button, state, x, y);
    glutPostRedisplay();
}

void BallViewer::motion(int x, int y) 
{
    y = getHeight() - y; 

    /*
     * call _obj->motion to handle texture movement
     */
    if(specKey == GLUT_ACTIVE_SHIFT){
        _obj->motion(x,y);
    }else{
        if(_uiAction)
            _uiAction->update(x,y);
    }
    glutPostRedisplay();
}

void BallViewer::key(unsigned char k, int , int) 
{
  
  switch(k) {
    case 'q':	 // quit the viewer
      exit(0);
      break;
    case 'c':	 // center the arc ball
      if(_uiAction) {
	delete _uiAction;
	_uiAction = 0;
      }
      centerArcBall();
      break;
    case 'w':
      writePPM();
      break;
    case 'u':
      {
      string camfname;
      cerr << "writing camera " << endl;
      camfname = string(_title) + ".cam";
      ofstream camf(camfname.c_str());
      getCamera()->write(camf);
      camf.close();
      break;
      }
    default:
      _obj->key(k);
      break;
    }
    glutPostRedisplay(); 
}

void BallViewer::specialKey(int k, int /* x */, int /* y */ ) 
{
    switch(k)     {
        case GLUT_KEY_HOME:
            positionCamera();
            positionArcBall();
            break;
        default:
            _obj->specialKey(k);
            break;
    }
    glutPostRedisplay();
}

void BallViewer::centerArcBall() 
{
    float z = getCamera()->project(centerPoint()).z();
  
    float top  =  getHeight() * 5.0f/6.0f;
    float right = getWidth()  * 5.0f/6.0f;
    float cx = getWidth() * 0.5f;
    float cy = getWidth() * 0.5;
  
    vec3f pos = getCamera()->unproject(vec3f(cx, cy, z));
    float rad = min((getCamera()->unproject(vec3f(right, cy, z)) - pos).l2(),
                    (getCamera()->unproject(vec3f(cx, top, z)) - pos).l2());
  
    getArcBall()->Init();
    getArcBall()->SetPosition(pos);
    getArcBall()->SetRadius(rad);
    getArcBall()->Update();
}


void BallViewer::positionArcBall() 
{
    vec3f center = centerPoint();
    float r = 0.5f * (minPoint() - maxPoint()).l2();
    getArcBall()->Init();
    getArcBall()->SetPosition(center);
    getArcBall()->SetRadius(r);
    getArcBall()->Update();
}

void BallViewer::setLight()
{
    GLfloat lgt1_diffuse[] =  { 0.5f, 0.5f, 0.5f, 1.0f };
 
    GLfloat lgt1_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    //GLfloat lgt1_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat lgt1_ambient[] =  {0.1f, 0.1f, 0.1f, 1.0f};
  
    GLfloat lgt2_diffuse[] =  { 0.5f, 0.5f, 0.5f, 1.0f };

    GLfloat lgt2_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    //GLfloat lgt2_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat lgt2_ambient[] =  {0.1f, 0.1f, 0.1f, 1.0f};
  
    GLfloat light_pos0[] =   { 1.0f,  1.0f, 1.0f, 0.0f };
    GLfloat light_pos1[] =   { 1.0f,  1.0f, 0.0f, 0.0f };
    GLfloat light_pos2[] =   { -1.0f, 1.0f, 0.0f, 0.0f };
  
    glEnable(GL_LIGHTING);
  
    glLightfv(GL_LIGHT0, GL_POSITION,light_pos0);
  
    glLightfv(GL_LIGHT1, GL_POSITION,light_pos1);
    glLightfv(GL_LIGHT1, GL_SPECULAR,lgt1_specular);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,lgt1_diffuse);
    glLightfv(GL_LIGHT1, GL_AMBIENT, lgt1_ambient);
  
    glLightfv(GL_LIGHT2, GL_POSITION,light_pos2);
    glLightfv(GL_LIGHT2, GL_SPECULAR,lgt2_specular);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, lgt2_diffuse);
    glLightfv(GL_LIGHT2, GL_AMBIENT, lgt2_ambient);
  
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
  
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
  
    glCheck();
}

//---------------------------------------------------
void BallViewer::writePPM()
{
    GLubyte *outputImage = (GLubyte*) malloc(_width*_height*3*sizeof(GLubyte));
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, _width, _height, GL_RGB, GL_UNSIGNED_BYTE, outputImage);
    char filename[100];
    sprintf(filename, "%s_%d.ppm", _title, _imageCnt); _imageCnt++; cerr<<filename<<endl;
    GLubyte* tmp = new GLubyte[3*_width];
    for(int i = 0; i < _height/2; ++i) {
        memcpy(tmp, &outputImage[3*_width*i], 3*_width);
        memcpy(&outputImage[3*_width*i], &outputImage[3*_width*(_height-1-i)], 3*_width);
        memcpy(&outputImage[3*_width*(_height-1-i)], tmp, 3*_width);
    }
    delete[] tmp;
    int max = 255;
    FILE* file;  file = fopen(filename,"w");
    if (file != 0) {
        fprintf(file,"P6\n");
        fprintf(file,"# created by subdivision software, NYU\n");
        fprintf(file, "%d ", _width);
        fprintf(file, "%d ", _height);
        fprintf(file, "%d\n", max);
        fwrite(outputImage, sizeof(unsigned char), 3*_width*_height, file);
        fclose(file);
    } else {
        printf("error writing ppm image: %s\n",filename);
    }
    delete outputImage;
}
