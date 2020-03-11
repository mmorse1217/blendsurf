//#include "CCcommon.h"
//#include "compat.h"
#include <GL/gl.h>
#include "psOpenGL.h"
#include <GL/glu.h>
#include <GL/glut.h>


#include "glcheck.hpp"
#include "viewer.hpp"

int Viewer::_modifiers;

vector<Viewer*> Viewer::_viewer;

Viewer::Viewer(char* t, int w, int h) : 
        _width(w), _height(h), _leftDown(false), _middleDown(false), _rightDown(false) {
  
    if(t == 0) 
        sprintf(_title, "Viewer#%d", _viewer.size());
    else
        strcpy(_title, t);

    _viewer.push_back(this);

    glutInitWindowPosition(100,100);
    glutInitWindowSize(_width, _height);

    _winId = glutCreateWindow(_title);
    glutSetWindow(getId());
  
    glutDisplayFunc(displayWrapper);
    glutIdleFunc(idleWrapper);
    glutMouseFunc(mouseWrapper);
    glutMotionFunc(motionWrapper);
    glutReshapeFunc(reshapeWrapper);
    glutKeyboardFunc(keyWrapper);
    glutSpecialFunc(specialKeyWrapper);
  
    glCheck();
}


void Viewer::initGL(int* argc, char** argv) {
    glutInit(argc, argv);
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
    //  glCheck();
}

Viewer::~Viewer() { ; }


void Viewer::setSize(int w, int h) {
    reshape(w, h);
    glCheck();
}


void Viewer::setPos(int x, int y) {
    glutSetWindow(getId());
    glutPositionWindow(x, y);
    glutPostRedisplay();
    glCheck();
}

void Viewer::idle() {
    //do nothing
}

void Viewer::reshape(int w, int h) {
    _width = w;
    _height = h;
}

void Viewer::mouse(int , int, int, int ) {
    glutPostRedisplay();
    glCheck();
}


void Viewer::motion(int , int) {
    glutPostRedisplay();
    glCheck();
}


void Viewer::key(unsigned char , int , int ) {
    glutPostRedisplay();
    glCheck();
}


void Viewer::specialKey(int, int, int) {
    glutPostRedisplay();
    glCheck();
}


void Viewer::redisplayAll() {
    for(int i=0; i<_viewer.size(); ++i) {
        glutSetWindow(_viewer[i]->getId());
        glutPostRedisplay();
    }
    glCheck();
}


void Viewer::setWindow() {
    glutSetWindow(getId());
    glutPostRedisplay();
}

//--------------------------------------------------------------
// dispatcher

Viewer* Viewer::getCurrentViewer() {
    int id = glutGetWindow();
    for(int i = 0; i < Viewer::_viewer.size(); ++i)
        if(_viewer[i]->getId() == id)
            return _viewer[i];
    return 0;
}

void Viewer::displayWrapper() {
    Viewer* v = getCurrentViewer();
    if(v) v->display();
}

void Viewer::idleWrapper() {
    Viewer* v = getCurrentViewer();
    if(v) v->idle();
}

void Viewer::reshapeWrapper(int x, int y) {
    Viewer* v = getCurrentViewer();
    if(v) v->reshape(x, y);
}

void Viewer::mouseWrapper(int button, int state, int x, int y) {
    _modifiers = glutGetModifiers();
    Viewer* v = getCurrentViewer();
    if(v) v->basemouse(button, state, x, y);
}

void Viewer::motionWrapper(int x, int y) {
    //_modifiers = glutGetModifiers();
    Viewer* v = getCurrentViewer();
    if(v) v->motion(x, y);
}

void Viewer::keyWrapper(unsigned char k, int x, int y) {
    Viewer* v = getCurrentViewer();
    if(v) v->key(k, x, y);
}

void Viewer::specialKeyWrapper(int k, int x, int y) {
    Viewer* v = getCurrentViewer();
    if(v) v->specialKey(k, x, y);
}

void Viewer::basemouse(int button, int state, int x, int y) {
    if(button == GLUT_LEFT_BUTTON) 
        _leftDown = (state == GLUT_DOWN);
    else if(button == GLUT_MIDDLE_BUTTON)
        _middleDown = (state == GLUT_DOWN);
    else if(button == GLUT_RIGHT_BUTTON)
        _rightDown = (state == GLUT_DOWN);
    mouse(button, state, x, y);
}
