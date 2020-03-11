#ifndef _BALLVIEWER_HPP_
#define _BALLVIEWER_HPP_

#include "viewer3d.hpp"
#include "arcball.hpp"

class UiAction;

class BallViewer:public Viewer3D
{
private:
    typedef Vec3T<float> vec3f;

public:
    BallViewer(char* name = 0, int w = 512, int h = 512);
    virtual ~BallViewer();  //virtual void addObject(GeoObject* geoObject);
    ArcBall* getArcBall() { return &_arcball; }
    // shoot a ray through pixel (x,y) and return object space coordinates
    //Ray objectRay(int x, int y);
    HMatrix worldSpaceTransform() {
        return
            HMatrix::Translation(_arcball.Position()) 
            * _arcball.Rotation()
            * HMatrix::Translation(-_arcball.Position());
    }
    void showBall() { _drawBall = true; }
    void hideBall() { _drawBall = false; }  

    virtual void display();  
protected:
    UiAction* _uiAction;
    ArcBall _arcball;
    int _imageCnt;
protected:

    virtual void renderObject();
    virtual void mouse(int button, int state, int x, int y);
    virtual void motion(int x, int y);
    virtual void key(unsigned char k, int x, int y);
    virtual void specialKey(int k, int x, int y);    //virtual void writePPM();  //void writePPM();
public:
    void setLight();
    void writePPM();
  
    void centerArcBall();
    void positionArcBall();  
  
    bool _drawBall;
};

#endif
