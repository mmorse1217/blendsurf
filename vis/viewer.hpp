#ifndef _VIEWER_HPP_
#define _VIEWER_HPP_

#include "common.hpp"

class Viewer 
{
protected:
    // window data
    int _winId;
    int _width;
    int _height;
    char _title[80];
    // mouse things
    bool _leftDown;
    bool _middleDown;
    bool _rightDown;
    // modifiers
    static int _modifiers;

public:  
    Viewer(char* name, int w = 512, int h = 512);
    virtual ~Viewer();

    // window management
    int getId() const { return _winId; }
    int getWidth() { return _width; }
    int getHeight() { return _height; }
    void setSize(int w, int h);
    void setPos(int x, int y);
    void setWindow();

    bool leftDown() const { return _leftDown; }
    bool middleDown() const { return _middleDown; }
    bool rightDown() const { return _rightDown; }
    int modifiers() const { return _modifiers; }

    static void initGL(int* argc, char** argv);

protected:
    virtual void display() = 0;
    virtual void idle();
    virtual void reshape(int w, int h);
    virtual void mouse(int button, int state, int x, int y);
    virtual void motion(int x, int y);
    virtual void key(unsigned char k, int x, int y);
    virtual void specialKey(int k, int x, int y);
  

    void basemouse(int button, int state, int x, int y);

    // some wrapper functions
    static Viewer* getCurrentViewer();
    static void displayWrapper();
    static void idleWrapper();
    static void reshapeWrapper(int x, int y);
    static void mouseWrapper(int button, int state, int x, int y);
    static void motionWrapper(int x, int y);
    static void keyWrapper(unsigned char k, int x, int y);
    static void specialKeyWrapper(int k, int x, int y);
    static vector<Viewer*> _viewer;

    //protected:
    //void positionCamera();

public:
    static void redisplayAll();
};

#endif

