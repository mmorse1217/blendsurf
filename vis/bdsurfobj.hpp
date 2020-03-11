#ifndef _BDSURFOBJ_HPP_
#define _BDSURFOBJ_HPP_

#include "geoobject.hpp"
#include "bdsurf.hpp"
#include "texture.hpp"

/**
 * Rendering class for a manifold surface.
 */
class BdSurfObj: public GeoObject {
private:
    BdSurf* _bdsurf;            /// manifold surface object
    int _lvl;                   /// rendering level
    int _gen;                   /// rendering type (FULL | HALF | EXTRA_FULL)
    int _alt;                   /// derivative evaluation mode (NORMAL | HIGHER_ORDER)
    //------------
    int _renderctrl;            /// rendering type variable
    int _surfctrl;              /// surface property variable
    int _activevert;            /// active vertex number for rendering
    GLuint _chkname;            /// check board texture name
    NumMatrix _colormap;        /// color map
    vector< vector< Matrix< Point2 > > > _Vfxy;     /// x,y coordinates of a point
    vector< vector< Matrix< Point3 > > > _Vfpos;    /// surface coordinates of a point in R^3
    vector< vector< Matrix< Point3 > > > _Vfnor;    /// surface normal at a point
  
    vector< Matrix<bool> >   _Vfcgd;    /// good or not
    vector< Matrix<Point3> > _Vfcps;    /// position
    vector< Matrix<Point3> > _Vfcl0;    /// curvature line 0
    vector< Matrix<Point3> > _Vfcl1;    /// curvature line 1
  
    vector< vector< Matrix< double > > > _Vfd1; //first order magnitude
    vector< vector< Matrix< double > > > _Vfd2; //second order magnitude
    vector< vector< Matrix< double > > > _Vfd3; //third order magnitude
    vector< vector< Matrix< double > > > _Vfd4; //fourth order magnitude
    vector< vector< Matrix< double > > > _Vfd5; //fifth order magnitude
    vector< vector< Matrix< double > > > _Vfd6; //6th order magnitude
    vector< vector< Matrix< double > > > _Vfd7; //7th order magnitude
    vector< vector< Matrix< double > > > _Vfd8; //8th order magnitude
public:
    /**
     * Constructor method
     */
  BdSurfObj(BdSurf* surf, int lvl, int gen, int alt, char* texFile, 
	    int renderctrl=RENDER_SURF,int surfctrl=SURF_NONE, int activevert=-1): 
            GeoObject(), _bdsurf(surf), _lvl(lvl), _gen(gen), _alt(alt),
            _renderctrl(renderctrl),	_surfctrl(surfctrl), _activevert(activevert),
            _movementctrl(TEX_MOVE_COLOR),   _texturectrl(TEX_TYPE_IMAGE)
            {selectionMode = false; _texture = new Texture(surf,lvl,texFile);}
    /**
     * glut render callback handler
     */
    void render();
    /**
     * Method to render the manifold surface object
     */
    void renderObject();
    /**
     * glut idle callback handler
     */
    void idle() {;}
    /**
     * glut key callback handler
     * @param unsigned char     k       character pressed by user
     */
    void key(unsigned char k);
    /**
     * returns the center point of the mesh.
     */
    vec3f centerPoint() { Point3 c = _bdsurf->ctr(); return vec3f(c[0],c[1],c[2]); }
    /**
     * returns the lower left corner (min) of the bounding box of mesh.
     */
    vec3f minPoint() { Point3 a,b; _bdsurf->bbox(a,b); return vec3f(a[0],a[1],a[2]); }
    /**
     * returns the upper right corner (max) of the bounding box of the mesh.
     */
    vec3f maxPoint() { Point3 a,b; _bdsurf->bbox(a,b); return vec3f(b[0],b[1],b[2]); }
    /**
     * sets the material properties
     */
    void setMaterial();

    /// Flags for rendering type
    enum { 
        RENDER_SURF = 1,        // render surface
        RENDER_FRAME = 2,       // render mesh
        RENDER_INFLBDRY = 4,    // render influence boundary
        RENDER_CVTLINE = 8,     // render curvature line
        RENDER_CPT = 16,        // render control mesh
        RENDER_CCPTS = 32,      // render limit points
        RENDER_TEXTURE = 64     // render texture
    };
    /// flags for rendering surface properties
    enum {
        SURF_NONE = 0,          // display no surface property
        SURF_CUBEMAP = 1,       // display cubemap
        SURF_CHECKBOARD = 2,    // display checkboard
        SURF_GD1 = 3,           // 1st order magnitude
        SURF_GD2 = 4,           // 2nd order magnitude
        SURF_GD3 = 5,           // 3rd order magnitude
        SURF_GD4 = 6,           // 4th order magnitude
        SURF_GD5 = 7,           // 5th order magnitude
        SURF_GD6 = 8,           // 6th order magnitude
        SURF_GD7 = 9,           // 7th order magnitude
        SURF_GD8 = 10,           // 8th order magnitude
        //...
        SURF_TTL = 11
    };

  
private:
    typedef pair<int,int> intpair;

    Texture* _texture;          /// texture object
    float fovy;                 /// fovy angle
    float aspect;               /// aspect ratio of window
    float znear;                /// distance between near face of viewing volume and camera
    float zfar;                 /// distance between far face of viewing volume and camera
    int base;                   /// base number used in encoding names in selection rendering
    int power;                  /// 2^power = base
    int textureRenderType;      /// clamp | repeat | clamp_to_border
    int _movementctrl;          /// control movement type(uv, selection buffer, color picking)
    int _texturectrl;           /// control texture implementation (image | transperancy | displacement)
    bool selectionMode;         /// true if user is _moving the testure in selection mode
    bool updateTexture;         /// true if texture calculation is needed
    int _texSize[2];             /// size of texture image in pixels.
    int mouseOldPosition[2];    /// previous position of mouse in screen coordinates
    int mousePosition[2];       /// current position of mouse in screen coordinates
    unsigned char* _texImage;    /// points to the texture image data

    /**
     * Makes surface point selection, updates the texture object with new
     * surface position and normal
     */
    void doSelection();
    /**
     * rendering function used in color picking mode.
     */
    void pickRender();
    /**
     * selection algorithm using color picking
     * @param int       x       x coordinate of mouse in pixel coordinates
     * @param int       y       y coordinate of mouse in pixel coordinates
     * @param int&      face    reference to face number on which the user clicked(x,y)
     * @param double    cd      cd coordinates of nearest point to the point clicked
     * @param Point3    p       surface coordinates of nearest point to the point clicked
     * @param Point3    n       surface normal at nearest point to the point clicked
     */
    bool selectColor(int x, int y, int& face, double cd[2], Point3& p, Point3& n);
    /**
     * selection algorithm using selection buffer
     * @param int       x       x coordinate of mouse in pixel coordinates
     * @param int       y       y coordinate of mouse in pixel coordinates
     * @param int&      face    reference to face number on which the user clicked(x,y)
     * @param double    cd      cd coordinates of nearest point to the point clicked
     * @param Point3    p       surface coordinates of nearest point to the point clicked
     * @param Point3    n       surface normal at nearest point to the point clicked
     */
    void selectSurfacePoint(int x, int y, int& face, double cd[2], Point3& p, Point3& n);
    /**
     * Draws mesh for selection algorithm using selection buffer
     */
    void drawRects(GLenum mode);
    /**
     * process outcome of selection buffer based selection
     */
    GLuint processHits(GLint hits, GLuint buffer[]);
    /**
     * glut special key callback handler
     */
    void specialKey(unsigned char k);
    /**
     * glut mouse callback handler
     */
    void mouse(int state, int button, int x, int y);
    /**
     * glut motion callback handler
     */
    void motion(int x, int y);
    /**
     * updates the perspective paramters of the view
     * @param float     f       fovy angle
     * @param float     a       aspect ratio
     * @param float     zn      distance between near face of viewing volume and camera
     * @param float     zf      distance between far face of viewing volume and camera
     */
    void setPerspectiveParam(float f, float a, float zn, float zf);
    /**
     * unpacks given number into r,g,b components for color picking
     * @param int       v       input number
     * @param int[]     r       r,g,b
     */
  int advanceRender( ) {   
    cerr << "surfctl" << _surfctrl << endl;
    if((_surfctrl >= SURF_GD1) && (_surfctrl < SURF_TTL-1))
      { _surfctrl++; return 1; }
    else return 0; 
  }

    void unpack(int v, int r[3]);
    /// direction of rotation
    enum {
        ROTATE_CW = 1,          // rotate clock-wise
        ROTATE_CCW = 2          // rotate conter clockwise
    };

    /// direction input for texture movement using arrow keys (not always coincide with screen direction)
    enum {
        MOVE_RIGHT = 1,         // move right
        MOVE_LEFT = 2,          // move left
        MOVE_UP = 3,            // move up
        MOVE_DOWN = 4           // move down
    };
    /// texture movement control flag
    enum {
        TEX_MOVE_UV = 0,        // movement using arrow keys
        TEX_MOVE_SELECT = 1,    // movement using selection buffer picking
        TEX_MOVE_COLOR = 2,     // movement using color picking
        TEX_MOVE_TTL = 3
    };
    /// flag for representation of tecture image
    enum {
        TEX_TYPE_IMAGE = 0,     // display as image
        TEX_TYPE_TRANS = 4,     // display as alpha component (not used)
        TEX_TYPE_DISP = 1,      // display as displacement map
        TEX_TYPE_TTL = 2
    };
    /// flags for rendering type
    enum {
        RENDER_FULL = 0,        // render full mesh
        RENDER_HALF = 1,        // render half mesh
        RENDER_EXTRA = 2        // render extraordinary full charts
    };
    /// flags for derivative evaluation mode
    enum {
        EVAL_NORMAL = 0,        // normal evaluation
        EVAL_HIGH_ORDER = 1,     // higher order evaluation
        EVAL_CURVATURE = 2     // curvature evaluation
    };
};
#endif
