#ifndef TEXTURE_HPP_
#define TEXTURE_HPP_

#include "common.hpp"
#include "nummatrix.hpp"
#include "mat2t.hpp"
#include "vec2t.hpp"
#include "vec3t.hpp"
#include "bdsurf.hpp"
#include <GL/glu.h>
#include <cstring>

/**
 * Texture class
 */
class Texture{
public:
    /**
     * constructor
     */
	Texture(BdSurf* surf, int lvl, char* texFile);
    /**
     * destructor
     */
	~Texture();
    /**
     * Accessor method for Vftxy
     * @param int   i   texture face number
     */
    Matrix < Point2 >& Vfxy(int i);
    /**
     * Accessor method for Vftpos | Vfdpos
     * Depending on the rendering type (image|displacement) texture coordinates
     * (vftpos) or displacement coordinates (Vfdpos) is returned
     * @param int   i   texture face number
     */
    Matrix < Point3 >& Vfpos(int i);
    /**
     * Accessor method for Vftnor | Vfdnor
     * Depending on the rendering type (image|displacement) texture coordinates
     * (vftnor) or displacement coordinates (Vfdnor) is returned
     * @param int   i   texture face number
     */
    Matrix < Point3 >& Vfnor(int i);
    /**
     * Accessor method for Vftxy(i)
     * @param int   i   texture face number
     */
    Matrix < Point2 >& Vftxy(int i){return _Vftxy[i];}
    /**
     * Accessor method for Vftpos(i)
     * @param int   i   texture face number
     */
    Matrix < Point3 >& Vftpos(int i){return _Vftpos[i];}
    /**
     * Accessor method for Vftnor(i)
     * @param int   i   texture face number
     */
    Matrix < Point3 >& Vftnor(int i){return _Vftnor[i];}
    /**
     * Accessor method for Vftxy
     */
    vector < Matrix < Point2 > >& Vftxy(){return _Vftxy;}
    /**
     * Accessor method for Vftpos
     */
    vector < Matrix < Point3 > >& Vftpos(){return _Vftpos;}
    /**
     * Accessor method for Vftnor
     */
    vector < Matrix < Point3 > >& Vftnor(){return _Vftnor;}
    /**
     * Accessor method for Vftxy
     */
    vector < Matrix < Point2 > >& Vfxy();
    /**
     * Accessor method for Vftpos | Vfdpos
     * Depending on the rendering type (image|displacement) texture coordinates
     * (vftpos) or displacement coordinates (Vfdpos) is returned
     */
    vector < Matrix < Point3 > >& Vfpos();
    /**
     * Accessor method for Vftnor | Vfdnor
     * Depending on the rendering type (image|displacement) texture coordinates
     * (vftnor) or displacement coordinates (Vfdnor) is returned
     */
    vector < Matrix < Point3 > >& Vfnor();
    /**
     * Scales the texture
     * @param double    s   scale
     */
    void scaleTexture(double s){ _textureScale += s; }
    /**
     * Scales the displacement map
     * @param double    s   scale
     */
    void scaleDispMap(double s){ _dispMapScale += s; }
    /**
     * Acessor method for number of texture faces
     * @return int Number of texture faces
     */
    int numFaces(){return _textureFaces.size();}
    /**
     * Accessor method for texture face numbers
     * @param  int      i       texture face number
     * @return int              global face number for corresponding face
     */
    int textFace(int i){return _textureFaces[i];}
    /**
     * Stores current position information as old position information
     */
    void storeOldPos();
    /**
     * Stores given position information as current position information
     * @param int       face                face number at which texture center is located
     * @param Point2&   textCenter          texture center in cd coordinate system
     * @param Point3&   textSurfaceCenter   texture center in R^3
     * @param Point3&   textSurfNormal      surface normal if texture center
     */
    void storePos(int face, Point2& textCenter, Point3& textSurfCenter, Point3& textSurfNormal);
    /**
     * calculate texture coordinates with the updated parameters
     */
    void updateTextureParam();
    /**
     * Read the texture image file
     */
    GLubyte* read(int* width, int* height);
    /**
     * rotate the texture in given direction (ccw / cw)
     * @param int       direction       ROTATE_CW | ROTATE_CCW
     */
    void rotateTexture(int direction);
    /**
     * moves the texture in given direction with predefined amount
     * @param int       direction       MOVE_RIGHT | LEFT | UP | DOWN
     */
    void moveTexture(int direction);
    /**
     * moves the texture center by the displacement vector d (in Fcd coordinates)
     * @param double[2] d               move texture by d
     */
    void moveTextureCenter(double d[2]);
    /**
     * sets necessary flags for texture mode.
     * @param int       faceNum         global face number of texture center
     * @param Point2&   textCenter      texture center in cd coordinate system
     * @param Point3&   textSurfaceCenter   texture center in R^3
     * @param Point3&   textSurfNormal      surface normal if texture center
     */
    void setTextureMode(int faceNum, Point2& textCenter, Point3& textSurfCenter, 
                        Point3& textSurfNormal);
    /**
     * Switch texture mode between texture image and displacement map
     */
    void switchTextureMode();
    /**
     * Convert texture image colors to displacement map
     */
    void calcDisplacementMap();
    /**
     * Load manifold surface cy, position and normals
     * @param ...       xy      xy coordinates of manifold surface
     * @param ...       pos     surface coordinates of manifold surface in R^3
     * @param ...       nor     surface normals of manifold surface
     */
    void initData(vector< vector< Matrix< Point2 > > >& xy,
                  vector< vector< Matrix< Point3 > > >& pos,
                  vector< vector< Matrix< Point3 > > >& nor);
    /**
     * remove texture
     */
    void remove();
private:
    typedef pair<int,int> intpair;

    BdSurf* _bdsurf;            // manifold surface object
    int _lvl;                   // rendering level
    vector< vector< Matrix< Point2 > > > _Vfxy;     // xy coordinates of manifold surface
    vector< vector< Matrix< Point3 > > > _Vfpos;    // surface coordinates of manifold surface
    vector< vector< Matrix< Point3 > > > _Vfnor;    // surface normals of manifold surface
    vector < Matrix < Point2 > > _Vftxy;            // texture coordinates
    vector < Matrix < Point3 > > _Vftpos;           // corresponding surface coordinates
    vector < Matrix < Point3 > > _Vftnor;           // corresponding surface normal
    vector < Matrix < double > > _alpha;            // alpha value
    vector < Matrix < Point2 > > _Vfdxy;            // displacement coordinates
    vector < Matrix < Point3 > > _Vfdpos;           // position value for displacement map
    vector < Matrix < Point3 > > _Vfdnor;           // surface normal for displacement map

    int _texturectrl;            // control texture implementation (image | transperancy | displacement)
    bool _moving;                // true if the texture is currently being dragged in u,v mode
    int _oldTextFaceNumber;      // the face number on which texture center was placed
    int _textureFaceNumber;      // the face number on which texture center is placed
    double _textureScale;        // the scale of the texture 0--1
    double _dispMapScale;        // the scale of the displacement map
    double _newTextureRotation;  // rotation angle of texture (relative to initial orientation)
    double _oldTextureRotation;  // previous rotation angle of texture (relative to initial orientation)
    int _texSize[2];             // size of texture image in pixels.
    Point2 _textureOrientation;  // texture orientation in Fcd coordinates
    Point3 _textureSurfOrient;   // texture orientation on the surface
    Point3 _oldTextureNormal;    // normal direction at texture center at previous location
    Point3 _textureSurfNormal;   // normal direction at texture center at current location
    Point3 _oldTextureSurfCenter;// previous center of texture in surface coordinates
    Point3 _textureSurfCenter;   // center of texture in surface coordinates
    Point2 _oldTextureCenter;    // previous texture center in Fcd coordinates
    Point2 _textureCenter;       // texture center in Fcd coordinates
    int _centerVertex;           // local vertex number according to which the displacement is defined
    double _dispMatrix[4];       // rotation matrix for displacement vector when changing faces
    double _orientMatrix[4];     // rotation matrix for orientation vector when changing faces
    vector<int> _textureFaces;   // stores the faces that are in the union of 4 charts sharing face F
    unsigned char* _texImage;    // points to the texture image data
    double* _dispMap;            // displacement map data
    GLuint _texname;             // texture name
    char _texFile[100];          // texture file name

    /**
     * applies the displacement map calculated from texture image to surface coordinate
     * and normals
     */
    void applyDisplacementMap();
    /**
     * Calculates the surface coordinates of a boundary point and its neightbor points.
     * Returns false if the given edge is a boundary edge.
     * @param int       i       x coordinate of a point in Fcd coordinates in of
     * @param int       j       x coordinate of a point in Fcd coordinates in of
     * @param int       of      texture face number
     * @param int       oe      edge number of the corresponding global face
     * @param vector<Point3>& pos   position coordinates of center and 4 neighbor points
     */
    bool calcBoundryPos(int i, int j, int of, int oe, vector<Point3>& pos);
    /**
     * updates the orientation of the texture when a face change occurs (of to nf)
     * @param int       nf      new global face number
     * @param int       of      old global face number
     */
    void updateOrientation(int nf, int of);
    /**
     * calculates orientation of texture using variation in normal direction of texture center
     */
    void getSurfOrientation();
    /**
     * calculates the texture coordinates
     * @param Point2    center  center of texture in Fcd coordinates
     * @param Point2    orient  orientation of texture in Fcd coordinates
     * @param int       face    global face number on which the texture center is located
     */
    void calculateTexture(Point2 center, Point2 orient, int face);
    /**
     * applies weighting coefficient w, to given corresponding x values
     * @param double[4] w       weight coefficients
     * @param double[4] x       values to be weighted
     */
    double weight (double w[4], double x[4]);
    /**
     * rotates a given vector in, by angle ang, and outputs new one to out
     * @param double    ang     rotation angle in radians
     * @param double[2] in      input vector
     * @param double[2] out     rotated vector
     */
    void rotateVector(double ang, double in[2], double* out);
    /**
     * rotates vector a, around vector v, by angle ang
     * @param Point3&   a       vector to be rotated
     * @param Point3&   v       rotation axis
     * @param double    ang     rotation angle in radians
     * @param Point3&   r       rotated vector
     */
    void rotateVector(Point3& a, Point3& v, double ang, Point3& r);
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

    /* 
     * From Nate Robins OpenGL tutorials
     *
     * glmReadPPM: read a PPM raw (type P6) file.  The PPM file has a header
     * that should look something like:
     *
     *    P6
     *    # comment
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
    GLubyte* readTexture(char* filename, int* width, int* height);

};

#endif // TEXTURE_HPP_
