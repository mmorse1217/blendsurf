/*
  psOpenGL, an OpenGL to PostScript conversion utility
  author: Gary I. Wu, gwu@cco.caltech.edu

  parts of this program was freely adapted from an old version of psgl.
  copyright notice for psgl follows at the end of this comment block.

    
  USAGE:  
  #include "psOpenGL.h" in all your files which uses OpenGL calls.
  Do this after the <GL/gl.h> header file have been included.

  Place the calls to psOpenGL_open() and psOpenGL_close() calls
  as deemed necessary.  All OpenGL output between these two calls
  will be trapped and converted to PostScript output and written to
  a file.  The file name is determined by the call to psOpenGL_open().

  At this point, your program should still behave as usual.  Then, 
  when you want to activate the conversion, 
  #define PSOPENGL 
  in all of your source codes before including psOpenGL.h.

  Be sure to compile and link in psOpenGL.c also (including the
  #define PSOPENGL as appropriate).

  Typically, you would call psOpenGL_open(), then generate a
  full screen refresh, and then call psOpenGL_close(), and then
  refresh the screen again.

  Restrictions:
  You should probably not change projections while doing
  conversions.
  You should probably not make calls to glRenderMode().

  If you're using OpenInventor, you may want to refresh the
  screen _before_ you call psOpenGL_open() also; this does away
  with the GLXBadContextState error.

  There may be some others. :)

  ProtoTypes:
  void psOpenGL_open( const char* filename );
  void psOpenGL_close( void );


  ----------------------------------------------------------------------------------
  NOTICE for parts of code based on psgl:
  psgl.c copyright (c) 1990, 1991, 1992  seth j. teller, seth@miro.berkeley.edu.

  permission is granted to modify and redistribute this code in any manner
  as long as this notice is preserved.  all standard disclaimers apply.
  ----------------------------------------------------------------------------------
*/    



//#ifdef WIN32
//#include <windows.h>
//#endif //WIN32

//#include <assert.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include "common.hpp"
#include "vec3t.hpp"

#include <GL/gl.h>


/* Do not remove the following two lines! */
#define PSOPENGL_CODE
#include "psOpenGL.h" 

#define TAG_DepthTest 0x0001;

#define OBJECTSIZE_INC (8 * 1024)

#define FALSE 0
#define TRUE 1

int psOpenGL_fill_gray=0;

static int             psOpenGL_BufferSize     = 4096 * 1024;
static int             psOpenGL_BufferMaxUsage = 0;
static GLfloat*        psOpenGL_Buffer         = NULL;

static int             psOpenGL_ObjectSize     = OBJECTSIZE_INC;
static int             psOpenGL_ObjectNum      = 0;
static struct OBJECT*  psOpenGL_Object         = NULL;

static int             psOpenGL_TimeStamp      = 0;

static GLint           psOpenGL_RenderMode     = GL_RENDER;

static int             psOpenGL_OriginX        = 0;
static int             psOpenGL_OriginY        = 0;
static int             psOpenGL_SizeX          = 0;
static int             psOpenGL_SizeY          = 0;

static FILE*           psOpenGL_File           = NULL;
static char            psOpenGL_FileName[512];

static char            psOpenGL_Epilogue[]     = "showpage";

static int             psOpenGL_DebugLevel     = 0;

static unsigned int    psOpenGL_Tag            = 0;

static int             psOpenGL_Enabled        = FALSE;


#ifdef _WIN32
#ifdef OBJ_BITMAP 
#undef OBJ_BITMAP 
#endif //OBJ_BITMAP 
#endif //_WIN32

enum OBJECT_TYPE
{
	OBJ_ERROR,
	OBJ_POINT,
	OBJ_LINE,
	OBJ_POLYGON,
	OBJ_BITMAP,
	OBJ_PIXEL,
	OBJ_PASSTHROUGH
};

struct COLOR
{
    GLfloat r;
    GLfloat g;
    GLfloat b;
    GLfloat a;
};

struct VERTEX
{
    GLfloat x;
    GLfloat y;
    GLfloat z;
};

struct OBJECT
{
    enum OBJECT_TYPE Type;
    int              NumVertices;
    struct COLOR     Color;
    struct VERTEX*   Vertex;
    int              TimeStamp;
    unsigned int     Tag;
};


 
static void initData       ( void );
static void setCurrTag     ( void );
static int  parseColor     ( struct COLOR* color, const float* buff );
static int  parseVertex    ( struct VERTEX* vertex, const float* buff );
static void parseFeedback  ( const float* buff, long n );
static void displayObjects ( struct OBJECT* obj, int objNum );
static void outputObjects  ( struct OBJECT* obj, int objNum );
static int  objectComp     ( const void* A, const void* B );
static void processData    ( void );

static void psgl_derive_coords ( void );
static void psgl_derive_bbox   ( void );
static void psgl_emit_prologue ( void );
static void psgl_emit_epilogue ( void );
static void psgl_sort          ( void );





static void initData( void )
{
    int int_buffer[4];
    int i = 0;

    if( psOpenGL_Buffer != NULL )
        free( psOpenGL_Buffer );
    psOpenGL_Buffer = (GLfloat*) malloc( psOpenGL_BufferSize * sizeof( GLfloat ) );
    if( psOpenGL_Buffer == NULL )
	{
	    fprintf( stderr, "psOpenGL: insufficient memory allocating feedback buffer;\n" );
	    fprintf( stderr, "          decrease psOpenGL_BufferSize\n" );
	    exit( 1 );
	}
    assert( psOpenGL_Buffer != NULL );

    if( psOpenGL_Object != NULL )
	{
	    for( i = 0; i < psOpenGL_ObjectNum; i++ )
            if( psOpenGL_Object[i].Vertex != NULL )
                free( psOpenGL_Object[i].Vertex );
	    free( psOpenGL_Object );
	}
    psOpenGL_ObjectNum = 0;
    psOpenGL_Object = (struct OBJECT*) malloc( psOpenGL_ObjectSize * sizeof( struct OBJECT ) );
    if( psOpenGL_Object == NULL )
	{
	    fprintf( stderr, "psOpenGL: insufficient memory allocating object buffer;\n" );
	    fprintf( stderr, "          decrease psOpenGL_ObjectSize\n" );
	    exit( 1 );
	}
    assert( psOpenGL_Object );


    
    glGetIntegerv( GL_VIEWPORT, (GLint *)int_buffer );
    
    psOpenGL_OriginX = int_buffer[0];
    psOpenGL_OriginY = int_buffer[1];
    psOpenGL_SizeX   = int_buffer[2];
    psOpenGL_SizeY   = int_buffer[3];  
  
}



static void setCurrTag( void )
{
    GLboolean bool_buffer;

    psOpenGL_Tag = 0;

    glGetBooleanv( GL_DEPTH_TEST, &bool_buffer );
    if( bool_buffer )
        psOpenGL_Tag |= TAG_DepthTest;
}





void psOpenGL_open( const char* filename )
{

#ifdef PSOPENGL
    psOpenGL_Enabled = TRUE;
#endif /* PSOPENGL */
  
    if( !psOpenGL_Enabled )
        return;

    glGetIntegerv( GL_RENDER_MODE, &psOpenGL_RenderMode );
    if( psOpenGL_RenderMode == GL_FEEDBACK )
        return;
  
  
    initData();
    setCurrTag();
  
    glFeedbackBuffer( psOpenGL_BufferSize, GL_3D_COLOR, psOpenGL_Buffer );
    (void) glRenderMode( GL_FEEDBACK );
  
    psOpenGL_File = fopen( filename, "w" );
    strcpy( psOpenGL_FileName, filename );
}


void psOpenGL_close( void )
{
    if( !psOpenGL_Enabled )
        return;
  
    psOpenGL_Enabled = FALSE;

    if( psOpenGL_RenderMode == GL_FEEDBACK )
        return;

    processData();

    (void) glRenderMode( psOpenGL_RenderMode );
    /* displayObjects( psOpenGL_Object, psOpenGL_ObjectNum );*/


    psgl_derive_coords();
    psgl_derive_bbox();

    psgl_sort(); 

    psgl_emit_prologue();
    outputObjects( psOpenGL_Object, psOpenGL_ObjectNum );
    psgl_emit_epilogue();



    printf( "psOpenGL.stat: feedback buffer usage: %d of %d floats\n", psOpenGL_BufferMaxUsage, psOpenGL_BufferSize );
    printf( "psOpenGL.stat: object buffer usage: %d of %d objects\n", psOpenGL_ObjectNum, psOpenGL_ObjectSize );
    printf( "\n" );
  
}






static void processData( void )
{
    int size = 0;

    glFinish();

    size = glRenderMode( GL_FEEDBACK );
    if( size < 0 )
	{
	    fprintf( stderr, "psOpenGL: feed back buffer overflow;\n" );
	    fprintf( stderr, "          increase psOpenGL_BufferSize\n" );
	    exit( 1 );
	}
    if( size > psOpenGL_BufferMaxUsage )
        psOpenGL_BufferMaxUsage = size;
    parseFeedback( psOpenGL_Buffer, size );
}





static int parseColor( struct COLOR* color, const float* buff )
{
    int k = 0;
    color->r = buff[k++];
    color->g = buff[k++];
    color->b = buff[k++];
    color->a = buff[k++];
    if (color->r>1.) color->r=1.;
    if (color->g>1.) color->g=1.;
    if (color->b>1.) color->b=1.;
    if (color->a>1.) color->a=1.;
    assert( 0. <= color->r && color->r <= 1. );
    assert( 0. <= color->g && color->g <= 1. );
    assert( 0. <= color->b && color->b <= 1. );
    assert( 0. <= color->a && color->a <= 1. );

    return k;
}

static int parseVertex( struct VERTEX* vertex, const float* buff )
{
    int k = 0;
    vertex->x = buff[k++];
    vertex->y = buff[k++];
    vertex->z = buff[k++];
    return k;
}


static void parseFeedback( const float* buff, long n )
{
    int k = 0; /* current index */
    int i = 0; /* temp variable */

    GLfloat floatType = 0.;
    int numVertices = 0;


    enum OBJECT_TYPE objType = OBJ_ERROR;
    struct COLOR color;
    struct VERTEX vertex;
    struct OBJECT obj;

    k = 0;

    while( k < n )
	{
	    while( psOpenGL_ObjectNum >= psOpenGL_ObjectSize )
		{
		    /* attempt to automatically increase object buffer size */
		    psOpenGL_ObjectSize += OBJECTSIZE_INC;
		    psOpenGL_Object = (OBJECT *) realloc( psOpenGL_Object, psOpenGL_ObjectSize * sizeof( struct OBJECT ) );

		    if( psOpenGL_Object == NULL )
			{
			    fprintf( stderr, "psOpenGL: insufficient memory reallocating object buffer;\n" );
			    fprintf( stderr, "          reduce output complexity\n" );
			    exit( 1 );
			}
		}



	    floatType = buff[k++];

	    if( floatType == GL_POINT_TOKEN ||
            floatType == GL_LINE_TOKEN ||
            floatType == GL_LINE_RESET_TOKEN || 
            floatType == GL_POLYGON_TOKEN )
		{
		    if( floatType == GL_POINT_TOKEN )
			{
			    numVertices = 1;
			    objType = OBJ_POINT;
			}
		    else if( floatType == GL_LINE_TOKEN ||
                     floatType == GL_LINE_RESET_TOKEN )
			{
			    numVertices = 2;
			    objType = OBJ_LINE;
			}
		    else if( floatType == GL_POLYGON_TOKEN )
			{
			    numVertices = (int) buff[k++];
			    objType = OBJ_POLYGON;
			}
		    else
                assert( 0 );

		    assert( numVertices > 0 );
		    obj.Type = objType;
		    obj.NumVertices = numVertices;
		    obj.Vertex = (struct VERTEX*) malloc( numVertices * sizeof( struct VERTEX ) );
		    if( obj.Vertex == NULL )
			{
			    fprintf( stderr, "psOpenGL: insufficient memory;\n" );
			    fprintf( stderr, "          decrease psOpenGL_BufferSize or psOpenGL_ObjectSize\n" );
			    exit( 1 );
			}
		    assert( obj.Vertex != NULL );

		    obj.Color.r = 0.;
		    obj.Color.g = 0.;
		    obj.Color.b = 0.;
		    obj.Color.a = 0.;
      
		    for( i = 0; i < numVertices; i++ )
			{
			    color.r = 0.;
			    color.g = 0.;
			    color.b = 0.;
			    color.a = 0.;
   
			    k += parseVertex( &obj.Vertex[i], &buff[k] );
			    k += parseColor( &color, &buff[k] );
   
			    obj.Color.r += color.r;
			    obj.Color.g += color.g;
			    obj.Color.b += color.b;
			    obj.Color.a += color.a;
			}
      
		    obj.Color.r /= numVertices;
		    obj.Color.g /= numVertices;
		    obj.Color.b /= numVertices;
		    obj.Color.a /= numVertices;

		    assert( 0. <= obj.Color.r && obj.Color.r <= 1. );
		    assert( 0. <= obj.Color.g && obj.Color.g <= 1. );
		    assert( 0. <= obj.Color.b && obj.Color.b <= 1. );
		    assert( 0. <= obj.Color.a && obj.Color.a <= 1. );

		    obj.TimeStamp = psOpenGL_TimeStamp++;
		    obj.Tag = psOpenGL_Tag;
		    assert( psOpenGL_ObjectNum < psOpenGL_ObjectSize );
		    psOpenGL_Object[ psOpenGL_ObjectNum ] = obj;
		}
	    else if( floatType == GL_BITMAP_TOKEN )
		{
		    objType = OBJ_BITMAP;
		    k += parseColor( &color, &buff[k] );
		    k += parseVertex( &vertex, &buff[k] );
		}
	    else if( floatType == GL_DRAW_PIXEL_TOKEN ||
                 floatType == GL_COPY_PIXEL_TOKEN )
		{
		    objType = OBJ_PIXEL;
		    k += parseColor( &color, &buff[k] );
		    k += parseVertex( &vertex, &buff[k] );
		}
	    else if( floatType == GL_PASS_THROUGH_TOKEN )
		{
		    objType = OBJ_PASSTHROUGH;
		    k++;
		}
	    else 
            assert( 0 );

    
	    psOpenGL_ObjectNum++;

	}
}




static void displayObjects( struct OBJECT* obj, int objNum )
{
    int i = 0;
    int j = 0;
    for( i = 0; i < objNum; i++ )
	{
	    switch( obj[i].Type )
		{
            case OBJ_POINT:
                printf( "POINT" );
                break;
            case OBJ_LINE:
                printf( "LINE" );
                break;
            case OBJ_POLYGON:
                printf( "POLYGON" );
                break;
            case OBJ_BITMAP:
                printf( "BITMAP" );
                break;
            case OBJ_PIXEL:
                printf( "PIXEL" );
                break;
            case OBJ_PASSTHROUGH:
                printf( "PASSTHROUGH" );
                break;
            default:
                assert( 0 );
		}

	    printf( "\n" );
	    printf( "  color: %4.2f %4.2f %4.2f\n", obj[i].Color.r, obj[i].Color.g, obj[i].Color.b );
	    printf( "  vertex:\n" );
	    for( j = 0; j < obj[i].NumVertices; j++ )
		{
		    printf( "    %4.2f %4.2f %4.2f\n", obj[i].Vertex[j].x, obj[i].Vertex[j].y, obj[i].Vertex[j].z );
		}
	}
}









/* postscript ops */
#define PPI 72.0
#define XMARGIN (0.25 * PPI)
#define YMARGIN (0.25 * PPI)
#define PAGE_WIDTH ((8.5 * PPI) - 2 * XMARGIN)
#define PAGE_HEIGHT ((11.0 * PPI) - 2 * YMARGIN)

/* map coords [0..psgl_xsize, 0..psgl_ysize] to [XMARGIN..XMARGIN+PAGE_WIDTH, YMARGIN..YMARGIN+PAGE_HEIGHT] */

static   float psgl_sx, psgl_sy, psgl_xmargin, psgl_ymargin;
static   float    psgl_xmin, psgl_ymin,
    psgl_xmax, psgl_ymax;       /* bounding box */

#define PSGLX(x) (psgl_xmargin + (psgl_sx * x))
#define PSGLY(y) (psgl_ymargin + (psgl_sy * y))



static void psgl_derive_coords ( void )
{
    float hw;
    float page_width, page_height;

    psgl_xmargin = XMARGIN;
    psgl_ymargin = YMARGIN;
    page_width = PAGE_WIDTH;
    page_height = PAGE_HEIGHT;
  
    /* an explanation of coordinates used.  */
    /* all cache entries are stored in (floating point) window coords:
       [0.0..psgl_xsize] x [0.0..psgl_ysize] x [psgl_zbnear .. psgl_zbfar]   */
    /* also, widths, heights, and depths are in these coordinates.  */
  
    /* we map this rectangle to a postscript page, preserving aspect ratio */
    /* the postscript page is in points, 72/inch; 8 1/2" wide x 11" tall */
  
    /* find the height/width ratio of current window */
    hw = (float) psOpenGL_SizeY / (float) psOpenGL_SizeX;
  
    if (hw > page_height/page_width)
	{
	    /* scale height to fit */
	    psgl_sy = page_height / (float) psOpenGL_SizeY;
	    psgl_sx = psgl_sy;
	}
    else    
	{
	    psgl_sx = page_width / (float) psOpenGL_SizeX;
	    psgl_sy = psgl_sx;
	}
}


static void psgl_derive_bbox( void )
{

    float psgl_w = 0.;
    float psgl_h = 0.;

    /* clip bounding box to graphics window */
    psgl_xmin = 0.0;
    psgl_ymin = 0.0;
    psgl_xmax = psOpenGL_SizeX;
    psgl_ymax = psOpenGL_SizeY;
  
    /* convert to postscript coords */
    psgl_xmin = PSGLX (psgl_xmin);
    psgl_ymin = PSGLY (psgl_ymin);
  
    psgl_xmax = PSGLX (psgl_xmax);
    psgl_ymax = PSGLY (psgl_ymax);
  
    /* expand by one percent for border objects */
    psgl_w = psgl_xmax - psgl_xmin;
    psgl_h = psgl_ymax - psgl_ymin;
  
    psgl_xmin -= psgl_w * 0.0125;
    psgl_xmax += psgl_w * 0.0125;
  
    psgl_ymin -= psgl_h * 0.0125;
    psgl_ymax += psgl_h * 0.0125;
}



static char moveto[] = "moveto";
static char lineto[] = "lineto";
static char rlineto[] = "rlineto";
static char newpath[] = "newpath";
static char closepath[] = "closepath";
static char stroke[] = "stroke";
static char fill_[] = "fill";
static char setrgbcolor[] = "setrgbcolor";

static void psgl_emit_prologue( void )
{
  
    printf( "psOpenGL.psgl_emit_prologue(): writing '%s'\n", psOpenGL_FileName);


    /* write EPS header */
    fprintf (psOpenGL_File, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf (psOpenGL_File, "%%%%Title: (%s)\n", psOpenGL_FileName);
    fprintf (psOpenGL_File, "%%%%Creator: psOpenGL (sgi OpenGL to postscript conversion)\n");
    fprintf (psOpenGL_File, "%%%%BoundingBox: %4.2f %4.2f %4.2f %4.2f\n", psgl_xmin, psgl_ymin, psgl_xmax, psgl_ymax);
    fprintf (psOpenGL_File, "%%%%DocumentFonts: Times-Roman Times-Italic Times-Bold Courier-Bold\n");
    fprintf (psOpenGL_File, "%%%%Pages: 1\n");
  
    fprintf (psOpenGL_File, "/mt {moveto} def\n");       /* define mt == moveto */
    fprintf (psOpenGL_File, "/lt {lineto} def\n");       /* define lt == lineto */
    fprintf (psOpenGL_File, "/rl {rlineto} def\n");      /* define rl == rlineto */
    fprintf (psOpenGL_File, "/np {newpath} def\n");      /* define np == newpath */
    fprintf (psOpenGL_File, "/cp {closepath} def\n");   /* define cp == closepath */
    fprintf (psOpenGL_File, "/st {stroke} def\n");       /* define st == stroke */
    //#ifndef PSOPENGL_GRAYSCALE
    if (!psOpenGL_fill_gray)
        fprintf (psOpenGL_File, "/fi {fill} def\n");      /* define fi == fill */
    //#else
    else
        fprintf (psOpenGL_File, "/fi {gsave 1 setgray fill grestore 1 setlinewidth 0 setgray stroke } def\n");     /* define fi == fill */
    //#endif
    fprintf (psOpenGL_File, "/sc {setrgbcolor} def\n"); /* define sc == setrgbcolor */
  
    sprintf (moveto,"mt");
    sprintf (lineto,"lt");
    sprintf (rlineto,"rt");
    sprintf (newpath,"np");
    sprintf (closepath,"cp");
    sprintf (stroke,"st");
    sprintf (fill_,"fi");
    sprintf (setrgbcolor,"sc");
  
    fprintf (psOpenGL_File, "\ngsave\n\n");
    /* invoke defaults */
    fprintf (psOpenGL_File, "1 setlinejoin\n");          /* round ends */
    fprintf (psOpenGL_File, "1 setlinecap\n");        /* round ends */
    fprintf (psOpenGL_File, "1.0 setlinewidth\n");
    fprintf (psOpenGL_File, "/Times-Bold findfont 18 scalefont setfont\n");
  
    /* clip to bounding box */
    fprintf (psOpenGL_File, "newpath\n");
    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmin, psgl_ymin, moveto);
    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmax, psgl_ymin, lineto);
    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmax, psgl_ymax, lineto);
    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmin, psgl_ymax, lineto);
    fprintf (psOpenGL_File, "closepath clip\n");
  
    /* show the clip path */
    if (psOpenGL_DebugLevel)  
	{
	    fprintf (psOpenGL_File, "[2 4] 0 setdash 3.0 setlinewidth\n");
	    fprintf (psOpenGL_File, "newpath\n");
	    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmin, psgl_ymin, moveto);
	    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmax, psgl_ymin, lineto);
	    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmax, psgl_ymax, lineto);
	    fprintf (psOpenGL_File, "%4.2f %4.2f %s\n", psgl_xmin, psgl_ymax, lineto);
	    fprintf (psOpenGL_File, "closepath stroke\n");
	    fprintf (psOpenGL_File, "[] 0 setdash 1.0 setlinewidth\n");
	}
}



static void outputObjects( struct OBJECT* obj, int objNum )
{
    int i = 0;
    int j = 0;

    struct COLOR currColor;

    for( i = 0; i < objNum; i++ )
	{
	    if ( i == 0 || 
             obj[i].Color.r != currColor.r ||
             obj[i].Color.g != currColor.g ||
             obj[i].Color.b != currColor.b )
		{
		    currColor = obj[i].Color;
		    fprintf( psOpenGL_File, "%4.2f %4.2f %4.2f %s\n", currColor.r, currColor.g, currColor.b, setrgbcolor );
		}

	    switch( obj[i].Type )
		{
            case OBJ_POINT:
                fprintf( psOpenGL_File, "%s\n", newpath);
                fprintf( psOpenGL_File, "%4.2f %4.2f %s %4.2f %4.2f %s\n",
                         PSGLX( obj[i].Vertex[0].x ), PSGLY( obj[i].Vertex[0].y ), moveto, 
                         PSGLX( obj[i].Vertex[0].x ), PSGLY( obj[i].Vertex[0].y ), lineto );
                fprintf( psOpenGL_File, "%s %s\n", closepath, stroke );
                break;

            case OBJ_LINE:
            case OBJ_POLYGON:
                fprintf( psOpenGL_File, "%s\n", newpath );
                for( j = 0; j < obj[i].NumVertices; j++ )
                    fprintf( psOpenGL_File, "%4.2f %4.2f %s\n",
                             PSGLX( obj[i].Vertex[j].x ), PSGLY( obj[i].Vertex[j].y ), (j == 0) ? moveto : lineto );
      
                if( obj[i].Type == OBJ_POLYGON )
                    fprintf( psOpenGL_File, "%s %s\n", closepath, fill_ );
                else
                    fprintf( psOpenGL_File, "%s\n", stroke );      
                break;

            case OBJ_BITMAP:
                break;
            case OBJ_PIXEL:
                break;
            case OBJ_PASSTHROUGH:
                break;
            default:
                assert( 0 );
		}

	}
}



static void psgl_emit_epilogue(void )
{
    fprintf (psOpenGL_File, "\ngrestore\n");
  
    if( psOpenGL_Epilogue != NULL )
        fprintf (psOpenGL_File, "%s\n", psOpenGL_Epilogue);

    printf( "psOpenGL.psgl_emit_epilogue(): closing '%s'\n", psOpenGL_FileName);
    fclose( psOpenGL_File );
    psOpenGL_File = NULL;

}





static void psgl_sort( void )
{
    printf( "psOpenGL.psgl_sort(): sorting %d items... ", psOpenGL_ObjectNum );

    /* sort the cache by z depth, time drawn, etc. */
    qsort( psOpenGL_Object, psOpenGL_ObjectNum, sizeof( struct OBJECT ), objectComp );

    printf( "done.\n" );
}

////typedef Vec3T<double> CVec;

#define DEPTH_EPS 0.01

static int compLinePoly( struct OBJECT* l, struct OBJECT* p) {
    /* use first three vertices of the polygon; 
       if they are collinear, give up */
  
    if ( p->NumVertices >= 3 ) {
        Vec3T<float> p0(p->Vertex[0].x, p->Vertex[0].y, p->Vertex[0].z);
        Vec3T<float> p1(p->Vertex[1].x, p->Vertex[1].y, p->Vertex[1].z);
        Vec3T<float> p2(p->Vertex[2].x, p->Vertex[2].y, p->Vertex[2].z);
        assert( l->NumVertices == 2);
        Vec3T<float> l0(l->Vertex[0].x, l->Vertex[0].y, l->Vertex[0].z);
        Vec3T<float> l1(l->Vertex[1].x, l->Vertex[1].y, l->Vertex[1].z);
        Vec3T<float> px(p1 - p0);
        Vec3T<float> py(p2 - p0);
        Vec3T<float> n = cross(px,py).dir();    //n.normalize();
        Vec3T<float> dpl0(l0 - p0);
        Vec3T<float> dpl1(l1 - p0);
        /* are we in the plane? then use time stamp */
        if( (fabs(dot(n,dpl0)) < DEPTH_EPS) && 
            (fabs(dot(n,dpl1)) < DEPTH_EPS)) {
            //      cerr << "in one plane";
            // cerr << "n = " << n << endl;
            // if( l->TimeStamp < p->TimeStamp )
            return -1;
            //else
            // return 1;
        } 
        // cerr << "not in one plane ";
        double det = px.x()*py.y() - px.y()*py.x();
        double t0,t1;
        if( fabs(det) > DEPTH_EPS) { 
            double idet = 1.0/det;
            double s00 = (py.y()*dpl0.x() - py.x()*dpl0.y())*idet;
            double s01 = -(px.x()*dpl0.y() - px.y()*dpl0.x())*idet;      
            t0 = px.z()*s00 + py.z()*s01 - dpl0.z();
            double s10 = (py.y()*dpl1.x() - py.x()*dpl1.y())*idet;
            double s11 = -(px.x()*dpl1.y() - px.y()*dpl1.x())*idet;      
            t1 = px.z()*s10 + py.z()*s11 - dpl1.z();
        } else { 
            cerr << "det = 0\n";
            return -1;
        }
        //    cerr << t0 << " " << t1 << " " << det;
        //cerr << "\nvecs: " << p0 << "\t" << p1 << "\t" << p2 << "\t\t" << l0 << 
        //  "\t" << l1 << endl;
        // cerr << "n = " << n << endl;
        // favor the lines
        // both ends of the line are in front
        if( (t0 >= -DEPTH_EPS) && (t1 >= -DEPTH_EPS) ) return -1;
        // both ends of the line are behind
        if( (t0 < -DEPTH_EPS) && (t1 < -DEPTH_EPS) ) return 1;
        // in other cases
        // should split it in two, but that would require restructuring 
        // of the code; for the time being compromize
        if( t0 + t1 > 0 ) return -1; else return 1;
    } else { 
        printf( "poly with 2 vertices\n"); return -1;
    }
}

/* if zbuffering was off either or both of A and B, later one wins */
/* otherwise, compare z depths: */
/* A farther than B: -1 */
/* A closer than B:   1 */
/* farthest items end up at beginning of list, closest at end. */
static int objectComp( const void* A, const void* B )
{
    struct OBJECT* a = (struct OBJECT*) A;
    struct OBJECT* b = (struct OBJECT*) B;

    int j = 0;
    int aDepthTest = 0;
    int bDepthTest = 0;
    GLfloat aDepth = 0.;
    GLfloat bDepth = 0.;

    aDepthTest = a->Tag & TAG_DepthTest;
    bDepthTest = b->Tag & TAG_DepthTest;

    if( aDepthTest && bDepthTest ) {
        if( ( a->Type == OBJ_LINE )   && ( b->Type == OBJ_POLYGON )) {
            return compLinePoly(a,b);

        } else if( ( b->Type == OBJ_LINE )   && ( a->Type == OBJ_POLYGON )) {

            return -compLinePoly(b,a);
      
        } else  {
            for( j = 0; j < a->NumVertices; j++ )
                aDepth += a->Vertex[j].z;
            aDepth /= a->NumVertices;
    
            for( j = 0; j < b->NumVertices; j++ )
                bDepth += b->Vertex[j].z;
            bDepth /= b->NumVertices;
        }
    }
    // give line priority, whenever at least one coord is in front 
    // of the polygon

    if( !(aDepthTest && bDepthTest) || (aDepth == bDepth) )
	{
	    if( a->TimeStamp < b->TimeStamp )
            return -1; /* draw A frist */
	    else
            return 1; /* draw B first */
	}
    else
	{
	    if( aDepth > bDepth )
            return -1; /* draw A frist */
	    else
            return 1; /* draw B first */
	}
}

/* REPLACEMENT OpenGL function calls */


void psOpenGL_End( void )
{
    glEnd();

    if( !psOpenGL_Enabled )
        return;

    processData();
}
