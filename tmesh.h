/* tmesh.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, January 2005
 * Version 2.6.4+
 */

typedef struct {
        int wx, wy, wz;
        float nx, ny, nz;
        float potential;
        int normal;
        Vert v;
    } VrtStruct;

typedef struct {
        int u, v, w;
        int normal;
        int col;
    } TriStruct;


#define SurfNone         0x00
#define SurfDots         0x01
#define SurfLines        0x02
#define SurfVectors      0x03
#define SurfNormal       0x04
#define SurfSolid        0x05
#define SurfTranslucent  0x06
#define SurfTransparent  0x07

/* FACES = 20*(4^ORDER)   */
/*       = 20<<(ORDER<<1) */
/* COUNT = FACES/2        */

#define NormOrder  3
#define NormFaces  1280
#define NormCount  640

/* ORDER  FACES  VERTS */
/*   0      20     12  */
/*   1      80     42  */
/*   2     320    162  */
/*   3    1280    642  */

#define ORDER   3
#define FACES   1280
#define VERTS   642
#define EDGES   20




#ifdef TMESH
int SwapBytes;

unsigned char NormInten[NormFaces];
Real Normal[NormCount][3];

/* Vertex co-ordinates of a unit icosahedron */
Real TV[VERTS][3] = {
    {  0.00000000, -0.85065080, -0.52573110 },  /*  0   9 */
    { -0.52573110,  0.00000000, -0.85065080 },  /*  1  10 */
    { -0.85065080, -0.52573110,  0.00000000 },  /*  2  11 */
    {  0.00000000, -0.85065080,  0.52573110 },  /*  3   6 */
    {  0.52573110,  0.00000000, -0.85065080 },  /*  4   7 */
    { -0.85065080,  0.52573110,  0.00000000 },  /*  5   8 */
    {  0.00000000,  0.85065080, -0.52573110 },  /*  6   3 */
    { -0.52573110,  0.00000000,  0.85065080 },  /*  7   4 */
    {  0.85065080, -0.52573110,  0.00000000 },  /*  8   5 */
    {  0.00000000,  0.85065080,  0.52573110 },  /*  9   0 */
    {  0.52573110,  0.00000000,  0.85065080 },  /* 10   1 */
    {  0.85065080,  0.52573110,  0.00000000 }   /* 11   2 */
        };

Real TN[FACES][3];
int TF[FACES][3];
Real TC[EDGES];
Real TS[EDGES];
unsigned int TEdges;
unsigned int TVert;
unsigned int TFace;

int SurfaceMode;

int HavePotentials;

VrtStruct *Vrt;
TriStruct *Tri;
int TriAlloc;
int TriCount;
int VrtAlloc;
int VrtCount;

#else
extern int SwapBytes;
extern unsigned char NormInten[NormFaces];
extern Real Normal[NormCount][3];
extern Real TV[VERTS][3];
extern Real TN[FACES][3];
extern int TF[FACES][3];
extern Real TC[EDGES];
extern Real TS[EDGES];
extern unsigned int TEdges;
extern unsigned int TVert;
extern unsigned int TFace;

extern int SurfaceMode;

extern int HavePotentials;

extern VrtStruct *Vrt;
extern TriStruct *Tri;
extern int TriCount;
extern int TriAlloc;
extern int VrtCount;
extern int VrtAlloc;
extern int TrnCount;
#endif


void InitialiseCAD( int );
void InitialiseTMesh( void );
void DeleteSurface( void );
void DisplaySurface( void );
void ColourSurfaceAttrib( int, int, int );
void ColourSurfacePotential( void );

void LoadGraspSurface( FILE*, int );
void CreateRibbons( void );

