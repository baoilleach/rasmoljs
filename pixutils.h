/* pixutils.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

typedef struct {
        int px, py, pz;      /* Spline Control Co-ordinate */
        int tx, ty, tz;      /* Spline Direction Vector    */
        int hnx, hny, hnz;   /* Horizontal Normal Vector   */
        int vnx, vny, vnz;   /* Vertical Normal Vector     */
        int dx, dy, dz;      /* Ribbon Height Vector       */
        int wx, wy, wz;      /* Ribbon Width Vector        */
        char hinten;         /* Horizontal Intensity       */
        char vinten;         /* Vertical Intensity         */
        short hsize;         /* Horizontal Vector Length   */
        short vsize;         /* Vertical Vector Length     */
        short wide;          /* Ribbon Width               */
    } Knot;

typedef struct {
        Pixel __huge *fbuf;
        short __huge *dbuf;
        int xmax, ymax;
        int yskip;
    } ViewStruct;


#define MAXVERT  6
typedef struct {
        int x, y, z;
        int inten;
    } Vert;

typedef struct {
        Vert v[MAXVERT];
        int count;
    } Poly;


#define ZValid(z)     ((!UseSlabPlane) || ((z)<SlabValue))
#define XValid(x)     (((x)>=0)&&((x)<View.xmax))
#define YValid(y)     (((y)>=0)&&((y)<View.ymax))


#ifdef PIXUTILS
ViewStruct View;
int SplineCount;
int FontSize;

#else
extern ViewStruct View;
extern int SplineCount;
extern int FontSize;
#endif


void PlotPoint( int, int, int, int );
void ClipPoint( int, int, int, int );
void PlotDeepPoint( int, int, int, int );
void ClipDeepPoint( int, int, int, int );

void ClipLine( int, int, int, int, int, int, int );
void DrawTwinLine( int, int, int, int, int, int, int, int );
void ClipTwinLine( int, int, int, int, int, int, int, int );
void DrawTwinVector( int, int, int, int, int, int, int, int );
void ClipTwinVector( int, int, int, int, int, int, int, int );
void DrawTwinVector2( int, int, int, int, int, int, int, int );
void ClipTwinVector2( int, int, int, int, int, int, int, int );
void ClipDashVector( int, int, int, int, int, int, int, int );

void DrawSphere( int, int, int, int, int );
void ClipSphere( int, int, int, int, int );

void DrawCylinder( int, int, int, int, int, int, int, int, int );
void ClipCylinder( int, int, int, int, int, int, int, int, int );

void DrawFlatTriangle( Vert *p, Vert *q, Vert *r, int col );
void DrawTriangle( Vert *p, Vert *q, Vert *r );
void ClipTriangle( Vert *p, Vert *q, Vert *r );

void DrawPolygonOutline( Poly* );
void ClipPolygonOutline( Poly* );
void DrawFlatPolygon( Poly* );
void ClipFlatPolygon( Poly* );
void DrawPolygon( Poly* );
void ClipPolygon( Poly* );

void DashRibbon( Knot __far*, Knot __far*, int, int );
void StrandRibbon( Knot __far*, Knot __far*, int, int );
void SolidRibbon2( Knot __far*, Knot __far*, int, int );
void SolidRibbon( Knot __far*, Knot __far*, int );
void RectRibbon( Knot __far*, Knot __far*, int, int, int );

void SetFontSize( int );
void DisplayString( int, int, int, char*, int );
void InitialisePixUtils( void );

