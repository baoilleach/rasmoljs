/* tmesh.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, January 2005
 * Version 2.6.4+
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define TMESH
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "cmndline.h"
#include "scripts.h"
#include "repres.h"
#include "render.h"
#include "pixutils.h"
#include "graphics.h"
#include "tmesh.h"

/*=======================*/
/*  Normal Quantization  */
/*=======================*/ 

#define RootSix       2.44948974278

/* Face list of a unit icosahedron */
static const int F[20][4] = {
    {  0,  1,  2 },  /*   0   19 */
    {  0,  1,  4 },  /*   1   18 */
    {  0,  2,  3 },  /*   2   17 */
    {  0,  3,  8 },  /*   3   16 */
    {  0,  4,  8 },  /*   4   15 */
    {  1,  2,  5 },  /*   5   14 */
    {  1,  4,  6 },  /*   6   13 */
    {  1,  5,  6 },  /*   7   12 */
    {  2,  3,  7 },  /*   8   11 */
    {  2,  5,  7 },  /*   9   10 */

    {  4,  8, 11 },  /*  10    9 */
    {  4,  6, 11 },  /*  11    8 */
    {  3,  8, 10 },  /*  12    7 */
    {  3,  7, 10 },  /*  13    6 */
    {  8, 10, 11 },  /*  14    5 */
    {  5,  7,  9 },  /*  15    4 */
    {  5,  6,  9 },  /*  16    3 */
    {  6,  9, 11 },  /*  17    2 */
    {  7,  9, 10 },  /*  18    1 */
    {  9, 10, 11 }   /*  19    0 */
            };

static int NrmCount;

static void Normalize( Real *v )
{
    register Real len;

    len = (Real)sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
    v[0] /= len;  v[1] /= len;  v[2] /= len;
}


static void Tesselate( const Real *p, const Real *q, const Real *r, int d )
{
    auto Real u[3], v[3], w[3];
    register Real *n;

    if( d-- ) {
        u[0] = p[0]+q[0];  v[0] = q[0]+r[0];  w[0] = r[0]+p[0];
        u[1] = p[1]+q[1];  v[1] = q[1]+r[1];  w[1] = r[1]+p[1];
        u[2] = p[2]+q[2];  v[2] = q[2]+r[2];  w[2] = r[2]+p[2];
        Normalize(u);
        Normalize(v);
        Normalize(w);

        Tesselate( u, v, w, d );
        Tesselate( p, u, w, d );
        Tesselate( u, q, v, d );
        Tesselate( w, v, r, d );
    } else {
        n = Normal[NrmCount++];
        n[0] = p[0]+q[0]+r[0];
        n[1] = p[1]+q[1]+r[1];
        n[2] = p[2]+q[2]+r[2];
        Normalize( n );
        n[0] *= (ColourDepth/RootSix);
        n[1] *= (ColourDepth/RootSix);
        n[2] *= (ColourDepth/RootSix);
    }
}


static void TesselateNormals( void )
{
    register int i;

    NrmCount = 0;
    for( i=0; i<10; i++ )
        Tesselate( TV[F[i][0]], TV[F[i][1]], TV[F[i][2]], NormOrder );
}


static int QuantizeNormal( Real vx, Real vy, Real vz )
{
    register Real *n;
    register Real p,bestp;
    register int i,d,res;
    register int flag;
    register Real len;

    /* Avoid compiler warnings! */
    flag = False;
    res = 0;

    /* Normalise input vector! */
    len = sqrt( vx*vx + vy*vy + vz*vz );
    vx /= len;  vy /= len;  vz /= len;

    bestp = 0.0;
    for( i=0; i<NormCount; i+=(1<<(NormOrder<<1)) ) {
        n = Normal[i];
        p = vx*n[0] + vy*n[1] + vz*n[2];

        if( p > bestp ) {
            flag = False;
            bestp = p;
            res = i;
        } else if( -p > bestp ) {
            flag = True;
            bestp = -p;
            res = i;
        }
    }

    if( flag ) {
        vx = -vx;
        vy = -vy;
        vz = -vz;
    }

    for( d = 1<<((NormOrder-1)<<1); d ; d>>=2 ) {
        i = res+d;  n = Normal[i];
        p = vx*n[0] + vy*n[1] + vz*n[2];
        if( p > bestp ) {
            bestp = p;
            res = i;
        }

        i += d;  n = Normal[i];
        p = vx*n[0] + vy*n[1] + vz*n[2];
        if( p > bestp ) {
            bestp = p;
            res = i;
        }

        i += d;  n = Normal[i];
        p = vx*n[0] + vy*n[1] + vz*n[2];
        if( p > bestp ) {
            bestp = p;
            res = i;
        }
    }

#ifdef DEBUG
    printf("%lg %lg %lg -> %lg %lg %lg (%lg)\n",vx,vy,vz,
           Normal[res][0], Normal[res][1], Normal[res][2], bestp );
#endif

    if( flag ) {
        return (res<<1)+1;
    } else return res<<1;
}


static void FatalTMeshError( char *ptr )
{
    char buffer[80];

    sprintf(buffer,"TMesh Error: Unable to allocate %s!",ptr);
    RasMolFatalExit(buffer);
}


/*========================================*/
/*  Triangle Display List Initialization  */
/*========================================*/

static void DetermineByteOrder( void )
{
    register char *ptr;
    static int test;

    test = 0x00000001;
    ptr = (char*)&test;
    SwapBytes = ptr[0];
}


void InitialiseTMesh( void )
{
  DetermineByteOrder();

  VrtCount = VrtAlloc = 0;
  TriCount = TriAlloc = 0;
  Vrt = (VrtStruct*)0;
  Tri = (TriStruct*)0;

  SurfaceMode = SurfSolid;
  HavePotentials = False;
  TesselateNormals();
}


void DeleteSurface( void )
{
    register int shade;
    register int i;

    for (i=0; i<TriCount; i++)
    {
      shade = Colour2Shade(Tri[i].col);
      Shade[shade].refcount--;
    }

    if( Vrt ) free(Vrt);
    Vrt = (VrtStruct*)0;
    if( Tri ) free(Tri);
    Tri = (TriStruct*)0;
    VrtCount = VrtAlloc = 0;
    TriCount = TriAlloc = 0;
    HavePotentials = False;
}


static int NewVertex( void )
{
    if( VrtCount == VrtAlloc )
    {
        if( VrtAlloc == 0 )
        {
            Vrt = (VrtStruct*)malloc(1024*sizeof(VrtStruct));
            VrtAlloc = 1024;
        } else {
            VrtAlloc += 1024;
            Vrt = (VrtStruct*)realloc(Vrt,VrtAlloc*sizeof(VrtStruct));
        }
        if( !Vrt ) FatalTMeshError("vertex");
    }
    return VrtCount++;
}


static int NewTriangle( void )
{
    if( TriCount == TriAlloc )
    {
        if( TriAlloc == 0 )
        {
            Tri = (TriStruct*)malloc(1024*sizeof(TriStruct));
            TriAlloc = 1024;
        } else {
            TriAlloc += 1024;
            Tri = (TriStruct*)realloc(Tri,TriAlloc*sizeof(TriStruct));
        }
        if( !Tri ) FatalTMeshError("triangle");
    }
    return TriCount++;
}


static void TransformVertices( void )
{
    register VrtStruct *v;
    register int xi,yi,zi;
    register Real x,y,z;
    register int i;

    for( i=0; i<VrtCount; i++ ) {
        v = &Vrt[i];
        x = v->wx;
        y = v->wy;
        z = v->wz;

        xi = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
        yi = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
        zi = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
        v->v.inten = ColourMask;
        v->v.x = xi;
        v->v.y = yi;
        v->v.z = zi;
    }
}


static void QuantizeInten( void )
{
    register unsigned char *ptr;
    register Real lx,ly,lz;
    register int i,inten;
    register Real *n;

#ifdef INVERT
    lx = RotX[0] + RotY[0] + (RotZ[0] + RotZ[0]);
    ly = RotX[1] + RotY[1] + (RotZ[1] + RotZ[1]);
    lz = RotX[2] + RotY[2] + (RotZ[2] + RotZ[2]);
#else
    lx = RotX[0] - RotY[0] + (RotZ[0] + RotZ[0]);
    ly = RotX[1] - RotY[1] + (RotZ[1] + RotZ[1]);
    lz = RotX[2] - RotY[2] + (RotZ[2] + RotZ[2]);
#endif

    ptr = NormInten;
    for( i=0; i<NormCount; i++ ) {
        n = Normal[i];
        inten = (int)(lx*n[0]+ly*n[1]+lz*n[2]);
        *ptr++ = (inten >  1)? (unsigned char)inten    : 1;
        *ptr++ = (inten < -1)? (unsigned char)(-inten) : 1;
    }
}


/*=======================*/
/*  Rendering Functions  */
/*=======================*/

static void DisplaySurfDots()
{
    register Vert *p,*q,*r;
    register TriStruct *t;
    register int i, col;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u].v;
        q = &Vrt[t->v].v;
        r = &Vrt[t->w].v;

        col = t->col;
        ClipDeepPoint(p->x, p->y, p->z, col);
        ClipDeepPoint(q->x, q->y, q->z, col);
        ClipDeepPoint(r->x, r->y, r->z, col);
    }
}


static void DisplaySurfLines()
{
    register Vert *p,*q,*r;
    register TriStruct *t;
    register int i, col;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u].v;
        q = &Vrt[t->v].v;
        r = &Vrt[t->w].v;

        col = t->col;
        ClipTwinVector( p->x,p->y,p->z, q->x,q->y,q->z, col,col);
        ClipTwinVector( q->x,q->y,q->z, r->x,r->y,r->z, col,col);
        ClipTwinVector( r->x,r->y,r->z, p->x,p->y,p->z, col,col);
    }
}


static void DisplaySurfVectors( void )
{
    register VrtStruct *p,*q,*r;
    register TriStruct *t;
    register int i, col;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u];
        q = &Vrt[t->v];
        r = &Vrt[t->w];

        if( (q->v.x-p->v.x)*(r->v.y-p->v.y) <=
            (r->v.x-p->v.x)*(q->v.y-p->v.y) )
            continue;

        col = t->col + (NormInten[p->normal]+NormInten[q->normal]+1)/2;
        ClipTwinLine( p->v.x,p->v.y,p->v.z, q->v.x,q->v.y,q->v.z, col,col);
        col = t->col + (NormInten[q->normal]+NormInten[r->normal]+1)/2;
        ClipTwinLine( q->v.x,q->v.y,q->v.z, r->v.x,r->v.y,r->v.z, col,col);
        col = t->col + (NormInten[r->normal]+NormInten[p->normal]+1)/2;
        ClipTwinLine( r->v.x,r->v.y,r->v.z, p->v.x,p->v.y,p->v.z, col,col);
    }
}


static void DisplaySurfNormal( void )
{
    register Vert *p,*q,*r;
    register TriStruct *t;
    register int i,inten;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u].v;
        q = &Vrt[t->v].v;
        r = &Vrt[t->w].v;

        if( (q->x-p->x)*(r->y-p->y) <= (r->x-p->x)*(q->y-p->y) )
            continue;
        inten = t->col+NormInten[t->normal];
        /* ClipFlatTriangle( p, q, r, inten ); */

        p->inten = inten;
        q->inten = inten;
        r->inten = inten;
        ClipTriangle( p, q, r );
    }
}


static void DisplaySurfSolid( void )
{
    register VrtStruct *p, *q, *r;
    register TriStruct *t;
    register int i;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u];
        q = &Vrt[t->v];
        r = &Vrt[t->w];

        if( (q->v.x-p->v.x)*(r->v.y-p->v.y) <=
            (r->v.x-p->v.x)*(q->v.y-p->v.y) )
            continue;

        p->v.inten = t->col+NormInten[p->normal];
        q->v.inten = t->col+NormInten[q->normal];
        r->v.inten = t->col+NormInten[r->normal];
        ClipTriangle( &p->v, &q->v, &r->v );
    }
}


void DisplaySurface( void )
{
    TransformVertices();
    QuantizeInten();

    switch (SurfaceMode) {
        case SurfDots:         DisplaySurfDots();
                               break;

        case SurfLines:        DisplaySurfLines();
                               break;

        case SurfVectors:      DisplaySurfVectors();
                               break;

        case SurfNormal:       DisplaySurfNormal();
                               break;

        case SurfSolid:        DisplaySurfSolid();
                               break;

        case SurfTransparent:  DisplaySurfSolid();
                               StencilBuffers(True);
                               break;

        case SurfTranslucent:  DisplaySurfSolid();
                               StencilBuffers(False);
                               break;
    }
}


/*===============================*/
/*  Surface Colouring Functions  */
/*===============================*/

void ColourSurfaceAttrib( int r, int g, int b )
{
    register int shade;
    register int col;
    register int i;

    for (i=0; i<TriCount; i++)
    {
      shade = Colour2Shade(Tri[i].col);
      Shade[shade].refcount--;
    }

    shade = DefineShade(r,g,b);
    col = Shade2Colour(shade);
    for (i=0; i<TriCount; i++)
    {
      Shade[shade].refcount++;
      Tri[i].col = col;
    }
}


void ColourSurfacePotential( void )
{
    register int s0, s1, s2, s3, s4, s5, s6;
    register double maxpot, minpot;
    register double median, pot;
    register TriStruct *t;
    register int shade;
    register int i,val;

    if( TriCount == 0 )
        return;

    if( !HavePotentials ) {
        ColourSurfaceAttrib(255,255,255);
        return;
    }

    maxpot = Vrt[0].potential;
    minpot = Vrt[0].potential;
    for( i=1; i<VrtCount; i++ ) {
        if( Vrt[i].potential > maxpot )
            maxpot = Vrt[i].potential;
        if( Vrt[i].potential < minpot )
            minpot = Vrt[i].potential;
    }

    if( minpot == maxpot ) {
        ColourSurfaceAttrib(255,255,255);
        HavePotentials = False;
        return;
    }

    for( i=0; i<TriCount; i++ )
    {
      shade = Colour2Shade(Tri[i].col);
      Shade[shade].refcount--;
    }

    /* Define in priority order! */
    s3 = DefineShade(255,255,255);  Shade[s3].refcount++;
    s0 = DefineShade(255,  0,  0);  Shade[s0].refcount++;
    s6 = DefineShade(  0,  0,255);  Shade[s6].refcount++;
    s2 = DefineShade(255,170,170);  Shade[s2].refcount++;
    s4 = DefineShade(170,170,255);  Shade[s4].refcount++;
    s1 = DefineShade(255, 85, 85);  Shade[s1].refcount++;
    s5 = DefineShade( 85, 85,255);  Shade[s5].refcount++;

    median = 0.5*(minpot+maxpot);
    maxpot -= median;
    minpot -= median;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        pot = (Vrt[t->u].potential +
               Vrt[t->v].potential +
               Vrt[t->w].potential -
               median);
        if( pot > 0.0 ) {
            val = (int)floor(3.5*pot/maxpot);
            switch (val)
            {
            case 0:  t->col = Shade2Colour(s3);
                     Shade[s3].refcount++;
                     break;
            case 1:  t->col = Shade2Colour(s4);
                     Shade[s4].refcount++;
                     break;
            case 2:  t->col = Shade2Colour(s5);
                     Shade[s5].refcount++;
                     break;
            default:
            case 3:  t->col = Shade2Colour(s6);
                     Shade[s6].refcount++;
                     break;
            }
        } else {
            val = (int)floor(3.5*pot/minpot);
            switch (val)
            {
            case 0:  t->col = Shade2Colour(s3);
                     Shade[s3].refcount++;
                     break;
            case 1:  t->col = Shade2Colour(s2);
                     Shade[s2].refcount++;
                     break;
            case 2:  t->col = Shade2Colour(s1);
                     Shade[s1].refcount++;
                     break;
            default:
            case 3:  t->col = Shade2Colour(s0);
                     Shade[s0].refcount++;
                     break;
            }
        }
    }

    Shade[s0].refcount--;
    Shade[s1].refcount--;
    Shade[s2].refcount--;
    Shade[s3].refcount--;
    Shade[s4].refcount--;
    Shade[s5].refcount--;
    Shade[s6].refcount--;
}


/*=============================*/
/*  Display List Manipulation  */
/*=============================*/

void TriangleNormal1( TriStruct *t )
{
    register VrtStruct *p, *q, *r;
    register Real nx, ny, nz;

    p = &Vrt[t->u];
    q = &Vrt[t->v];
    r = &Vrt[t->w];

    nx = p->nx + q->nx + r->nx;
    ny = p->ny + q->ny + r->ny;
    nz = p->nz + q->nz + r->nz;

    t->normal = QuantizeNormal( nx, ny, nz );
}


void TriangleNormals1( void )
{
    register Real nx,ny,nz;
    register VrtStruct *p,*q,*r;
    register TriStruct *t;
    register int i;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u];
        q = &Vrt[t->v];
        r = &Vrt[t->w];

        nx = p->nx + q->nx + r->nx;
        ny = p->ny + q->ny + r->ny;
        nz = p->nz + q->nz + r->nz;

	t->normal = QuantizeNormal( nx, ny, nz );
    }
}


void TriangleNormal2( TriStruct *t )
{
    register VrtStruct *p, *q, *r;
    register Real nx, ny, nz;

    p = &Vrt[t->u];
    q = &Vrt[t->v];
    r = &Vrt[t->w];

    nx = (q->wy-p->wy)*(r->wz-p->wz) - (r->wy-p->wy)*(q->wz-p->wz);
    ny = (r->wx-p->wx)*(q->wz-p->wz) - (q->wx-p->wx)*(r->wz-p->wz);
    nz = (q->wx-p->wx)*(r->wy-p->wy) - (r->wx-p->wx)*(q->wy-p->wy);

    t->normal = QuantizeNormal( nx, ny, nz );
}


void TriangleNormals2( void )
{
    register Real nx,ny,nz;
    register VrtStruct *p,*q,*r;
    register TriStruct *t;
    register int i;

    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
        p = &Vrt[t->u];
        q = &Vrt[t->v];
        r = &Vrt[t->w];

        nx = (q->wy-p->wy)*(r->wz-p->wz) - (r->wy-p->wy)*(q->wz-p->wz);
        ny = (r->wx-p->wx)*(q->wz-p->wz) - (q->wx-p->wx)*(r->wz-p->wz);
        nz = (q->wx-p->wx)*(r->wy-p->wy) - (r->wx-p->wx)*(q->wy-p->wy);

	t->normal = QuantizeNormal( nx, ny, nz );
    }
}


void NormalizeTriangle( TriStruct *t )
{
    register VrtStruct *p, *q, *r;
    register int normal;
    Real n[3];

    p = &Vrt[t->u];
    q = &Vrt[t->v];
    r = &Vrt[t->w];

    n[0] = (q->wy-p->wy)*(r->wz-p->wz) - (r->wy-p->wy)*(q->wz-p->wz);
    n[1] = (r->wx-p->wx)*(q->wz-p->wz) - (q->wx-p->wx)*(r->wz-p->wz);
    n[2] = (q->wx-p->wx)*(r->wy-p->wy) - (r->wx-p->wx)*(q->wy-p->wy);
    Normalize(n);

    normal = QuantizeNormal( n[0], n[1], n[2] );
    t->normal = normal;

    p->nx = (float)n[0];
    p->ny = (float)n[1];
    p->nz = (float)n[2];
    p->normal = normal;

    q->nx = (float)n[0];
    q->ny = (float)n[1];
    q->nz = (float)n[2];
    q->normal = normal;

    r->nx = (float)n[0];
    r->ny = (float)n[1];
    r->nz = (float)n[2];
    r->normal = normal;
}


/*=====================*/
/*  Grasp Surface I/O  */
/*=====================*/

static void ReadGraspLine( FILE *fp, char *buffer )
{
    register char *ptr;
    register int i;

    for( i=0; i<4; i++ )
        if( getc(fp) == EOF )
            return;

    ptr = buffer;
    for( i=0; i<80; i++ )
        *ptr++ = fgetc(fp);
    *ptr = '\0';

    for( i=0; i<4; i++ )
        if( getc(fp) == EOF )
            return;
}



static float ReadGraspFloat( FILE *fp )
{
    register char *ptr;
    register char ch;
    static float result;

    ptr = (char*)&result;
    fread(&result,sizeof(float),1,fp);

    if( SwapBytes )
    {   ch = ptr[0];  ptr[0] = ptr[3];  ptr[3] = ch;
        ch = ptr[1];  ptr[1] = ptr[2];  ptr[2] = ch;
    }
    return result;
}


static int ReadGraspInt( FILE *fp )
{
    register char *ptr;
    register char ch;
    static int result;

    ptr = (char*)&result;
    fread(&result,sizeof(int),1,fp);
    if( SwapBytes )
    {   ch = ptr[0];  ptr[0] = ptr[3];  ptr[3] = ch;
        ch = ptr[1];  ptr[1] = ptr[2];  ptr[2] = ch;
    }
    return result;
}


static void BeginGraspBlock( FILE *fp )
{
    static int result;
    fread(&result,sizeof(int),1,fp);
}


static void EndGraspBlock( FILE *fp )
{
    static int result;
    fread(&result,sizeof(int),1,fp);
}


static int GraspAttributeIndex( char *name, char *buffer )
{
    register int index;
    register char *src;
    register char *dst;
    char temp[64];

    index = 0;
    src = buffer;
    while( *src && (*src!=' ') )
    {   if( *src == ',') src++; 

        dst = temp;
        while( *src && (*src!=',') && (*src!=' ') )
            *dst++ = *src++;
        *dst = '\0';

        index++;
        if( !strcmp(temp,name) )
            return index;
    }
    return 0;
}


void LoadGraspSurface( FILE *fp, int info )
{
    register unsigned long size;
    register int ipot,iacc;
    register int shade,col;
    register float x,y,z;
    register int u,v,w;
    register int i,j;
    char buffer[82];

    DeleteSurface();

    /* Assume "format=2" */
    ReadGraspLine(fp,buffer);

    ReadGraspLine(fp,buffer);
    iacc = GraspAttributeIndex("accessibles",buffer);

    ReadGraspLine(fp,buffer);
    ipot = GraspAttributeIndex("potentials",buffer);

    ReadGraspLine(fp,buffer);
    sscanf(buffer,"%d %d",&VrtCount,&TriCount);
    ReadGraspLine(fp,buffer);

    /* AllocateTriMesh */
    size = VrtCount*sizeof(VrtStruct);
    Vrt = (VrtStruct*)malloc(size);
    if( !Vrt ) {
        VrtCount = 0;
        TriCount = 0;
        return;
    }
    
    size = TriCount*sizeof(TriStruct);
    Tri = (TriStruct*)malloc(size);
    if( !Tri ) {
        free(Vrt);
        VrtCount = 0;
        TriCount = 0;
        return;
    }
    TriAlloc = TriCount;
    VrtAlloc = VrtCount;

    ReDrawFlag |= RFColour | RFRefresh;
    shade = DefineShade(255,255,255);
    col = Shade2Colour(shade);

    /* Vertices */
    BeginGraspBlock(fp);
    for( i=0; i<VrtCount; i++ ) {
        x =  ReadGraspFloat(fp);
        y =  ReadGraspFloat(fp);
        z = -ReadGraspFloat(fp);

#ifdef INVERT
        y = -y;
#endif

        Vrt[i].wx = (int)(x * 250.0) - OrigCX;
        Vrt[i].wy = (int)(y * 250.0) - OrigCY;
        Vrt[i].wz = (int)(z * 250.0) - OrigCZ;
    }
    EndGraspBlock(fp);

    if( iacc ) {
        /* Accessibles */
        BeginGraspBlock(fp);
        for( i=0; i<VrtCount; i++ ) {
            ReadGraspFloat(fp);
            ReadGraspFloat(fp);
            ReadGraspFloat(fp);
        }
        EndGraspBlock(fp);
    }

    /* Normals */
    BeginGraspBlock(fp);
    for( i=0; i<VrtCount; i++ ) {
        x =  ReadGraspFloat(fp);
        y =  ReadGraspFloat(fp);
        z = -ReadGraspFloat(fp);

#ifdef INVERT
        y = -y;
#endif

        Vrt[i].nx = x;
        Vrt[i].ny = y;
        Vrt[i].nz = z;
        Vrt[i].normal = QuantizeNormal( x, y, z );
    }
    EndGraspBlock(fp);

    /* Triangles */
    BeginGraspBlock(fp);
    for( i=0; i<TriCount; i++ ) {
        u = ReadGraspInt(fp)-1;
        v = ReadGraspInt(fp)-1;
        w = ReadGraspInt(fp)-1;
#ifdef INVERT
        Tri[i].u = u;
        Tri[i].v = v;
        Tri[i].w = w;
#else
        Tri[i].u = u;
        Tri[i].v = w;
        Tri[i].w = v;
#endif
        Shade[shade].refcount++;
        Tri[i].col = col;
    }
    EndGraspBlock(fp);

    if( ipot ) {
        /* Skip Preceeding Properties! */
        for( j=1; j<ipot; j++ ) {
            BeginGraspBlock(fp);
            for( i=0; i<VrtCount; i++ )
                ReadGraspFloat(fp);
            EndGraspBlock(fp);
        }

        BeginGraspBlock(fp);
        for( i=0; i<VrtCount; i++ )
             Vrt[i].potential = ReadGraspFloat(fp);
        EndGraspBlock(fp);
        HavePotentials = True;
        
    } else {
        for( i=0; i<VrtCount; i++ )
            Vrt[i].potential = 0.0;
        HavePotentials = False;
    }

    TriangleNormals1();
    SurfaceMode = SurfSolid;

    if( info )
    {   InvalidateCmndLine();

        sprintf(buffer,"%d triangles, %d vertices\n",TriCount,VrtCount);
        WriteString(buffer);
    }
}

/*===================*/
/*  Protein Ribbons  */
/*===================*/

static void CreateSphere( Real *c, double rad, int col )
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register unsigned int i;
    register int idx,shade;
    register int u, v, w;
    register int *f;
    int tmp[VERTS];

    for( i=0; i<TVert; i++ ) {
      idx = NewVertex();
      vrt = &Vrt[idx];
      nx = TV[i][0];
      ny = TV[i][1];
      nz = TV[i][2];

      vrt->potential = 0.0;
      vrt->normal = QuantizeNormal(nx,ny,nz);
      vrt->wx = (int)(nx*rad + c[0]);
      vrt->wy = (int)(ny*rad + c[1]);
      vrt->wz = (int)(nz*rad + c[2]);
      vrt->nx = (float)nx;
      vrt->ny = (float)ny;
      vrt->nz = (float)nz;
      tmp[i] = idx;
    }

    for( i=0; i<TFace; i++ )
    {
      f = TF[i];
      idx = NewTriangle();
      t = &Tri[idx];
      u = tmp[f[0]];
      v = tmp[f[1]];
      w = tmp[f[2]];

      vrt = &Vrt[u];
      nx = vrt->nx;
      ny = vrt->ny;
      nz = vrt->nz;

      vrt = &Vrt[v];
      nx += vrt->nx;
      ny += vrt->ny;
      nz += vrt->nz;

      vrt = &Vrt[w];
      nx += vrt->nx;
      ny += vrt->ny;
      nz += vrt->nz;

      t->normal = QuantizeNormal(nx,ny,nz);
      t->col = col;
      t->u = u;
      t->v = v;
      t->w = w;
    }

    shade = Colour2Shade(col);
    Shade[shade].refcount += TFace;
}

#ifdef UNUSED
static void CreateCylinder( Real *src, Real *dst, double rad,
                            int cap1, int cap2, int col )
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register Real rx, ry, rz;
    register int idx, curr, prev;
    register int normal, shade;
    register unsigned int i;
    register Real tc, ts;
    register int cen;

    Real d[3], x[3], y[3];
    Real top[EDGES][3];
    Real bot[EDGES][3];
    int tidx[EDGES];
    int bidx[EDGES];

    d[0] = dst[0] - src[0];
    d[1] = dst[1] - src[1];
    d[2] = dst[2] - src[2];
    Normalize(d);

    if( fabs(d[0]) < fabs(d[1]) ) {
        x[0] = 0.0;
        x[1] = d[2];
        x[2] = -d[1];
    } else {
        x[0] = d[2];
        x[1] = 0.0;
        x[2] = -d[0];
    }
    Normalize(x);

    y[0] = d[1]*x[2] - d[2]*x[1];
    y[1] = d[2]*x[0] - d[0]*x[2];
    y[2] = d[0]*x[1] - d[1]*x[0];
    Normalize(y);

    for( i=0; i<TEdges; i++ ) {
        nx = TC[i]*x[0] + TS[i]*y[0];
        ny = TC[i]*x[1] + TS[i]*y[1];
        nz = TC[i]*x[2] + TS[i]*y[2];
        normal = QuantizeNormal(nx,ny,nz);
  
        rx = rad*nx;
        ry = rad*ny;
        rz = rad*nz;

        bot[i][0] = src[0]+rx;
        bot[i][1] = src[1]+ry;
        bot[i][2] = src[2]+rz;

        top[i][0] = dst[0]+rx;
        top[i][1] = dst[1]+ry;
        top[i][2] = dst[2]+rz;

        idx = NewVertex();
        vrt = &Vrt[idx];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)bot[i][0];
        vrt->wy = (int)bot[i][1];
        vrt->wz = (int)bot[i][2];
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;
        bidx[i] = idx;

        idx = NewVertex();
        vrt = &Vrt[idx];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)top[i][0];
        vrt->wy = (int)top[i][1];
        vrt->wz = (int)top[i][2];
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;
        tidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
        tc = TC[curr]+TC[prev];
        ts = TS[curr]+TS[prev];

        nx = tc*x[0] + ts*y[0];
        ny = tc*x[1] + ts*y[1];
        nz = tc*x[2] + ts*y[2];

        normal = QuantizeNormal(nx,ny,nz);

        idx = NewTriangle();
        t = &Tri[idx];
        t->normal = normal;
        t->col = col;
        t->u = bidx[curr];
        t->v = tidx[curr];
        t->w = tidx[prev];

        idx = NewTriangle();
        t = &Tri[idx];
        t->normal = normal;
        t->col = col;
        t->u = bidx[curr];
        t->v = tidx[prev];
        t->w = bidx[prev];

        prev = curr++;
    }

    shade = Colour2Shade(col);
    Shade[shade].refcount += TEdges*2;

    if (cap1) {
        nx = -d[0];
        ny = -d[1];
        nz = -d[2];
        normal = QuantizeNormal(nx,ny,nz);
  
        cen = NewVertex();
        vrt = &Vrt[cen];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)src[0];
        vrt->wy = (int)src[1];
        vrt->wz = (int)src[2];
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;

        for( i=0; i<TEdges; i++ ) {
            idx = NewVertex();
            vrt = &Vrt[idx];
            vrt->potential = 0.0;
            vrt->normal = normal;
            vrt->wx = (int)bot[i][0];
            vrt->wy = (int)bot[i][1];
            vrt->wz = (int)bot[i][2];
            vrt->nx = (float)nx;
            vrt->ny = (float)ny;
            vrt->nz = (float)nz;
            bidx[i] = idx;
        }

        curr = 0;
        prev = TEdges-1;
        for( i=0; i<TEdges; i++ ) {
            idx = NewTriangle();
            t = &Tri[idx];
            t->normal = normal;
            t->col = col;
            t->u = cen;
            t->v = bidx[curr];
            t->w = bidx[prev];
            prev = curr++;
        }
        Shade[shade].refcount += TEdges;
    }

    if (cap2) {
        nx = d[0];
        ny = d[1];
        nz = d[2];
        normal = QuantizeNormal(nx,ny,nz);
  
        cen = NewVertex();
        vrt = &Vrt[cen];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)dst[0];
        vrt->wy = (int)dst[1];
        vrt->wz = (int)dst[2];
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;

        for( i=0; i<TEdges; i++ ) {
            idx = NewVertex();
            vrt = &Vrt[idx];
            vrt->potential = 0.0;
            vrt->normal = normal;
            vrt->wx = (int)top[i][0];
            vrt->wy = (int)top[i][1];
            vrt->wz = (int)top[i][2];
            vrt->nx = (float)nx;
            vrt->ny = (float)ny;
            vrt->nz = (float)nz;
            tidx[i] = idx;
        }

        curr = 0;
        prev = TEdges-1;
        for( i=0; i<TEdges; i++ ) {
            idx = NewTriangle();
            t = &Tri[idx];
            t->normal = normal;
            t->col = col;
            t->u = cen;
            t->v = tidx[prev];
            t->w = tidx[curr];
            prev = curr++;
        }
        Shade[shade].refcount += TEdges;
    }
}
#endif


typedef struct {
      Real p[3];  /* Spline Control Co-ordinate */
      Real t[3];  /* Spline Control Vector      */
      Real h[3];  /* Horizontal Normal Vector   */
      Real v[3];  /* Vertical Normal Vector     */
      Real wide;  /* Ribbon Width               */
    } TKnot;


static void Interpolate(TKnot *src, TKnot *dst, TKnot *mid,
                        double p, double q)
{
    mid->t[0] = p*src->t[0] + q*dst->t[0];
    mid->t[1] = p*src->t[1] + q*dst->t[1];
    mid->t[2] = p*src->t[2] + q*dst->t[2];
    Normalize(mid->t);

    mid->h[0] = p*src->h[0] + q*dst->h[0];
    mid->h[1] = p*src->h[1] + q*dst->h[1];
    mid->h[2] = p*src->h[2] + q*dst->h[2];
    Normalize(mid->h);

    mid->v[0] = p*src->v[0] + q*dst->v[0];
    mid->v[1] = p*src->v[1] + q*dst->v[1];
    mid->v[2] = p*src->v[2] + q*dst->v[2];
    Normalize(mid->v);

    mid->wide = p*src->wide + q*dst->wide;
}


static void TraceBegCap(TKnot *ptr, int wide, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register int normal, shade;
    register int idx, curr, prev;
    register unsigned int i;
    register Real x, y, z;
    register int cen;
    int cidx[EDGES];

    nx = -ptr->t[0];
    ny = -ptr->t[1];
    nz = -ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    cen = NewVertex();
    vrt = &Vrt[cen];
    vrt->potential = 0.0;
    vrt->normal = normal;
    vrt->wx = (int)ptr->p[0];
    vrt->wy = (int)ptr->p[1];
    vrt->wz = (int)ptr->p[2];
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    for( i=0; i<TEdges; i++ ) {
        x = TC[i]*ptr->h[0] + TS[i]*ptr->v[0];
        y = TC[i]*ptr->h[1] + TS[i]*ptr->v[1];
        z = TC[i]*ptr->h[2] + TS[i]*ptr->v[2];

        idx = NewVertex();
        vrt = &Vrt[idx];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)(ptr->p[0]+wide*x);
        vrt->wy = (int)(ptr->p[1]+wide*y);
        vrt->wz = (int)(ptr->p[2]+wide*z);
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;
        cidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
        idx = NewTriangle();
        t = &Tri[idx];
        t->normal = normal;
        t->col = col;
        t->u = cen;
        t->v = cidx[curr];
        t->w = cidx[prev];
        prev = curr++;
    }
    shade = Colour2Shade(col);
    Shade[shade].refcount += TEdges;
}


static void TraceEndCap(TKnot *ptr, int wide, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register int normal, shade;
    register int idx, curr, prev;
    register unsigned int i;
    register Real x, y, z;
    register int cen;
    int cidx[EDGES];

    nx = ptr->t[0];
    ny = ptr->t[1];
    nz = ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    cen = NewVertex();
    vrt = &Vrt[cen];
    vrt->potential = 0.0;
    vrt->normal = normal;
    vrt->wx = (int)ptr->p[0];
    vrt->wy = (int)ptr->p[1];
    vrt->wz = (int)ptr->p[2];
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    for( i=0; i<TEdges; i++ ) {
        x = TC[i]*ptr->h[0] + TS[i]*ptr->v[0];
        y = TC[i]*ptr->h[1] + TS[i]*ptr->v[1];
        z = TC[i]*ptr->h[2] + TS[i]*ptr->v[2];

        idx = NewVertex();
        vrt = &Vrt[idx];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)(ptr->p[0]+wide*x);
        vrt->wy = (int)(ptr->p[1]+wide*y);
        vrt->wz = (int)(ptr->p[2]+wide*z);
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;
        cidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
        idx = NewTriangle();
        t = &Tri[idx];
        t->normal = normal;
        t->col = col;
        t->u = cen;
        t->v = cidx[prev];
        t->w = cidx[curr];
        prev = curr++;
    }
    shade = Colour2Shade(col);
    Shade[shade].refcount += TEdges;
}


static void CreateTrace(TKnot *src, TKnot *dst,
                        int wide, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register int idx, curr, prev;
    register int normal, shade;
    register unsigned int i;

    Real top[EDGES][3];
    Real bot[EDGES][3];
    int tidx[EDGES];
    int bidx[EDGES];

    for( i=0; i<TEdges; i++ ) {
      nx = TC[i]*src->h[0] + TS[i]*src->v[0];
      ny = TC[i]*src->h[1] + TS[i]*src->v[1];
      nz = TC[i]*src->h[2] + TS[i]*src->v[2];
      normal = QuantizeNormal(nx,ny,nz);

      bot[i][0] = src->p[0]+wide*nx;
      bot[i][1] = src->p[1]+wide*ny;
      bot[i][2] = src->p[2]+wide*nz;

      idx = NewVertex();
      vrt = &Vrt[idx];
      vrt->potential = 0.0;
      vrt->normal = normal;
      vrt->wx = (int)bot[i][0];
      vrt->wy = (int)bot[i][1];
      vrt->wz = (int)bot[i][2];
      vrt->nx = (float)nx;
      vrt->ny = (float)ny;
      vrt->nz = (float)nz;
      bidx[i] = idx;

      nx = TC[i]*dst->h[0] + TS[i]*dst->v[0];
      ny = TC[i]*dst->h[1] + TS[i]*dst->v[1];
      nz = TC[i]*dst->h[2] + TS[i]*dst->v[2];
      normal = QuantizeNormal(nx,ny,nz);

      top[i][0] = dst->p[0]+wide*nx;
      top[i][1] = dst->p[1]+wide*ny;
      top[i][2] = dst->p[2]+wide*nz;

      idx = NewVertex();
      vrt = &Vrt[idx];
      vrt->potential = 0.0;
      vrt->normal = normal;
      vrt->wx = (int)top[i][0];
      vrt->wy = (int)top[i][1];
      vrt->wz = (int)top[i][2];
      vrt->nx = (float)nx;
      vrt->ny = (float)ny;
      vrt->nz = (float)nz;
      tidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
      idx = NewTriangle();
      t = &Tri[idx];
      t->u = bidx[curr];
      t->v = tidx[curr];
      t->w = tidx[prev];
      TriangleNormal1(t);
      t->col = col;

      idx = NewTriangle();
      t = &Tri[idx];
      t->u = bidx[curr];
      t->v = tidx[prev];
      t->w = bidx[prev];
      TriangleNormal1(t);
      t->col = col;

      prev = curr++;
    }

    shade = Colour2Shade(col);
    Shade[shade].refcount += 2*TEdges;
}


static void CartoonBegCap(TKnot *ptr, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register int normal, shade;
    register int p, q, r, s;
    register Real nx,ny,nz;
    register int idx;
    Real h[3], v[3];

    nx = -ptr->t[0];
    ny = -ptr->t[1];
    nz = -ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    h[0] = ptr->h[0]*ptr->wide;
    h[1] = ptr->h[1]*ptr->wide;
    h[2] = ptr->h[2]*ptr->wide;
    v[0] = ptr->v[0]*CartoonHeight;
    v[1] = ptr->v[1]*CartoonHeight;
    v[2] = ptr->v[2]*CartoonHeight;

    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + h[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] + h[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] + h[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + h[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] + h[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] + h[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - h[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] - h[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] - h[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - h[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] - h[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] - h[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = q;
    t->w = r;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = s;

    shade = Colour2Shade(col);
    Shade[shade].refcount += 2;
}


static void CartoonEndCap(TKnot *ptr, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register int normal, shade;
    register int p, q, r, s;
    register Real nx,ny,nz;
    register int idx;
    Real h[3], v[3];
     
    nx = ptr->t[0];
    ny = ptr->t[1];
    nz = ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    h[0] = ptr->h[0]*ptr->wide;
    h[1] = ptr->h[1]*ptr->wide;
    h[2] = ptr->h[2]*ptr->wide;
    v[0] = ptr->v[0]*CartoonHeight;
    v[1] = ptr->v[1]*CartoonHeight;
    v[2] = ptr->v[2]*CartoonHeight;

    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + h[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] + h[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] + h[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + h[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] + h[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] + h[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - h[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] - h[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] - h[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - h[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] - h[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] - h[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = s;
    t->w = r;

    shade = Colour2Shade(col);
    Shade[shade].refcount += 2;
}


static void CreateCartoon(TKnot *src, TKnot *dst, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register int normal, shade;
    register int p, q, r, s;
    register int idx;

    Real sh[3], sv[3];
    Real dh[3], dv[3];
    Real crn[8][3];
    Real n[3];

    sh[0] = src->h[0]*src->wide;
    sh[1] = src->h[1]*src->wide;
    sh[2] = src->h[2]*src->wide;
    sv[0] = src->v[0]*CartoonHeight;
    sv[1] = src->v[1]*CartoonHeight;
    sv[2] = src->v[2]*CartoonHeight;

    dh[0] = dst->h[0]*dst->wide;
    dh[1] = dst->h[1]*dst->wide;
    dh[2] = dst->h[2]*dst->wide;
    dv[0] = dst->v[0]*CartoonHeight;
    dv[1] = dst->v[1]*CartoonHeight;
    dv[2] = dst->v[2]*CartoonHeight;

    crn[0][0] = src->p[0] + sh[0] + sv[0];
    crn[0][1] = src->p[1] + sh[1] + sv[1];
    crn[0][2] = src->p[2] + sh[2] + sv[2];
    crn[1][0] = src->p[0] + sh[0] - sv[0];
    crn[1][1] = src->p[1] + sh[1] - sv[1];
    crn[1][2] = src->p[2] + sh[2] - sv[2];
    crn[2][0] = src->p[0] - sh[0] - sv[0];
    crn[2][1] = src->p[1] - sh[1] - sv[1];
    crn[2][2] = src->p[2] - sh[2] - sv[2];
    crn[3][0] = src->p[0] - sh[0] + sv[0];
    crn[3][1] = src->p[1] - sh[1] + sv[1];
    crn[3][2] = src->p[2] - sh[2] + sv[2];

    crn[4][0] = dst->p[0] + dh[0] + dv[0];
    crn[4][1] = dst->p[1] + dh[1] + dv[1];
    crn[4][2] = dst->p[2] + dh[2] + dv[2];
    crn[5][0] = dst->p[0] + dh[0] - dv[0];
    crn[5][1] = dst->p[1] + dh[1] - dv[1];
    crn[5][2] = dst->p[2] + dh[2] - dv[2];
    crn[6][0] = dst->p[0] - dh[0] - dv[0];
    crn[6][1] = dst->p[1] - dh[1] - dv[1];
    crn[6][2] = dst->p[2] - dh[2] - dv[2];
    crn[7][0] = dst->p[0] - dh[0] + dv[0];
    crn[7][1] = dst->p[1] - dh[1] + dv[1];
    crn[7][2] = dst->p[2] - dh[2] + dv[2];

    Normalize(sh);
    Normalize(sv);
    Normalize(dh);
    Normalize(dv);

    /* Right */
    normal = QuantizeNormal(sh[0],sh[1],sh[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[0][0];
    vrt->wy = (int)crn[0][1];
    vrt->wz = (int)crn[0][2];
    vrt->nx = (float)sh[0];
    vrt->ny = (float)sh[1];
    vrt->nz = (float)sh[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[1][0];
    vrt->wy = (int)crn[1][1];
    vrt->wz = (int)crn[1][2];
    vrt->nx = (float)sh[0];
    vrt->ny = (float)sh[1];
    vrt->nz = (float)sh[2];

    normal = QuantizeNormal(dh[0],dh[1],dh[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[4][0];
    vrt->wy = (int)crn[4][1];
    vrt->wz = (int)crn[4][2];
    vrt->nx = (float)dh[0];
    vrt->ny = (float)dh[1];
    vrt->nz = (float)dh[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[5][0];
    vrt->wy = (int)crn[5][1];
    vrt->wz = (int)crn[5][2];
    vrt->nx = (float)dh[0];
    vrt->ny = (float)dh[1];
    vrt->nz = (float)dh[2];

    n[0] = sh[0] + dh[0];
    n[1] = sh[1] + dh[1];
    n[2] = sh[2] + dh[2];
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Bottom */
    normal = QuantizeNormal(-sv[0],-sv[1],-sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[1][0];
    vrt->wy = (int)crn[1][1];
    vrt->wz = (int)crn[1][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[2][0];
    vrt->wy = (int)crn[2][1];
    vrt->wz = (int)crn[2][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    normal = QuantizeNormal(-dv[0],-dv[1],-dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[5][0];
    vrt->wy = (int)crn[5][1];
    vrt->wz = (int)crn[5][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[6][0];
    vrt->wy = (int)crn[6][1];
    vrt->wz = (int)crn[6][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    n[0] = -(sv[0] + dv[0]);
    n[1] = -(sv[1] + dv[1]);
    n[2] = -(sv[2] + dv[2]);
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Left */
    normal = QuantizeNormal(-sh[0],-sh[1],-sh[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[2][0];
    vrt->wy = (int)crn[2][1];
    vrt->wz = (int)crn[2][2];
    vrt->nx = (float)-sh[0];
    vrt->ny = (float)-sh[1];
    vrt->nz = (float)-sh[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[3][0];
    vrt->wy = (int)crn[3][1];
    vrt->wz = (int)crn[3][2];
    vrt->nx = (float)-sh[0];
    vrt->ny = (float)-sh[1];
    vrt->nz = (float)-sh[2];

    normal = QuantizeNormal(-dh[0],-dh[1],-dh[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[6][0];
    vrt->wy = (int)crn[6][1];
    vrt->wz = (int)crn[6][2];
    vrt->nx = (float)-dh[0];
    vrt->ny = (float)-dh[1];
    vrt->nz = (float)-dh[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[7][0];
    vrt->wy = (int)crn[7][1];
    vrt->wz = (int)crn[7][2];
    vrt->nx = (float)-dh[0];
    vrt->ny = (float)-dh[1];
    vrt->nz = (float)-dh[2];

    n[0] = -(sh[0] + dh[0]);
    n[1] = -(sh[1] + dh[1]);
    n[2] = -(sh[2] + dh[2]);
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Top */
    normal = QuantizeNormal(sv[0],sv[1],sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[3][0];
    vrt->wy = (int)crn[3][1];
    vrt->wz = (int)crn[3][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[0][0];
    vrt->wy = (int)crn[0][1];
    vrt->wz = (int)crn[0][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    normal = QuantizeNormal(dv[0],dv[1],dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[7][0];
    vrt->wy = (int)crn[7][1];
    vrt->wz = (int)crn[7][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[4][0];
    vrt->wy = (int)crn[4][1];
    vrt->wz = (int)crn[4][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    n[0] = sv[0] + dv[0];
    n[1] = sv[1] + dv[1];
    n[2] = sv[2] + dv[2];
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;


    shade = Colour2Shade(col);
    Shade[shade].refcount += 8;
}


static void ArrowBegCap(TKnot *ptr, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register int normal, shade;
    register int p, q, r, s;
    register Real nx,ny,nz;
    register Real wide;
    register int idx;
    Real h[3], v[3], f[3];

    nx = -ptr->t[0];
    ny = -ptr->t[1];
    nz = -ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    wide = ptr->wide*1.5;
    v[0] = ptr->v[0]*CartoonHeight;
    v[1] = ptr->v[1]*CartoonHeight;
    v[2] = ptr->v[2]*CartoonHeight;
    h[0] = ptr->h[0]*ptr->wide;
    h[1] = ptr->h[1]*ptr->wide;
    h[2] = ptr->h[2]*ptr->wide;
    f[0] = ptr->h[0]*wide;
    f[1] = ptr->h[1]*wide;
    f[2] = ptr->h[2]*wide;

    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + f[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] + f[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] + f[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + f[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] + f[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] + f[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + h[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] + h[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] + h[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] + h[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] + h[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] + h[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = q;
    t->w = r;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = s;


    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - h[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] - h[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] - h[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - h[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] - h[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] - h[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - f[0] - v[0]);
    vrt->wy = (int)(ptr->p[1] - f[1] - v[1]);
    vrt->wz = (int)(ptr->p[2] - f[2] - v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)(ptr->p[0] - f[0] + v[0]);
    vrt->wy = (int)(ptr->p[1] - f[1] + v[1]);
    vrt->wz = (int)(ptr->p[2] - f[2] + v[2]);
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = q;
    t->w = r;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = s;

    shade = Colour2Shade(col);
    Shade[shade].refcount += 4;
}


static void CreateArrow(TKnot *src, TKnot *dst, int col)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register int normal, shade;
    register int p, q, r, s;
    register Real wide;
    register int idx;

    Real sh[3], sv[3], sf[3];
    Real dh[3], dv[3], df[3];
    Real crn[16][3];
    Real n[3];

    wide = src->wide*1.5;
    sv[0] = src->v[0]*CartoonHeight;
    sv[1] = src->v[1]*CartoonHeight;
    sv[2] = src->v[2]*CartoonHeight;
    sh[0] = src->h[0]*src->wide;
    sh[1] = src->h[1]*src->wide;
    sh[2] = src->h[2]*src->wide;
    sf[0] = src->h[0]*wide;
    sf[1] = src->h[1]*wide;
    sf[2] = src->h[2]*wide;

    wide = dst->wide*1.5;
    dv[0] = dst->v[0]*CartoonHeight;
    dv[1] = dst->v[1]*CartoonHeight;
    dv[2] = dst->v[2]*CartoonHeight;
    dh[0] = dst->h[0]*dst->wide;
    dh[1] = dst->h[1]*dst->wide;
    dh[2] = dst->h[2]*dst->wide;
    df[0] = dst->h[0]*wide;
    df[1] = dst->h[1]*wide;
    df[2] = dst->h[2]*wide;

    crn[0][0] = src->p[0] + sh[0] + sv[0];
    crn[0][1] = src->p[1] + sh[1] + sv[1];
    crn[0][2] = src->p[2] + sh[2] + sv[2];
    crn[1][0] = src->p[0] + sf[0] + sv[0];
    crn[1][1] = src->p[1] + sf[1] + sv[1];
    crn[1][2] = src->p[2] + sf[2] + sv[2];
    crn[2][0] = src->p[0] + sf[0] - sv[0];
    crn[2][1] = src->p[1] + sf[1] - sv[1];
    crn[2][2] = src->p[2] + sf[2] - sv[2];
    crn[3][0] = src->p[0] + sh[0] - sv[0];
    crn[3][1] = src->p[1] + sh[1] - sv[1];
    crn[3][2] = src->p[2] + sh[2] - sv[2];
    crn[4][0] = src->p[0] - sh[0] - sv[0];
    crn[4][1] = src->p[1] - sh[1] - sv[1];
    crn[4][2] = src->p[2] - sh[2] - sv[2];
    crn[5][0] = src->p[0] - sf[0] - sv[0];
    crn[5][1] = src->p[1] - sf[1] - sv[1];
    crn[5][2] = src->p[2] - sf[2] - sv[2];
    crn[6][0] = src->p[0] - sf[0] + sv[0];
    crn[6][1] = src->p[1] - sf[1] + sv[1];
    crn[6][2] = src->p[2] - sf[2] + sv[2];
    crn[7][0] = src->p[0] - sh[0] + sv[0];
    crn[7][1] = src->p[1] - sh[1] + sv[1];
    crn[7][2] = src->p[2] - sh[2] + sv[2];

    crn[8][0] = dst->p[0] + dh[0] + dv[0];
    crn[8][1] = dst->p[1] + dh[1] + dv[1];
    crn[8][2] = dst->p[2] + dh[2] + dv[2];
    crn[9][0] = dst->p[0] + df[0] + dv[0];
    crn[9][1] = dst->p[1] + df[1] + dv[1];
    crn[9][2] = dst->p[2] + df[2] + dv[2];
    crn[10][0] = dst->p[0] + df[0] - dv[0];
    crn[10][1] = dst->p[1] + df[1] - dv[1];
    crn[10][2] = dst->p[2] + df[2] - dv[2];
    crn[11][0] = dst->p[0] + dh[0] - dv[0];
    crn[11][1] = dst->p[1] + dh[1] - dv[1];
    crn[11][2] = dst->p[2] + dh[2] - dv[2];
    crn[12][0] = dst->p[0] - dh[0] - dv[0];
    crn[12][1] = dst->p[1] - dh[1] - dv[1];
    crn[12][2] = dst->p[2] - dh[2] - dv[2];
    crn[13][0] = dst->p[0] - df[0] - dv[0];
    crn[13][1] = dst->p[1] - df[1] - dv[1];
    crn[13][2] = dst->p[2] - df[2] - dv[2];
    crn[14][0] = dst->p[0] - df[0] + dv[0];
    crn[14][1] = dst->p[1] - df[1] + dv[1];
    crn[14][2] = dst->p[2] - df[2] + dv[2];
    crn[15][0] = dst->p[0] - dh[0] + dv[0];
    crn[15][1] = dst->p[1] - dh[1] + dv[1];
    crn[15][2] = dst->p[2] - dh[2] + dv[2];

    Normalize(sh);
    Normalize(sv);
    Normalize(dh);
    Normalize(dv);

    /* Top Right */
    normal = QuantizeNormal(sv[0],sv[1],sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[0][0];
    vrt->wy = (int)crn[0][1];
    vrt->wz = (int)crn[0][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[1][0];
    vrt->wy = (int)crn[1][1];
    vrt->wz = (int)crn[1][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    normal = QuantizeNormal(dv[0],dv[1],dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[8][0];
    vrt->wy = (int)crn[8][1];
    vrt->wz = (int)crn[8][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[9][0];
    vrt->wy = (int)crn[9][1];
    vrt->wz = (int)crn[9][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    n[0] = sv[0] + dv[0];
    n[1] = sv[1] + dv[1];
    n[2] = sv[2] + dv[2];
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Right */
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[1][0];
    vrt->wy = (int)crn[1][1];
    vrt->wz = (int)crn[1][2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[2][0];
    vrt->wy = (int)crn[2][1];
    vrt->wz = (int)crn[2][2];

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[9][0];
    vrt->wy = (int)crn[9][1];
    vrt->wz = (int)crn[9][2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[10][0];
    vrt->wy = (int)crn[10][1];
    vrt->wz = (int)crn[10][2];

    idx = NewTriangle();
    t = &Tri[idx];
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;
    NormalizeTriangle(t);

    idx = NewTriangle();
    t = &Tri[idx];
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;
    NormalizeTriangle(t);

    /* Bottom Right */
    normal = QuantizeNormal(-sv[0],-sv[1],-sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[2][0];
    vrt->wy = (int)crn[2][1];
    vrt->wz = (int)crn[2][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[3][0];
    vrt->wy = (int)crn[3][1];
    vrt->wz = (int)crn[3][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    normal = QuantizeNormal(-dv[0],-dv[1],-dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[10][0];
    vrt->wy = (int)crn[10][1];
    vrt->wz = (int)crn[10][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[11][0];
    vrt->wy = (int)crn[11][1];
    vrt->wz = (int)crn[11][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    n[0] = -(sv[0] + dv[0]);
    n[1] = -(sv[1] + dv[1]);
    n[2] = -(sv[2] + dv[2]);
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Bottom */
    normal = QuantizeNormal(-sv[0],-sv[1],-sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[3][0];
    vrt->wy = (int)crn[3][1];
    vrt->wz = (int)crn[3][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[4][0];
    vrt->wy = (int)crn[4][1];
    vrt->wz = (int)crn[4][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    normal = QuantizeNormal(-dv[0],-dv[1],-dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[11][0];
    vrt->wy = (int)crn[11][1];
    vrt->wz = (int)crn[11][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[12][0];
    vrt->wy = (int)crn[12][1];
    vrt->wz = (int)crn[12][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    n[0] = -(sv[0] + dv[0]);
    n[1] = -(sv[1] + dv[1]);
    n[2] = -(sv[2] + dv[2]);
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Bottom Left */
    normal = QuantizeNormal(-sv[0],-sv[1],-sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[4][0];
    vrt->wy = (int)crn[4][1];
    vrt->wz = (int)crn[4][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[5][0];
    vrt->wy = (int)crn[5][1];
    vrt->wz = (int)crn[5][2];
    vrt->nx = (float)-sv[0];
    vrt->ny = (float)-sv[1];
    vrt->nz = (float)-sv[2];

    normal = QuantizeNormal(-dv[0],-dv[1],-dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[12][0];
    vrt->wy = (int)crn[12][1];
    vrt->wz = (int)crn[12][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[13][0];
    vrt->wy = (int)crn[13][1];
    vrt->wz = (int)crn[13][2];
    vrt->nx = (float)-dv[0];
    vrt->ny = (float)-dv[1];
    vrt->nz = (float)-dv[2];

    n[0] = -(sv[0] + dv[0]);
    n[1] = -(sv[1] + dv[1]);
    n[2] = -(sv[2] + dv[2]);
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Left */
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[5][0];
    vrt->wy = (int)crn[5][1];
    vrt->wz = (int)crn[5][2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[6][0];
    vrt->wy = (int)crn[6][1];
    vrt->wz = (int)crn[6][2];

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[13][0];
    vrt->wy = (int)crn[13][1];
    vrt->wz = (int)crn[13][2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->potential = 0.0;
    vrt->wx = (int)crn[14][0];
    vrt->wy = (int)crn[14][1];
    vrt->wz = (int)crn[14][2];

    idx = NewTriangle();
    t = &Tri[idx];
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;
    NormalizeTriangle(t);

    idx = NewTriangle();
    t = &Tri[idx];
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;
    NormalizeTriangle(t);

    /* Top Left */
    normal = QuantizeNormal(sv[0],sv[1],sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[6][0];
    vrt->wy = (int)crn[6][1];
    vrt->wz = (int)crn[6][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[7][0];
    vrt->wy = (int)crn[7][1];
    vrt->wz = (int)crn[7][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    normal = QuantizeNormal(dv[0],dv[1],dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[14][0];
    vrt->wy = (int)crn[14][1];
    vrt->wz = (int)crn[14][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[15][0];
    vrt->wy = (int)crn[15][1];
    vrt->wz = (int)crn[15][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    n[0] = sv[0] + dv[0];
    n[1] = sv[1] + dv[1];
    n[2] = sv[2] + dv[2];
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;

    /* Top */
    normal = QuantizeNormal(sv[0],sv[1],sv[2]);
 
    p = NewVertex();
    vrt = &Vrt[p];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[7][0];
    vrt->wy = (int)crn[7][1];
    vrt->wz = (int)crn[7][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    q = NewVertex();
    vrt = &Vrt[q];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[0][0];
    vrt->wy = (int)crn[0][1];
    vrt->wz = (int)crn[0][2];
    vrt->nx = (float)sv[0];
    vrt->ny = (float)sv[1];
    vrt->nz = (float)sv[2];

    normal = QuantizeNormal(dv[0],dv[1],dv[2]);

    r = NewVertex();
    vrt = &Vrt[r];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[15][0];
    vrt->wy = (int)crn[15][1];
    vrt->wz = (int)crn[15][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    s = NewVertex();
    vrt = &Vrt[s];
    vrt->normal = normal;
    vrt->potential = 0.0;
    vrt->wx = (int)crn[8][0];
    vrt->wy = (int)crn[8][1];
    vrt->wz = (int)crn[8][2];
    vrt->nx = (float)dv[0];
    vrt->ny = (float)dv[1];
    vrt->nz = (float)dv[2];

    n[0] = sv[0] + dv[0];
    n[1] = sv[1] + dv[1];
    n[2] = sv[2] + dv[2];
    normal = QuantizeNormal(n[0],n[1],n[2]);

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = p;
    t->v = r;
    t->w = q;

    idx = NewTriangle();
    t = &Tri[idx];
    t->normal = normal;
    t->col = col;
    t->u = r;
    t->v = s;
    t->w = q;


    shade = Colour2Shade(col);
    Shade[shade].refcount += 16;
}


static void ElipseBegCap(TKnot *ptr, int col1, int col2)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register int normal, shade;
    register int idx, curr, prev;
    register unsigned int i;
    register Real tc, ts;
    register int cen;
    int cidx[EDGES];

    nx = -ptr->t[0];
    ny = -ptr->t[1];
    nz = -ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    cen = NewVertex();
    vrt = &Vrt[cen];
    vrt->potential = 0.0;
    vrt->normal = normal;
    vrt->wx = (int)ptr->p[0];
    vrt->wy = (int)ptr->p[1];
    vrt->wz = (int)ptr->p[2];
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    for( i=0; i<TEdges; i++ ) {
        tc = TC[i]*ptr->wide;
        ts = TS[i]*CartoonHeight;

        idx = NewVertex();
        vrt = &Vrt[idx];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)(ptr->p[0]+(tc*ptr->h[0] + ts*ptr->v[0]));
        vrt->wy = (int)(ptr->p[1]+(tc*ptr->h[1] + ts*ptr->v[1]));
        vrt->wz = (int)(ptr->p[2]+(tc*ptr->h[2] + ts*ptr->v[2]));
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;
        cidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
        idx = NewTriangle();
        t = &Tri[idx];
        t->col = ((i+i)<TEdges) ? col1 : col2;
        t->normal = normal;
        t->u = cen;
        t->v = cidx[curr];
        t->w = cidx[prev];
        prev = curr++;
    }
    shade = Colour2Shade(col1);
    Shade[shade].refcount += TEdges/2;
    shade = Colour2Shade(col2);
    Shade[shade].refcount += TEdges/2;
}


static void ElipseEndCap(TKnot *ptr, int col1, int col2)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register int normal, shade;
    register int idx, curr, prev;
    register unsigned int i;
    register Real tc, ts;
    register int cen;
    int cidx[EDGES];

    nx = ptr->t[0];
    ny = ptr->t[1];
    nz = ptr->t[2];
    normal = QuantizeNormal(nx,ny,nz);

    cen = NewVertex();
    vrt = &Vrt[cen];
    vrt->potential = 0.0;
    vrt->normal = normal;
    vrt->wx = (int)ptr->p[0];
    vrt->wy = (int)ptr->p[1];
    vrt->wz = (int)ptr->p[2];
    vrt->nx = (float)nx;
    vrt->ny = (float)ny;
    vrt->nz = (float)nz;

    for( i=0; i<TEdges; i++ ) {
        tc = TC[i]*ptr->wide;
        ts = TS[i]*CartoonHeight;

        idx = NewVertex();
        vrt = &Vrt[idx];
        vrt->potential = 0.0;
        vrt->normal = normal;
        vrt->wx = (int)(ptr->p[0]+(tc*ptr->h[0] + ts*ptr->v[0]));
        vrt->wy = (int)(ptr->p[1]+(tc*ptr->h[1] + ts*ptr->v[1]));
        vrt->wz = (int)(ptr->p[2]+(tc*ptr->h[2] + ts*ptr->v[2]));
        vrt->nx = (float)nx;
        vrt->ny = (float)ny;
        vrt->nz = (float)nz;
        cidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
        idx = NewTriangle();
        t = &Tri[idx];
        t->col = ((i+i)<TEdges) ? col1 : col2;
        t->normal = normal;
        t->u = cen;
        t->v = cidx[prev];
        t->w = cidx[curr];
        prev = curr++;
    }
    shade = Colour2Shade(col1);
    Shade[shade].refcount += TEdges/2;
    shade = Colour2Shade(col2);
    Shade[shade].refcount += TEdges/2;
}


static void CreateElipse(TKnot *src, TKnot *dst,
                         int col1, int col2)
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register int idx, curr, prev;
    register int normal, shade;
    register unsigned int i;
    register Real tc, ts;

    Real top[EDGES][3];
    Real bot[EDGES][3];
    int tidx[EDGES];
    int bidx[EDGES];
    Real n[3];

    for( i=0; i<TEdges; i++ ) {
      tc = TC[i]*CartoonHeight;
      ts = TS[i]*src->wide;
      n[0] = tc*src->h[0] + ts*src->v[0];
      n[1] = tc*src->h[1] + ts*src->v[1];
      n[2] = tc*src->h[2] + ts*src->v[2];
      normal = QuantizeNormal(n[0],n[1],n[2]);
      Normalize(n);

      tc = TC[i]*src->wide;
      ts = TS[i]*CartoonHeight;

      bot[i][0] = src->p[0]+(tc*src->h[0] + ts*src->v[0]);
      bot[i][1] = src->p[1]+(tc*src->h[1] + ts*src->v[1]);
      bot[i][2] = src->p[2]+(tc*src->h[2] + ts*src->v[2]);

      idx = NewVertex();
      vrt = &Vrt[idx];
      vrt->potential = 0.0;
      vrt->normal = normal;
      vrt->wx = (int)bot[i][0];
      vrt->wy = (int)bot[i][1];
      vrt->wz = (int)bot[i][2];
      vrt->nx = (float)n[0];
      vrt->ny = (float)n[1];
      vrt->nz = (float)n[2];
      bidx[i] = idx;

      tc = TC[i]*CartoonHeight;
      ts = TS[i]*dst->wide;
      n[0] = tc*dst->h[0] + ts*dst->v[0];
      n[1] = tc*dst->h[1] + ts*dst->v[1];
      n[2] = tc*dst->h[2] + ts*dst->v[2];
      normal = QuantizeNormal(n[0],n[1],n[2]);
      Normalize(n);

      tc = TC[i]*dst->wide;
      ts = TS[i]*CartoonHeight;

      top[i][0] = dst->p[0]+(tc*dst->h[0] + ts*dst->v[0]);
      top[i][1] = dst->p[1]+(tc*dst->h[1] + ts*dst->v[1]);
      top[i][2] = dst->p[2]+(tc*dst->h[2] + ts*dst->v[2]);

      idx = NewVertex();
      vrt = &Vrt[idx];
      vrt->potential = 0.0;
      vrt->normal = normal;
      vrt->wx = (int)top[i][0];
      vrt->wy = (int)top[i][1];
      vrt->wz = (int)top[i][2];
      vrt->nx = (float)n[0];
      vrt->ny = (float)n[1];
      vrt->nz = (float)n[2];
      tidx[i] = idx;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
      idx = NewTriangle();
      t = &Tri[idx];
      t->col = ((i+i)<TEdges) ? col1 : col2;
      t->u = bidx[curr];
      t->v = tidx[curr];
      t->w = tidx[prev];
      TriangleNormal1(t);

      idx = NewTriangle();
      t = &Tri[idx];
      t->col = ((i+i)<TEdges) ? col1 : col2;
      t->u = bidx[curr];
      t->v = tidx[prev];
      t->w = bidx[prev];
      TriangleNormal1(t);

      prev = curr++;
    }

    shade = Colour2Shade(col1);
    Shade[shade].refcount += TEdges;
    shade = Colour2Shade(col2);
    Shade[shade].refcount += TEdges;
}


#define STATE_INIT     0x00
#define STATE_NONE     0x01
#define STATE_DOTS     0x02
#define STATE_TRACE    0x03
#define STATE_CARTOON  0x04
#define STATE_ELIPSE   0x05

static void CreateRibbon( Chain __far *chain )
{
    register Group __far *group;
    register Atom __far *captr;
    register Atom __far *o1ptr;
    register Atom __far *o2ptr;
    register Atom __far *next;

    register int pcol1,pcol2;
    register int wide,arrow;
    register int col1,col2;
    register int pwide;
    register int prev;

    TKnot mid1, mid2, mid3, mid4, mid5, mid6, mid7;
    TKnot knot1, knot2;
    Real st[3], dt[3];
    Real b[3], d[3];

    pcol1 = 0;
    pcol2 = 0;
    pwide = 0;

    prev = STATE_INIT;
    group = chain->glist;
    if( IsProtein(group->refno) )
    {   captr = FindGroupAtom(group,1);
    } else captr = FindGroupAtom(group,7);

    while( group->gnext )
    {   
        if( IsProtein(group->gnext->refno) )
        {   next = FindGroupAtom(group->gnext,1);
            o1ptr = FindGroupAtom(group,3);
        } else /* Nucleic Acid */
        {   next = FindGroupAtom(group->gnext,7);
            o1ptr = FindGroupAtom(group->gnext,10);
        }

        if( !next || !captr || !o1ptr || (next->flag&BreakFlag) ||
            !((group->flag|group->gnext->flag)&DrawKnotFlag) )
        {
          switch (prev)
          {
          case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                               CreateSphere(knot1.p,pwide,pcol1);
                               break;
          case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                               break;
          case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                               break;
          }
          prev = STATE_INIT;
          group = group->gnext;
          captr = next;
          continue;
        }

        if( IsProtein(group->refno) )
        {   b[0] = o1ptr->xorg - captr->xorg;
            b[1] = o1ptr->yorg - captr->yorg;
            b[2] = o1ptr->zorg - captr->zorg;

        } else /* Nucleic Acid */
        {   o2ptr = FindGroupAtom(group,8);
            if( o2ptr && !FindGroupAtom(group,17) )
            {   /* Deoxyribonucleic Acid */
                b[0] = 0.5*(o1ptr->xorg + o2ptr->xorg) - captr->xorg;
                b[1] = 0.5*(o1ptr->yorg + o2ptr->yorg) - captr->yorg;
                b[2] = 0.5*(o1ptr->zorg + o2ptr->zorg) - captr->zorg;
            } else /* Ribonucleic Acid */
            {   b[0] = o1ptr->xorg - captr->xorg;
                b[1] = o1ptr->yorg - captr->yorg;
                b[2] = o1ptr->zorg - captr->zorg;
            }
        }
        Normalize(b);

        knot2.p[0] = 0.5*(captr->xorg + next->xorg);
        knot2.p[1] = 0.5*(captr->yorg + next->yorg);
        knot2.p[2] = 0.5*(captr->zorg + next->zorg);

        knot2.t[0] = next->xorg - captr->xorg;
        knot2.t[1] = next->yorg - captr->yorg;
        knot2.t[2] = next->zorg - captr->zorg;
        Normalize(knot2.t);

        /* c := a x b */
        knot2.v[0] = knot2.t[1]*b[2] - knot2.t[2]*b[1];
        knot2.v[1] = knot2.t[2]*b[0] - knot2.t[0]*b[2];
        knot2.v[2] = knot2.t[0]*b[1] - knot2.t[1]*b[0];
        Normalize(knot2.v);

        /* d := c x a */
        d[0] = knot2.v[1]*knot2.t[2] - knot2.v[2]*knot2.t[1];
        d[1] = knot2.v[2]*knot2.t[0] - knot2.v[0]*knot2.t[2];
        d[2] = knot2.v[0]*knot2.t[1] - knot2.v[1]*knot2.t[0];
        Normalize(d);

        /* Compensate for narrowing of helices! */
        if( (group->struc&group->gnext->struc) & HelixFlag )
        {   
            /* 1.00 Angstrom Displacement */
#ifdef INVERT
            knot2.p[0] += 250.0*knot2.v[0];
            knot2.p[1] += 250.0*knot2.v[1];
            knot2.p[2] += 250.0*knot2.v[2];
#else
            knot2.p[0] -= 250.0*knot2.v[0];
            knot2.p[1] -= 250.0*knot2.v[1];
            knot2.p[2] -= 250.0*knot2.v[2];
#endif
        }

        /* Handle Carbonyl Oxygen Flip */
        if( prev && ((knot1.h[0]*d[0] +
                      knot1.h[1]*d[1] +
                      knot1.h[2]*d[2])<0.0) )
        {   knot2.h[0] = -d[0];   knot2.v[0] = -knot2.v[0];
            knot2.h[1] = -d[1];   knot2.v[1] = -knot2.v[1];
            knot2.h[2] = -d[2];   knot2.v[2] = -knot2.v[2];

            if( !(col2 = group->col1) )
                col2 = captr->col;
            if( group->struc & HelixFlag )
            {   if( !(col1 = group->col2) )
                    col1 = captr->col;
            } else col1 = col2;
        } else
        {   knot2.h[0] = d[0];
            knot2.h[1] = d[1];
            knot2.h[2] = d[2];

            if( !(col1 = group->col1) )
                col1 = captr->col;
            if( group->struc & HelixFlag )
            {   if( !(col2 = group->col2) )
                    col2 = captr->col;
            } else col2 = col1;
        }

        arrow = False;
        if( group->flag&CartoonFlag )
        {   if( DrawBetaArrows && (group->struc&SheetFlag) &&
                !(group->gnext->struc&SheetFlag) )
            {   arrow = True;
                wide = 0;
            } else wide = group->width;
        } else if( group->flag & WideKnotFlag )
        {   /* Average Ribbon Width */
            if( group->gnext->flag & WideKnotFlag )
            {   wide = (group->width+group->gnext->width)>>1;
            } else if( group->gnext->flag & CartoonFlag )
            {   wide = group->gnext->width;
            } else wide = group->width;
        } else wide = group->gnext->width;
        knot2.wide = wide;

        if( prev )
        {
            /* Set the Hermite Spline Tension */
            if( (group->struc|group->gnext->struc) & SheetFlag ) {
                st[0] = 400.0 * knot1.t[0];
                st[1] = 400.0 * knot1.t[1];
                st[2] = 400.0 * knot1.t[2];

                dt[0] = 400.0 * knot2.t[0];
                dt[1] = 400.0 * knot2.t[1];
                dt[2] = 400.0 * knot2.t[2];
            } else {
                st[0] = 1400.0 * knot1.t[0];
                st[1] = 1400.0 * knot1.t[1];
                st[2] = 1400.0 * knot1.t[2];

                dt[0] = 1400.0 * knot2.t[0];
                dt[1] = 1400.0 * knot2.t[1];
                dt[2] = 1400.0 * knot2.t[2];
            }

            /* Calculate Hermite Spline Points */
            mid1.p[0] = 0.9570310*knot1.p[0] + 0.0957031*st[0]
                      + 0.0429688*knot2.p[0] - 0.0136719*dt[0];
            mid1.p[1] = 0.9570310*knot1.p[1] + 0.0957031*st[1]
                      + 0.0429688*knot2.p[1] - 0.0136719*dt[1];
            mid1.p[2] = 0.9570310*knot1.p[2] + 0.0957031*st[2]
                      + 0.0429688*knot2.p[2] - 0.0136719*dt[2];
            Interpolate(&knot1,&knot2,&mid1,0.9570310,0.0429688);

            mid2.p[0] = 0.84375*knot1.p[0] + 0.140625*st[0]
                      + 0.15625*knot2.p[0] - 0.046875*dt[0];
            mid2.p[1] = 0.84375*knot1.p[1] + 0.140625*st[1]
                      + 0.15625*knot2.p[1] - 0.046875*dt[1];
            mid2.p[2] = 0.84375*knot1.p[2] + 0.140625*st[2]
                      + 0.15625*knot2.p[2] - 0.046875*dt[2];
            Interpolate(&knot1,&knot2,&mid2,0.84375,0.15625);

            mid3.p[0] = 0.683594*knot1.p[0] + 0.1464840*st[0]
                      + 0.316406*knot2.p[0] - 0.0878906*dt[0];
            mid3.p[1] = 0.683594*knot1.p[1] + 0.1464840*st[1]
                      + 0.316406*knot2.p[1] - 0.0878906*dt[1];
            mid3.p[2] = 0.683594*knot1.p[2] + 0.1464840*st[2]
                      + 0.316406*knot2.p[2] - 0.0878906*dt[2];
            Interpolate(&knot1,&knot2,&mid3,0.683594,0.316406);

            mid4.p[0] = 0.5*knot1.p[0] + 0.125*st[0]
                      + 0.5*knot2.p[0] - 0.125*dt[0];
            mid4.p[1] = 0.5*knot1.p[1] + 0.125*st[1]
                      + 0.5*knot2.p[1] - 0.125*dt[1];
            mid4.p[2] = 0.5*knot1.p[2] + 0.125*st[2]
                      + 0.5*knot2.p[2] - 0.125*dt[2];
            Interpolate(&knot1,&knot2,&mid4,0.5,0.5);

            mid5.p[0] = 0.316406*knot1.p[0] + 0.0878906*st[0]
                      + 0.683594*knot2.p[0] - 0.1464840*dt[0];
            mid5.p[1] = 0.316406*knot1.p[1] + 0.0878906*st[1]
                      + 0.683594*knot2.p[1] - 0.1464840*dt[1];
            mid5.p[2] = 0.316406*knot1.p[2] + 0.0878906*st[2]
                      + 0.683594*knot2.p[2] - 0.1464840*dt[2];
            Interpolate(&knot1,&knot2,&mid5,0.316406,0.683594);

            mid6.p[0] = 0.15625*knot1.p[0] + 0.046875*st[0]
                      + 0.84375*knot2.p[0] - 0.140625*dt[0];
            mid6.p[1] = 0.15625*knot1.p[1] + 0.046875*st[1]
                      + 0.84375*knot2.p[1] - 0.140625*dt[1];
            mid6.p[2] = 0.15625*knot1.p[2] + 0.046875*st[2]
                      + 0.84375*knot2.p[2] - 0.140625*dt[2];
            Interpolate(&knot1,&knot2,&mid6,0.15625,0.84375);

            mid7.p[0] = 0.0429688*knot1.p[0] + 0.0136719*st[0]
                      + 0.9570310*knot2.p[0] - 0.0957031*dt[0];
            mid7.p[1] = 0.0429688*knot1.p[1] + 0.0136719*st[1]
                      + 0.9570310*knot2.p[1] - 0.0957031*dt[1];
            mid7.p[2] = 0.0429688*knot1.p[2] + 0.0136719*st[2]
                      + 0.9570310*knot2.p[2] - 0.0957031*dt[2];
            Interpolate(&knot1,&knot2,&mid7,0.0429688,0.9570310);


            if( group->flag & DotsFlag ) {
                switch (prev)
                {
                case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                     if (group->width >= pwide)
                                       CreateSphere(knot1.p,group->width,col1);
                                     else CreateSphere(knot1.p,pwide,pcol1);
                                     break;
                case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                                     CreateSphere(knot1.p,group->width,col1);
                                     break;
                case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                                     CreateSphere(knot1.p,group->width,col1);
                                     break;
                default:             CreateSphere(knot1.p,group->width,col1);
                }
                CreateSphere(mid4.p, group->width, col1);
                prev = STATE_DOTS;

            } else if( group->flag & TraceFlag ) {
                switch (prev)
                {
                case STATE_TRACE:    if (pwide != group->width) {
                                       TraceEndCap(&knot1,pwide,pcol1);
                                       if (group->width > pwide)
                                         CreateSphere(knot1.p,group->width,
                                                      col1);
                                       else CreateSphere(knot1.p,pwide,pcol1);
                                       TraceBegCap(&knot1,group->width,col1);
                                     }
                                     break;
                case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                                     CreateSphere(knot1.p,group->width,col1);
                                     TraceBegCap(&knot1,group->width,col1);
                                     break;
                case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                                     CreateSphere(knot1.p,group->width,col1);
                                     TraceBegCap(&knot1,group->width,col1);
                                     break;
                default:             CreateSphere(knot1.p,group->width,col1);
                                     TraceBegCap(&knot1,group->width,col1);
                }
                CreateTrace(&knot1, &mid1,  group->width, col1);
                CreateTrace(&mid1,  &mid2,  group->width, col1);
                CreateTrace(&mid2,  &mid3,  group->width, col1);
                CreateTrace(&mid3,  &mid4,  group->width, col1);
                CreateTrace(&mid4,  &mid5,  group->width, col1);
                CreateTrace(&mid5,  &mid6,  group->width, col1);
                CreateTrace(&mid6,  &mid7,  group->width, col1);
                CreateTrace(&mid7,  &knot2, group->width, col1);
                prev = STATE_TRACE;
                pwide = group->width;
                pcol1 = col1;

            } else if( group->flag & CartoonFlag ) {
                switch (prev)
                {
                case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                     CreateSphere(knot1.p,pwide,pcol1);
                                     CartoonBegCap(&knot1,col1);
                                     break;
                case STATE_CARTOON:  break;
                case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                                     CartoonBegCap(&knot1,col1);
                                     break;
                default:             CartoonBegCap(&knot1,col1);
                }
                if (arrow) {
                    ArrowBegCap(&knot1, col1);
                    CreateArrow(&knot1, &mid1,  col1);
                    CreateArrow(&mid1,  &mid2,  col1);
                    CreateArrow(&mid2,  &mid3,  col1);
                    CreateArrow(&mid3,  &mid4,  col1);
                    CreateArrow(&mid4,  &mid5,  col1);
                    CreateArrow(&mid5,  &mid6,  col1);
                    CreateArrow(&mid6,  &mid7,  col1);
                    CreateArrow(&mid7,  &knot2, col1);
                    prev = STATE_NONE;
                } else {
                    CreateCartoon(&knot1, &mid1,  col1);
                    CreateCartoon(&mid1,  &mid2,  col1);
                    CreateCartoon(&mid2,  &mid3,  col1);
                    CreateCartoon(&mid3,  &mid4,  col1);
                    CreateCartoon(&mid4,  &mid5,  col1);
                    CreateCartoon(&mid5,  &mid6,  col1);
                    CreateCartoon(&mid6,  &mid7,  col1);
                    CreateCartoon(&mid7,  &knot2, col1);
                    prev = STATE_CARTOON;
                    pcol1 = col1;
                }

            } else if( group->flag & DrawKnotFlag ) {
                switch (prev)
                {
                case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                     CreateSphere(knot1.p,pwide,pcol1);
                                     ElipseBegCap(&knot1,col1,col2);
                                     break;
                case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                                     ElipseBegCap(&knot1,col1,col2);
                                     break;
                case STATE_ELIPSE:   break;
                default:             ElipseBegCap(&knot1,col1,col2);
                }
                CreateElipse(&knot1, &mid1,  col1, col2);
                CreateElipse(&mid1,  &mid2,  col1, col2);
                CreateElipse(&mid2,  &mid3,  col1, col2);
                CreateElipse(&mid3,  &mid4,  col1, col2);
                CreateElipse(&mid4,  &mid5,  col1, col2);
                CreateElipse(&mid5,  &mid6,  col1, col2);
                CreateElipse(&mid6,  &mid7,  col1, col2);
                CreateElipse(&mid7,  &knot2, col1, col2);
                prev = STATE_ELIPSE;
                pcol1 = col1;
                pcol2 = col2;

            } else {
                switch (prev)
                {
                case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                     CreateSphere(knot1.p,pwide,pcol1);
                                     break;
                case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                                     break;
                case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                                     break;
                }
                prev = STATE_NONE;
            }
              
        } else if (group == chain->glist )
        {   knot1 = knot2;
            knot1.p[0] = captr->xorg;
            knot1.p[1] = captr->yorg;
            knot1.p[2] = captr->zorg;

            if( group->flag & DotsFlag ) {
                CreateSphere(knot1.p, group->width, col1);
                prev = STATE_DOTS;

            } else if( group->flag & TraceFlag ) {
                CreateSphere(knot1.p, group->width, col1);
                TraceBegCap(&knot1, group->width, col1);
                CreateTrace(&knot1, &knot2, group->width, col1);
                prev = STATE_TRACE;
                pwide = group->width;
                pcol1 = col1;

            } else if( group->flag & CartoonFlag ) {
                if (arrow) {
                     knot1.wide = group->width*1.5;
                     CartoonBegCap(&knot1, col1);
                     CreateCartoon(&knot1, &knot2, col1);
                     prev = STATE_NONE;
                } else {
                     CartoonBegCap(&knot1, col1);
                     CreateCartoon(&knot1, &knot2, col1);
                     prev = STATE_CARTOON;
                     pcol1 = col1;
                }

            } else if( group->flag & DrawKnotFlag ) {
                ElipseBegCap(&knot1, col1, col2);
                CreateElipse(&knot1, &knot2, col1, col2);
                prev = STATE_ELIPSE;
                pcol1 = col1;
                pcol2 = col2;

            } else {
                prev = STATE_NONE;
            }
        } else {
            prev = STATE_NONE;
        }
        group = group->gnext;
        captr = next;

        knot1 = knot2;
    }

    if( prev )
    {   if( !(col1 = group->col1) )
            col1 = captr->col;
        if( group->struc & HelixFlag )
        {   if( !(col2 = group->col2) )
                col2 = captr->col;
        } else col2 = col1;

        knot2.p[0] = captr->xorg;
        knot2.p[1] = captr->yorg;
        knot2.p[2] = captr->zorg;

        if( group->flag & DotsFlag ) {
            switch (prev)
            {
            case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                 if (group->width >= pwide ) {
                                   CreateSphere(knot1.p, group->width, col1);
                                 } else CreateSphere(knot1.p, pwide, pcol1);
                                 break;
            case STATE_CARTOON:  CartoonEndCap(&knot1, pcol1);
                                 CreateSphere(knot1.p, group->width, col1);
                                 break;
            case STATE_ELIPSE:   ElipseEndCap(&knot1, pcol1, pcol2);
                                 CreateSphere(knot1.p, group->width, col1);
                                 break;
            default:             CreateSphere(knot1.p, group->width, col1);
            }
            CreateSphere(knot2.p, group->width, col1);

        } else if( group->flag & TraceFlag ) {
            switch (prev)
            {
            case STATE_TRACE:   if( pwide != group->width ) {
                                  TraceEndCap(&knot1,pwide,pcol1);
                                  if (group->width >= pwide ) {
                                    CreateSphere(knot1.p, group->width, col1);
                                  } else CreateSphere(knot1.p, pwide, pcol1);
                                  TraceBegCap(&knot1,group->width,col1);
                                }
                                break;
            case STATE_CARTOON: CartoonEndCap(&knot1,pcol1);
                                CreateSphere(knot1.p,group->width,col1);
                                TraceBegCap(&knot1,group->width,col1);
                                break;
            case STATE_ELIPSE:  ElipseEndCap(&knot1,pcol1,pcol2);
                                CreateSphere(knot1.p,group->width,col1);
                                TraceBegCap(&knot1,group->width,col1);
                                break;
            default:            CreateSphere(knot1.p,group->width,col1);
                                TraceBegCap(&knot1,group->width,col1);
            }
            CreateTrace(&knot1, &knot2, group->width, col1);
            TraceEndCap(&knot2, group->width, col1);
            CreateSphere(knot2.p, group->width, col1);

        } else if( group->flag & CartoonFlag ) {
            switch (prev)
            {
            case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                 CreateSphere(knot1.p,pwide,pcol1);
                                 CartoonBegCap(&knot1,col1);
                                 break;
            case STATE_CARTOON:  break;
            case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                                 CartoonBegCap(&knot1,col1);
                                 break;
            default:             CartoonBegCap(&knot1,col1);
            }
            CreateCartoon(&knot1, &knot2, col1);
            CartoonEndCap(&knot2, col1);

        } else if( group->flag & DrawKnotFlag ) {
            switch (prev)
            {
            case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                 CreateSphere(knot1.p,pwide,pcol1);
                                 ElipseBegCap(&knot1, col1, col2);
                                 break;
            case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                                 ElipseBegCap(&knot1, col1, col2);
                                 break;
            case STATE_ELIPSE:   break;
            default:             ElipseBegCap(&knot1, col1, col2);
            }
            CreateElipse(&knot1, &knot2, col1, col2);
            ElipseEndCap(&knot2, col1, col2);

        } else {
            switch (prev)
            {
            case STATE_TRACE:    TraceEndCap(&knot1,pwide,pcol1);
                                 CreateSphere(knot1.p,pwide,pcol1);
                                 break;
            case STATE_CARTOON:  CartoonEndCap(&knot1,pcol1);
                                 break;
            case STATE_ELIPSE:   ElipseEndCap(&knot1,pcol1,pcol2);
                                 break;
            }
        }
    }
}


void CreateRibbons()
{
    register Chain __far *chain;
    register int res;

    if( ResolutionFlag != ResDefault ) {
        res = ResolutionFlag;
    } else res = ResMedium;
    InitialiseCAD(res);

    DeleteSurface();
    for( chain=Database->clist; chain; chain=chain->cnext )
        if( chain->glist )
            CreateRibbon( chain );
    ReDrawFlag |= RFRefresh;
}

