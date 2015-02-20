/* cad.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include <stdio.h>
#include <math.h>

#include "rasmol.h"
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "cmndline.h"
#include "scripts.h"
#include "render.h"
#include "repres.h"
#include "pixutils.h"
#include "tmesh.h"


#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif


/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)



/* FACES = 20*(4^ORDER)   */
/*       = 20<<(ORDER<<1) */

/* ORDER  FACES  VERTS */
/*   0      20     12  */
/*   1      80     42  */
/*   2     320    162  */
/*   3    1280    642  */

#define ORDER   3
#define FACES   1280
#define VERTS   642
#define EDGES   20

#define DOUBLE_DELTA  0.75
#define TRIPLE_DELTA  1.50

/* Face list of a unit icosahedron */
static int F[20][3] = {
    {  0,  2,  1 },  /*   0   19 */
    {  0,  1,  4 },  /*   1   18 */
    {  0,  3,  2 },  /*   2   17 */
    {  0,  8,  3 },  /*   3   16 */
    {  0,  4,  8 },  /*   4   15 */
    {  1,  2,  5 },  /*   5   14 */
    {  1,  6,  4 },  /*   6   13 */
    {  1,  5,  6 },  /*   7   12 */
    {  2,  3,  7 },  /*   8   11 */
    {  2,  7,  5 },  /*   9   10 */

    {  4, 11,  8 },  /*  10    9 */
    {  4,  6, 11 },  /*  11    8 */
    {  3,  8, 10 },  /*  12    7 */
    {  3, 10,  7 },  /*  13    6 */
    {  8, 11, 10 },  /*  14    5 */
    {  5,  7,  9 },  /*  15    4 */
    {  5,  9,  6 },  /*  16    3 */
    {  6,  9, 11 },  /*  17    2 */
    {  7, 10,  9 },  /*  18    1 */
    {  9, 10, 11 }   /*  19    0 */
            };


static double Modulus( Real *v )
{
    return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}


static void Normalize( Real *v )
{
    register double len;

    len = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
    v[0] /= len;  v[1] /= len;  v[2] /= len;
}


static int CreateVertex( Real *v )
{
    register unsigned int i;

    /* Reuse an existing vertex? */
    for( i=0; i<TVert; i++ )
        if( (fabs(TV[i][0]-v[0])<1.0e-5) &&
            (fabs(TV[i][1]-v[1])<1.0e-5) &&
            (fabs(TV[i][2]-v[2])<1.0e-5) )
            return i;

    TV[TVert][0] = v[0];
    TV[TVert][1] = v[1];
    TV[TVert][2] = v[2];
    return TVert++;
}


static void Tesselate( int pi, int qi, int ri, int d )
{
    auto Real u[3], v[3], w[3];
    register Real *p,*q,*r;
    register int ui,vi,wi;
    register Real *n;

    p = TV[pi];
    q = TV[qi];
    r = TV[ri];

    if( d-- ) {
        u[0] = p[0]+q[0];  v[0] = q[0]+r[0];  w[0] = r[0]+p[0];
        u[1] = p[1]+q[1];  v[1] = q[1]+r[1];  w[1] = r[1]+p[1];
        u[2] = p[2]+q[2];  v[2] = q[2]+r[2];  w[2] = r[2]+p[2];

        Normalize(u);  ui = CreateVertex(u);
        Normalize(v);  vi = CreateVertex(v);
        Normalize(w);  wi = CreateVertex(w);

        Tesselate( ui, vi, wi, d );
        Tesselate( pi, ui, wi, d );
        Tesselate( ui, qi, vi, d );
        Tesselate( wi, vi, ri, d );
    } else {
        TF[TFace][0] = pi;
        TF[TFace][1] = qi;
        TF[TFace][2] = ri;

        n = TN[TFace++];
        n[0] = (q[1]-p[1])*(r[2]-p[2]) - (r[1]-p[1])*(q[2]-p[2]);
        n[1] = (r[0]-p[0])*(q[2]-p[2]) - (q[0]-p[0])*(r[2]-p[2]);
        n[2] = (q[0]-p[0])*(r[1]-p[1]) - (r[0]-p[0])*(q[1]-p[1]);
        Normalize(n);
    }
}


void InitialiseCAD( int res )
{
    register unsigned int i;
    register double theta;

    TFace = 0;
    TVert = 12;

    for( i=0; i<20; i++ )
        Tesselate( F[i][0], F[i][1], F[i][2], res );

    TEdges = res*4+8;
    for( i=0; i<TEdges; i++ ) {
        theta = (2.0*PI*i)/TEdges;
        TC[i] = (Real)cos(theta);
        TS[i] = (Real)sin(theta);
    }
}


static void FatalCADError( char *ptr )
{
    InvalidateCmndLine();
    WriteString("CAD Error: Unable to create file `");
    WriteString( ptr );  WriteString("'!\n");
}


/*==========================*/
/*  Multiple Bond Handling  */
/*==========================*/

static int IsDoubleBond( Bond __far *ptr )
{
    if( !DrawDoubleBonds ) return False;
    return ptr->flag & DoubBondFlag;
}


static void MultipleBondVector( Bond __far *ptr, Real *v )
{
    register Bond __far *bptr;
    register Atom __far *src;
    register Atom __far *dst;
    auto Real s[3], d[3];
    auto Real r[3], x[3];
    register double len;
    register int init;

    src = ptr->srcatom;
    dst = ptr->dstatom;

    s[0] = src->xorg;
    s[1] = src->yorg;
    s[2] = src->zorg;

    d[0] = dst->xorg;
    d[1] = dst->yorg;
    d[2] = dst->zorg;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    x[0] = d[0] - s[0];
    x[1] = d[1] - s[1];
    x[2] = d[2] - s[2];
    Normalize(x);

    init = False;

    ForEachBond
        if( (bptr->srcatom==src) || (bptr->dstatom==src) ||
            (bptr->srcatom==dst) || (bptr->dstatom==dst) ) {
            if( bptr != ptr ) {
                if( bptr->srcatom == src ) {
                    r[0] = bptr->dstatom->xorg - src->xorg;
                    r[1] = bptr->dstatom->yorg - src->yorg;
                    r[2] = bptr->dstatom->zorg - src->zorg;
                } else if( bptr->dstatom == src ) {
                    r[0] = bptr->srcatom->xorg - src->xorg;
                    r[1] = bptr->srcatom->yorg - src->yorg;
                    r[2] = bptr->srcatom->zorg - src->zorg;
                } else if( bptr->srcatom == dst ) {
                    r[0] = bptr->dstatom->xorg - dst->xorg;
                    r[1] = bptr->dstatom->yorg - dst->yorg;
                    r[2] = bptr->dstatom->zorg - dst->zorg;
                } else /* bptr->dstatom == dst */ {
                    r[0] = bptr->srcatom->xorg - dst->xorg;
                    r[1] = bptr->srcatom->yorg - dst->yorg;
                    r[2] = bptr->srcatom->zorg - dst->zorg;
                }

#ifndef INVERT
                r[1] = -r[1];
#endif

                len = r[0]*x[0] + r[1]*x[1] + r[2]*x[2];

                r[0] -= len*x[0];
                r[1] -= len*x[1];
                r[2] -= len*x[2];

                Normalize(r);

                if( init ) {
                    len = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
                    if( len >= 0.0 ) {
                        v[0] += r[0];
                        v[1] += r[1];
                        v[2] += r[2];
                    } else {
                        v[0] -= r[0];
                        v[1] -= r[1];
                        v[2] -= r[2];
                    }
                } else {
                    v[0] = r[0];
                    v[1] = r[1];
                    v[2] = r[2];
                    init = True;
                }
            }
        }

    if( !init ) {
        /* Choose arbitrary plane */
        if( fabs(x[0]) < fabs(x[1]) ) {
            v[0] =  0.0;
            v[1] =  x[2];
            v[2] = -x[1];
        } else {
            v[0] =  x[2];
            v[1] =  0.0;
            v[2] = -x[0];
        }
    }
    Normalize(v);
}



/*===================*/
/*  Mark Atom Radii  */
/*===================*/

static unsigned int MarkHBonds( HBond __far *list, int mode )
{
    register unsigned int count;
    register HBond __far *ptr;
    register Atom __far *src;
    register Atom __far *dst;
    register short radius;

    count = 0;
    for( ptr=list; ptr; ptr=ptr->hnext )
        if( ptr->flag & CylinderFlag ) {
            radius = ptr->radius;
            if( radius > 0 ) {
                if( mode ) {
                    src = ptr->srcCA;
                    dst = ptr->dstCA;
                    if( !src || !dst )
                        continue;
                } else {
                    src = ptr->src;
                    dst = ptr->dst;
                }

                if( radius > src->mbox )
                    src->mbox = radius;
                if( radius > dst->mbox )
                    dst->mbox = radius;
                count++;
            }
        }
    return count;
}


static unsigned int MarkAtomRadii( void )
{
    register unsigned int count;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Bond __far *bptr;
    register Monitor *mptr;

    ForEachAtom
        if( aptr->flag & SphereFlag ) {
            aptr->mbox = aptr->radius;
        } else aptr->mbox = 0;

    count = 0;
    ForEachBond
        if( (bptr->flag&CylinderFlag) && !IsDoubleBond(bptr) ) {
            if( bptr->radius > bptr->srcatom->mbox )
                bptr->srcatom->mbox = bptr->radius;
            if( bptr->radius > bptr->dstatom->mbox )
                bptr->dstatom->mbox = bptr->radius;
            count++;
        }

    ForEachBack
        if( bptr->flag & CylinderFlag ) {
            if( bptr->radius > bptr->srcatom->mbox )
                bptr->srcatom->mbox = bptr->radius;
            if( bptr->radius > bptr->dstatom->mbox )
                bptr->dstatom->mbox = bptr->radius;
            count++;
        }

    count += MarkHBonds(Database->slist,SSBondMode);
    count += MarkHBonds(Database->hlist,HBondMode);

    if( MonitRadius ) {
        for( mptr=MonitList; mptr; mptr=mptr->next ) {
            if( mptr->src->mbox < MonitRadius )
                mptr->src->mbox = MonitRadius;
            if( mptr->dst->mbox < MonitRadius )
                mptr->dst->mbox = MonitRadius;
            count++;
        }
    }
    /* cylinder count! */
    return count;
}


/*===================*/
/*  Binary File I/O  */
/*===================*/

static void WriteBinaryInt( unsigned int x, FILE *fp )
{
    static unsigned int buffer;
    register unsigned char *ptr;
    register unsigned char ch;

    buffer = x;
    ptr = (unsigned char*)&buffer;
    if( !SwapBytes ) {
        ch = ptr[0];  ptr[0] = ptr[3];  ptr[3] = ch;
        ch = ptr[1];  ptr[1] = ptr[2];  ptr[2] = ch;
    }
    fwrite(ptr,4,1,fp);
}


static void WriteBinaryFloat( double x, FILE *fp )
{
    register unsigned char *ptr;
    register unsigned char ch;
    static float buffer;

    buffer = (float)x;
    ptr = (unsigned char*)&buffer;
    if( !SwapBytes ) {
        ch = ptr[0];  ptr[0] = ptr[3];  ptr[3] = ch;
        ch = ptr[1];  ptr[1] = ptr[2];  ptr[2] = ch;
    }
    fwrite(ptr,4,1,fp);
}


/*===============================*/
/*  STL Lithography File Format  */
/*===============================*/

static void WriteSTLTriangle( Real *p, Real *q, Real *r, Real *n, FILE *fp )
{
    if( BinaryFlag ) {
        WriteBinaryFloat(n[0],fp);
        WriteBinaryFloat(n[1],fp);
        WriteBinaryFloat(n[2],fp);

        WriteBinaryFloat(p[0],fp);
        WriteBinaryFloat(p[1],fp);
        WriteBinaryFloat(p[2],fp);

        WriteBinaryFloat(q[0],fp);
        WriteBinaryFloat(q[1],fp);
        WriteBinaryFloat(q[2],fp);

        WriteBinaryFloat(r[0],fp);
        WriteBinaryFloat(r[1],fp);
        WriteBinaryFloat(r[2],fp);

        fputc(0x00,fp);
        fputc(0x00,fp);
    } else /* ASCII */ {
        fprintf(fp,"  facet normal %g %g %g\n",n[0],n[1],n[2]);
        fputs("    outer loop\n",fp);

        fprintf(fp,"      vertex %g %g %g\n",p[0],p[1],p[2]);
        fprintf(fp,"      vertex %g %g %g\n",q[0],q[1],q[2]);
        fprintf(fp,"      vertex %g %g %g\n",r[0],r[1],r[2]);
        fputs("    endloop\n  endfacet\n",fp);
    }
}


static void WriteSTLTriangle2( Real *p, Real *q, Real *r, FILE *fp )
{
    auto Real n[3];

    n[0] = (q[1]-p[1])*(r[2]-p[2]) - (r[1]-p[1])*(q[2]-p[2]);
    n[1] = (r[0]-p[0])*(q[2]-p[2]) - (q[0]-p[0])*(r[2]-p[2]);
    n[2] = (q[0]-p[0])*(r[1]-p[1]) - (r[0]-p[0])*(q[1]-p[1]);
    Normalize(n);

    WriteSTLTriangle(p,q,r,n,fp);
}


static void WriteSTLSphere( double cx, double cy, double cz,
                            double rad, FILE *fp )
{
    static Real tmp[VERTS][3];
    register unsigned int i;
    register int *f;

    for( i=0; i<TVert; i++ ) {
        tmp[i][0] = TV[i][0]*rad + cx;
        tmp[i][1] = TV[i][1]*rad + cy;
        tmp[i][2] = TV[i][2]*rad + cz;
    }

    for( i=0; i<TFace; i++ ) {
         f = TF[i];
         WriteSTLTriangle(tmp[f[0]],tmp[f[1]],tmp[f[2]],TN[i],fp);
    }
}


static void WriteSTLCylinder( Real *src, Real *dst, double rad, FILE *fp )
{
    auto Real x[3], y[3];
    auto Real d[3], u[3];
    static Real top[EDGES][3];
    static Real bot[EDGES][3];
    register unsigned int curr,prev;
    register unsigned int i;
    register double rx,ry,rz;

    d[0] = dst[0] - src[0];
    d[1] = dst[1] - src[1];
    d[2] = dst[2] - src[2];
    Normalize(d);

    u[0] = -d[0];
    u[1] = -d[1];
    u[2] = -d[2];

    if( fabs(d[0]) < fabs(d[1]) ) {
        x[0] =  0.0;
        x[1] =  d[2];
        x[2] = -d[1];
    } else {
        x[0] =  d[2];
        x[1] =  0.0;
        x[2] = -d[0];
    }
    Normalize(x);

    y[0] = d[1]*x[2] - d[2]*x[1];
    y[1] = d[2]*x[0] - d[0]*x[2];
    y[2] = d[0]*x[1] - d[1]*x[0];
    Normalize(y);

    for( i=0; i<TEdges; i++ ) {
        rx = rad*(TC[i]*x[0]+TS[i]*y[0]);
        ry = rad*(TC[i]*x[1]+TS[i]*y[1]);
        rz = rad*(TC[i]*x[2]+TS[i]*y[2]);

        bot[i][0] = src[0]+rx;
        bot[i][1] = src[1]+ry;
        bot[i][2] = src[2]+rz;

        top[i][0] = dst[0]+rx;
        top[i][1] = dst[1]+ry;
        top[i][2] = dst[2]+rz;
    }

    curr = 0;
    prev = TEdges-1;
    for( i=0; i<TEdges; i++ ) {
         WriteSTLTriangle2(bot[curr],top[curr],top[prev],fp);
         WriteSTLTriangle2(bot[curr],top[prev],bot[prev],fp);
         WriteSTLTriangle(src,bot[curr],bot[prev],u,fp);
         WriteSTLTriangle(dst,top[prev],top[curr],d,fp);
         prev = curr++;
    }
}


static unsigned int CountSTLTriangles( void )
{
    register unsigned int count;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register double rad;

    count = MarkAtomRadii();

    /* Triangles per cylinder */
    count *= 4*TEdges;

    ForEachAtom
        if( aptr->mbox )
            count += TFace;

    if( DrawDoubleBonds ) {
        ForEachBond
            if( bptr->flag & CylinderFlag ) {
                if( bptr->flag & DoubBondFlag ) {
                    count += 8*TEdges;
                    rad = (DOUBLE_DELTA+1.0)*bptr->radius;
                    if( rad > (double)bptr->srcatom->mbox )
                        count += 2*TFace;
                    if( rad > (double)bptr->dstatom->mbox )
                        count += 2*TFace;
                } else if( bptr->flag & TripBondFlag ) {
                    count += 8*TEdges;
                    rad = (TRIPLE_DELTA+1.0)*bptr->radius;
                    if( rad > (double)bptr->srcatom->mbox )
                        count += 2*TFace;
                    if( rad > (double)bptr->dstatom->mbox )
                        count += 2*TFace;
                }
            }
    }

    /* Triangle Mesh */
    count += TriCount;
    return count;
}


static void WriteSTLAtoms( FILE *fp )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register double x,y,z,r;

    ForEachAtom
        if( aptr->mbox ) {
            x = (double)(aptr->xorg+OrigCX)/250.0;
            y = (double)(aptr->yorg+OrigCY)/250.0;
            z = (double)(aptr->zorg+OrigCZ)/250.0;
            r = (double)aptr->mbox/250.0;

            WriteSTLSphere(x,InvertY(y),z,r,fp);
        }
}


static void WriteSTLBond( Atom __far *src, Atom __far *dst,
                          int radius, FILE *fp )
{
    auto Real s[3];
    auto Real d[3];

    s[0] = (src->xorg+OrigCX)/250.0;
    s[1] = (src->yorg+OrigCY)/250.0;
    s[2] = (src->zorg+OrigCZ)/250.0;

    d[0] = (dst->xorg+OrigCX)/250.0;
    d[1] = (dst->yorg+OrigCY)/250.0;
    d[2] = (dst->zorg+OrigCZ)/250.0;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    WriteSTLCylinder(s,d,radius/250.0,fp);
}


static void WriteSTLMultipleBond( Bond __far *ptr, double delta, FILE *fp )
{
    register Atom __far *src;
    register Atom __far *dst;
    register double rad;
    auto Real s[3],os[3];
    auto Real d[3],od[3];
    auto Real v[3];

    src = ptr->srcatom;
    s[0] = (src->xorg+OrigCX)/250.0;
    s[1] = (src->yorg+OrigCY)/250.0;
    s[2] = (src->zorg+OrigCZ)/250.0;

    dst = ptr->dstatom;
    d[0] = (dst->xorg+OrigCX)/250.0;
    d[1] = (dst->yorg+OrigCY)/250.0;
    d[2] = (dst->zorg+OrigCZ)/250.0;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    MultipleBondVector(ptr,v);
    rad = ptr->radius/250.0;

    v[0] *= delta*rad;
    v[1] *= delta*rad;
    v[2] *= delta*rad;

    os[0] = s[0]+v[0];  od[0] = d[0]+v[0];
    os[1] = s[1]+v[1];  od[1] = d[1]+v[1];
    os[2] = s[2]+v[2];  od[2] = d[2]+v[2];

    if( (delta+1.0)*ptr->radius > (double)ptr->srcatom->mbox )
        WriteSTLSphere(os[0],os[1],os[2],rad,fp);
    WriteSTLCylinder(os,od,rad,fp);
    if( (delta+1.0)*ptr->radius > (double)ptr->dstatom->mbox )
        WriteSTLSphere(od[0],od[1],od[2],rad,fp);

    os[0] = s[0]-v[0];  od[0] = d[0]-v[0];
    os[1] = s[1]-v[1];  od[1] = d[1]-v[1];
    os[2] = s[2]-v[2];  od[2] = d[2]-v[2];

    if( (delta+1.0)*ptr->radius > (double)ptr->srcatom->mbox )
        WriteSTLSphere(os[0],os[1],os[2],rad,fp);
    WriteSTLCylinder(os,od,rad,fp);
    if( (delta+1.0)*ptr->radius > (double)ptr->dstatom->mbox )
        WriteSTLSphere(od[0],od[1],od[2],rad,fp);
}


static void WriteSTLHBonds( HBond __far *list, int mode, FILE *fp )
{
    register HBond __far *ptr;

    for( ptr=list; ptr; ptr=ptr->hnext )
        if( (ptr->flag&CylinderFlag) && (ptr->radius>0) ) {
            if( mode ) {
                if( ptr->srcCA && ptr->dstCA )
                    WriteSTLBond(ptr->srcCA,ptr->dstCA,ptr->radius,fp);
            } else WriteSTLBond(ptr->src,ptr->dst,ptr->radius,fp);
        }
}


static void WriteSTLBonds( FILE *fp )
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register Monitor *mptr;

    ForEachBond
        if( (bptr->flag&CylinderFlag) && !IsDoubleBond(bptr) )
            WriteSTLBond(bptr->srcatom,bptr->dstatom,bptr->radius,fp);

    if( DrawDoubleBonds ) {
        ForEachBond
            if( bptr->flag & CylinderFlag ) {
                if( bptr->flag & DoubBondFlag ) {
                    WriteSTLMultipleBond(bptr,DOUBLE_DELTA,fp);
                } else if( bptr->flag & TripBondFlag )
                    WriteSTLMultipleBond(bptr,TRIPLE_DELTA,fp);
            }
    }

    ForEachBack
        if( bptr->flag & CylinderFlag )
            WriteSTLBond(bptr->srcatom,bptr->dstatom,bptr->radius,fp);

    WriteSTLHBonds(Database->slist,SSBondMode,fp);
    WriteSTLHBonds(Database->hlist,HBondMode,fp);

    if( MonitRadius )
        for( mptr=MonitList; mptr; mptr=mptr->next )
            WriteSTLBond(mptr->src,mptr->dst,MonitRadius,fp);
}


static void WriteSTLTMesh( FILE *fp )
{
    Real pc[3],qc[3],rc[3];
    register VrtStruct *p, *q, *r;
    register TriStruct *t;
    register int i;

    for (i=0; i<TriCount; i++) {
        t = &Tri[i];
        p = &Vrt[t->u];
        q = &Vrt[t->v];
        r = &Vrt[t->w];

        pc[0] = (p->wx+OrigCX)/250.0;
        pc[1] = (p->wy+OrigCY)/250.0;
        pc[2] = (p->wz+OrigCZ)/250.0;

        qc[0] = (q->wx+OrigCX)/250.0;
        qc[1] = (q->wy+OrigCY)/250.0;
        qc[2] = (q->wz+OrigCZ)/250.0;

        rc[0] = (r->wx+OrigCX)/250.0;
        rc[1] = (r->wy+OrigCY)/250.0;
        rc[2] = (r->wz+OrigCZ)/250.0;

#ifndef INVERT
        pc[1] = -pc[1];
        qc[1] = -qc[1];
        rc[1] = -rc[1];
#endif

        WriteSTLTriangle2(pc,qc,rc,fp);
    }
}


int WriteSTLFile( char *name )
{
    register unsigned int count;
    register FILE *fp;
    register int res;
    register int i;

    if( !Database )
        return True;

    fp = fopen(name,BinaryFlag?"wb":"w");
    if( !fp ) {
        FatalCADError(name);
        return False;
    }

    if( ResolutionFlag != ResDefault ) {
        res = ResolutionFlag;
    } else res = ResMedium;

    InitialiseCAD(res);

    if( BinaryFlag ) {
        for( i=0; i<80; i++ )
            fputc(0x00,fp);
        count = CountSTLTriangles();
        WriteBinaryInt(count,fp);
    } else {
        fputs("solid\n",fp);
        MarkAtomRadii();
    }

    WriteSTLAtoms(fp);
    WriteSTLBonds(fp);
    WriteSTLTMesh(fp);

    if( !BinaryFlag )
        fputs("endsolid\n",fp);

    fclose(fp);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','TEXT',133);
#endif
    return True;
}


/*===========================*/
/*  AutoCAD DXF File Format  */
/*===========================*/

static void WriteDXFString( int code, char *ptr, FILE *fp )
{
    if( BinaryFlag ) {
        fputc(code,fp);
        while( *ptr ) {
            fputc(*ptr,fp);
            ptr++;
        }
        fputc(0x00,fp);
    } else fprintf(fp,"%d\n%s\n",code,ptr);
}


static void WriteDXFInteger( int code, int value, FILE *fp )
{
    if( BinaryFlag ) {
        fputc(code,fp);
        fputc(value&0xff,fp);
        fputc((value>>8)&0xff,fp);
    } else fprintf(fp,"%d\n%d\n",code,value);
}


static void WriteDXFReal( int code, double value, FILE *fp )
{
    register unsigned char *ptr;
    register unsigned char ch;


    if( BinaryFlag ) {
        fputc(code,fp);

        ptr = (unsigned char*)&value;

        if( !SwapBytes ) {
            ch = ptr[0];  ptr[0] = ptr[7];  ptr[7] = ch;
            ch = ptr[1];  ptr[1] = ptr[6];  ptr[6] = ch;
            ch = ptr[2];  ptr[2] = ptr[5];  ptr[5] = ch;
            ch = ptr[3];  ptr[3] = ptr[4];  ptr[4] = ch;
        }
        fwrite(ptr,1,8,fp);
    } else fprintf(fp,"%d\n%f\n",code,value);
}


static void WriteDXFEntities( FILE *fp )
{
    register VrtStruct *ptr;
    register Real x, y, z;
    register int i;

    for( i=0; i<TriCount; i++ ) {
        WriteDXFString(0,"3DFACE",fp);
        WriteDXFString(8,"0",fp);

        /* WriteDXFInteger(62,col,fp); */

        ptr = &Vrt[Tri[i].u];
        x = (ptr->wx+OrigCX)/250.0;
        y = (ptr->wy+OrigCY)/250.0;
        z = (ptr->wz+OrigCZ)/250.0;
#ifdef INVERT
        y = -y;
#endif
        WriteDXFReal(10,x,fp);
        WriteDXFReal(20,y,fp);
        WriteDXFReal(30,z,fp);

        ptr = &Vrt[Tri[i].v];
        x = (ptr->wx+OrigCX)/250.0;
        y = (ptr->wy+OrigCY)/250.0;
        z = (ptr->wz+OrigCZ)/250.0;
#ifdef INVERT
        y = -y;
#endif
        WriteDXFReal(11,x,fp);
        WriteDXFReal(21,y,fp);
        WriteDXFReal(31,z,fp);

        ptr = &Vrt[Tri[i].w];
        x = (ptr->wx+OrigCX)/250.0;
        y = (ptr->wy+OrigCY)/250.0;
        z = (ptr->wz+OrigCZ)/250.0;
#ifdef INVERT
        y = -y;
#endif
        WriteDXFReal(12,x,fp);
        WriteDXFReal(22,y,fp);
        WriteDXFReal(32,z,fp);

        WriteDXFReal(13,x,fp);
        WriteDXFReal(23,y,fp);
        WriteDXFReal(33,z,fp);
    }
}


int WriteDXFFile( char *name )
{
    register FILE *fp;

    if( !Database )
        return True;

    fp = fopen(name,BinaryFlag?"wb":"w");
    if( !fp ) {
        FatalCADError(name);
        return False;
    }

    if( BinaryFlag ) {
        fwrite("AutoCAD Binary DXF\015\012\032",1,22,fp);
    } /* else fputs("999\nCreated by RasMol v2.6\n",fp); */

    /* Header Section */
    WriteDXFString(0,"SECTION",fp);
    WriteDXFString(2,"HEADER",fp);
    WriteDXFString(0,"ENDSEC",fp);

    /* Table Section */
    WriteDXFString(0,"SECTION",fp);
    WriteDXFString(2,"TABLES",fp);

    WriteDXFString(0,"TABLE",fp);      /* LAYER TABLE WITH 1 ENTRY */
    WriteDXFString(2,"LAYER",fp);
    WriteDXFInteger(70,1,fp);

    WriteDXFString(0,"LAYER",fp);      /* DEFINE LAYER "0"         */
    WriteDXFString(2,"0",fp);
    WriteDXFInteger(70,0,fp);

    WriteDXFInteger(62,7,fp);          /* DEFAULT COLOR=7 (WHITE)  */
    WriteDXFString(6,"CONTINUOUS",fp); /* LINETYPE=CONTINUOUS      */

    WriteDXFString(0,"ENDTAB",fp);
    WriteDXFString(0,"ENDSEC",fp);

    /* Blocks Section */
    WriteDXFString(0,"SECTION",fp);
    WriteDXFString(2,"BLOCKS",fp);
    WriteDXFString(0,"ENDSEC",fp);

    /* Entities Section */
    WriteDXFString(0,"SECTION",fp);
    WriteDXFString(2,"ENTITIES",fp);

    WriteDXFEntities(fp);

    WriteDXFString(0,"ENDSEC",fp);
    WriteDXFString(0,"EOF",fp);

    fclose(fp);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','TEXT',133);
#endif
    return True;
}


/*===========================*/
/*  PLY Polygon File Format  */
/*===========================*/

static unsigned int CountPLYMultipleBondVertices( Bond __far *ptr,
                                                  double delta )
{
    register unsigned int count;
    register Atom __far *src;
    register Atom __far *dst;

    src = ptr->srcatom;
    dst = ptr->dstatom;

    count = 4*TEdges + 4;
    if( (delta+1.0)*ptr->radius > (double)src->mbox )
        count += 2*TVert;
    if( (delta+1.0)*ptr->radius > (double)dst->mbox )
        count += 2*TVert;
    if( !ptr->col && (src->col!=dst->col) )
        count += 2*TEdges;
    return count;
}


static unsigned int CountPLYHBondVertices( HBond __far *list, int mode )
{
    register unsigned int count;
    register HBond __far *ptr;

    count = 0;
    for( ptr=list; ptr; ptr=ptr->hnext )
        if( (ptr->flag&CylinderFlag) && (ptr->radius>0) ) {
            if( mode ) {
                if( ptr->srcCA && ptr->dstCA ) {
                    count += 2*TEdges + 2;
                    if( !ptr->col && (ptr->srcCA->col!=ptr->dstCA->col) )
                        count += TEdges;
                }
            } else {
                count += 2*TEdges + 2;
                if( !ptr->col && (ptr->src->col!=ptr->dst->col) )
                    count += TEdges;
            }
        }
    return count;
}


static unsigned int CountPLYVertices( void )
{
    register unsigned int count;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register Monitor *mptr;

    count = 0;
    ForEachAtom
        if( aptr->mbox )
            count += TVert;

    ForEachBond
        if( (bptr->flag&CylinderFlag) && !IsDoubleBond(bptr) ) {
            count += 2*TEdges+2;
            if( !bptr->col && (bptr->srcatom->col!=bptr->dstatom->col) )
                count += TEdges;
        }

    if( DrawDoubleBonds ) {
        ForEachBond
            if( bptr->flag & CylinderFlag ) {
                if( bptr->flag & DoubBondFlag ) {
                    count += CountPLYMultipleBondVertices(bptr,DOUBLE_DELTA);
                } else if( bptr->flag & TripBondFlag )
                    count += CountPLYMultipleBondVertices(bptr,TRIPLE_DELTA);
            }
    }

    ForEachBack
        if( bptr->flag & CylinderFlag ) {
            count += 2*TEdges+2;
            if( !bptr->col && (bptr->srcatom->col!=bptr->dstatom->col) )
                count += TEdges;
        }

    count += CountPLYHBondVertices(Database->slist,SSBondMode);
    count += CountPLYHBondVertices(Database->hlist,HBondMode);

    if( MonitRadius )
        for( mptr=MonitList; mptr; mptr=mptr->next ) {
            count += 2*TEdges+2;
            if( !mptr->col && (mptr->src->col!=mptr->dst->col) )
                count += TEdges;
        }

    /* Triangle Mesh */
    count += VrtCount;

    return count;
}


static unsigned int CountPLYMultipleBondFaces( Bond __far *ptr,
                                               double delta )
{
    register unsigned int count;
    register Atom __far *src;
    register Atom __far *dst;

    src = ptr->srcatom;
    dst = ptr->dstatom;

    count = 8*TEdges;
    if( (delta+1.0)*ptr->radius > (double)src->mbox )
        count += 2*TFace;
    if( (delta+1.0)*ptr->radius > (double)dst->mbox )
        count += 2*TFace;
    if( !ptr->col && (src->col!=dst->col) )
        count += 4*TEdges;
    return count;
}


static unsigned int CountPLYHBondFaces( HBond __far *list, int mode )
{
    register unsigned int count;
    register HBond __far *ptr;

    count = 0;
    for( ptr=list; ptr; ptr=ptr->hnext )
        if( (ptr->flag&CylinderFlag) && (ptr->radius>0) ) {
            if( mode ) {
                if( ptr->srcCA && ptr->dstCA ) {
                    count += 4*TEdges;
                    if( !ptr->col && (ptr->srcCA->col!=ptr->dstCA->col) )
                        count += 2*TEdges;
                }
            } else {
                count += 4*TEdges;
                if( !ptr->col && (ptr->src->col!=ptr->dst->col) )
                    count += 2*TEdges;
            }
        }
    return count;
}


static unsigned int CountPLYTriangles( void )
{
    register unsigned int count;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register Monitor *mptr;

    count = 0;
    ForEachAtom
        if( aptr->mbox )
            count += TFace;

    ForEachBond
        if( (bptr->flag&CylinderFlag) && !IsDoubleBond(bptr) ) {
            count += 4*TEdges;
            if( !bptr->col && (bptr->srcatom->col!=bptr->dstatom->col) )
                count += 2*TEdges;
        }

    if( DrawDoubleBonds ) {
        ForEachBond
            if( bptr->flag & CylinderFlag ) {
                if( bptr->flag & DoubBondFlag ) {
                    count += CountPLYMultipleBondFaces(bptr,DOUBLE_DELTA);
                } else if( bptr->flag & TripBondFlag )
                    count += CountPLYMultipleBondFaces(bptr,TRIPLE_DELTA);
            }
    }

    ForEachBack
        if( bptr->flag & CylinderFlag ) {
            count += 4*TEdges;
            if( !bptr->col && (bptr->srcatom->col!=bptr->dstatom->col) )
                count += 2*TEdges;
        }

    count += CountPLYHBondFaces(Database->slist,SSBondMode);
    count += CountPLYHBondFaces(Database->hlist,HBondMode);

    if( MonitRadius )
        for( mptr=MonitList; mptr; mptr=mptr->next ) {
            count += 4*TEdges;
            if( !mptr->col && (mptr->src->col!=mptr->dst->col) )
                count += 2*TEdges;
        }

    /* Triangle Mesh */
    count += TriCount;

    return count;
}


static void WritePLYVertex( double x, double y, double z, FILE *fp )
{
    if( BinaryFlag ) {
        WriteBinaryFloat(x,fp);
        WriteBinaryFloat(y,fp);
        WriteBinaryFloat(z,fp);
    } else fprintf(fp,"%g %g %g\n",x,y,z);
}


static void WritePLYTriangle( unsigned int u, unsigned int v, unsigned int w,
                              int col, FILE *fp )
{
    register ShadeDesc *shade;

    shade = Shade + Colour2Shade(col);
    if( BinaryFlag ) {
        fputc(shade->r,fp);
        fputc(shade->g,fp);
        fputc(shade->b,fp);

        fputc(0x03,fp);
        WriteBinaryInt(u,fp);
        WriteBinaryInt(v,fp);
        WriteBinaryInt(w,fp);
    } else
        fprintf(fp,"%d %d %d 3 %u %u %u\n",shade->r,shade->g,shade->b,u,v,w);
}


static void WritePLYSphere( double cx, double cy, double cz, double rad,
                            int col, unsigned int *index, FILE *fp )
{
    register unsigned int base;
    register unsigned int i;
    register double x,y,z;
    register int *f;

    if( index ) {
        base = *index;
        for( i=0; i<TFace; i++ ) {
            f = TF[i];
            WritePLYTriangle(f[0]+base,f[1]+base,f[2]+base,col,fp);
        }
        *index = base+TVert;
    } else {
        for( i=0; i<TVert; i++ ) {
            x = TV[i][0]*rad + cx;
            y = TV[i][1]*rad + cy;
            z = TV[i][2]*rad + cz;
            WritePLYVertex(x,y,z,fp);
        }
    }
}


static void WritePLYCylinder( Real *src, Real *dst, double rad,
                              int col1, int col2,
                              unsigned int *index, FILE *fp )
{
    auto Real d[3], m[3];
    auto Real x[3], y[3];
    static Real r[EDGES][3];
    register unsigned int curr,prev;
    register unsigned int base;
    register unsigned int top;
    register unsigned int mid;
    register unsigned int bot;
    register unsigned int i;

    if( index ) {
        curr = 0;
        prev = TEdges-1;

        base = *index;
        bot = base+2;
        top = bot+TEdges;
        if( col1 != col2 ) {
            mid = top+TEdges;
            *index = mid+TEdges;

            for( i=0; i<TEdges; i++ ) {
                WritePLYTriangle(bot+curr,mid+curr,mid+prev,col1,fp);
                WritePLYTriangle(bot+curr,mid+prev,bot+prev,col1,fp);
                WritePLYTriangle(mid+curr,top+curr,top+prev,col2,fp);
                WritePLYTriangle(mid+curr,top+prev,mid+prev,col2,fp);
                WritePLYTriangle(base+0,bot+curr,bot+prev,col1,fp);
                WritePLYTriangle(base+1,top+prev,top+curr,col2,fp);
                prev = curr++;
            }
        } else {
            for( i=0; i<TEdges; i++ ) {
                WritePLYTriangle(bot+curr,top+curr,top+prev,col1,fp);
                WritePLYTriangle(bot+curr,top+prev,bot+prev,col1,fp);
                WritePLYTriangle(base+0,bot+curr,bot+prev,col1,fp);
                WritePLYTriangle(base+1,top+prev,top+curr,col1,fp);
                prev = curr++;
            }
            *index = top+TEdges;
        }
    } else {
        d[0] = dst[0] - src[0];
        d[1] = dst[1] - src[1];
        d[2] = dst[2] - src[2];
        Normalize(d);

        if( fabs(d[0]) < fabs(d[1]) ) {
            x[0] =  0.0;
            x[1] =  d[2];
            x[2] = -d[1];
        } else {
            x[0] =  d[2];
            x[1] =  0.0;
            x[2] = -d[0];
        }
        Normalize(x);

        y[0] = d[1]*x[2] - d[2]*x[1];
        y[1] = d[2]*x[0] - d[0]*x[2];
        y[2] = d[0]*x[1] - d[1]*x[0];
        Normalize(y);

        for( i=0; i<TEdges; i++ ) {
            r[i][0] = rad*(TC[i]*x[0]+TS[i]*y[0]);
            r[i][1] = rad*(TC[i]*x[1]+TS[i]*y[1]);
            r[i][2] = rad*(TC[i]*x[2]+TS[i]*y[2]);
        }

        WritePLYVertex(src[0],src[1],src[2],fp);
        WritePLYVertex(dst[0],dst[1],dst[2],fp);

        for( i=0; i<TEdges; i++ )
            WritePLYVertex(src[0]+r[i][0],src[1]+r[i][1],src[2]+r[i][2],fp);
        for( i=0; i<TEdges; i++ )
            WritePLYVertex(dst[0]+r[i][0],dst[1]+r[i][1],dst[2]+r[i][2],fp);

        if( col1 != col2 ) {
            m[0] = 0.5*(src[0]+dst[0]);
            m[1] = 0.5*(src[1]+dst[1]);
            m[2] = 0.5*(src[2]+dst[2]);

            for( i=0; i<TEdges; i++ )
                WritePLYVertex(m[0]+r[i][0],m[1]+r[i][1],m[2]+r[i][2],fp);
        }
    }
}


static void WritePLYAtoms( FILE *fp, unsigned int *index )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register double x,y,z,r;

    ForEachAtom
        if( aptr->mbox ) {
            x = (double)(aptr->xorg+OrigCX)/250.0;
            y = (double)(aptr->yorg+OrigCY)/250.0;
            z = (double)(aptr->zorg+OrigCZ)/250.0;
            r = (double)aptr->mbox/250.0;

            WritePLYSphere(x,InvertY(y),z,r,aptr->col,index,fp);
        }
}


static void WritePLYBond( Atom __far *src, Atom __far *dst, int radius,
                          int col, unsigned int *index, FILE *fp )
{
    auto Real s[3];
    auto Real d[3];

    s[0] = (src->xorg+OrigCX)/250.0;
    s[1] = (src->yorg+OrigCY)/250.0;
    s[2] = (src->zorg+OrigCZ)/250.0;

    d[0] = (dst->xorg+OrigCX)/250.0;
    d[1] = (dst->yorg+OrigCY)/250.0;
    d[2] = (dst->zorg+OrigCZ)/250.0;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    if( !col ) {
        WritePLYCylinder(s,d,radius/250.0,src->col,dst->col,index,fp);
    } else WritePLYCylinder(s,d,radius/250.0,col,col,index,fp);
}


static void WritePLYMultipleBond( Bond __far *ptr, double delta,
                                  unsigned int *index, FILE *fp )
{
    register Atom __far *src;
    register Atom __far *dst;
    register int col1,col2;
    register double rad;
    auto Real s[3],os[3];
    auto Real d[3],od[3];
    auto Real v[3];

    src = ptr->srcatom;
    s[0] = (src->xorg+OrigCX)/250.0;
    s[1] = (src->yorg+OrigCY)/250.0;
    s[2] = (src->zorg+OrigCZ)/250.0;

    dst = ptr->dstatom;
    d[0] = (dst->xorg+OrigCX)/250.0;
    d[1] = (dst->yorg+OrigCY)/250.0;
    d[2] = (dst->zorg+OrigCZ)/250.0;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    MultipleBondVector(ptr,v);
    rad = ptr->radius/250.0;

    v[0] *= delta*rad;
    v[1] *= delta*rad;
    v[2] *= delta*rad;

    if( ptr->col ) {
        col1 = ptr->col;
        col2 = ptr->col;
    } else {
        col1 = src->col;
        col2 = dst->col;
    }

    os[0] = s[0]+v[0];  od[0] = d[0]+v[0];
    os[1] = s[1]+v[1];  od[1] = d[1]+v[1];
    os[2] = s[2]+v[2];  od[2] = d[2]+v[2];

    if( (delta+1.0)*ptr->radius > (double)src->mbox )
        WritePLYSphere(os[0],os[1],os[2],rad,col1,index,fp);
    WritePLYCylinder(os,od,rad,col1,col2,index,fp);
    if( (delta+1.0)*ptr->radius > (double)dst->mbox )
        WritePLYSphere(od[0],od[1],od[2],rad,col2,index,fp);

    os[0] = s[0]-v[0];  od[0] = d[0]-v[0];
    os[1] = s[1]-v[1];  od[1] = d[1]-v[1];
    os[2] = s[2]-v[2];  od[2] = d[2]-v[2];

    if( (delta+1.0)*ptr->radius > (double)src->mbox )
        WritePLYSphere(os[0],os[1],os[2],rad,col1,index,fp);
    WritePLYCylinder(os,od,rad,col1,col2,index,fp);
    if( (delta+1.0)*ptr->radius > (double)dst->mbox )
        WritePLYSphere(od[0],od[1],od[2],rad,col2,index,fp);
}


static void WritePLYHBonds( HBond __far *list, int mode,
                            unsigned int *index, FILE *fp )
{
    register HBond __far *ptr;

    for( ptr=list; ptr; ptr=ptr->hnext )
        if( (ptr->flag&CylinderFlag) && (ptr->radius>0) ) {
            if( mode ) {
                if( ptr->srcCA && ptr->dstCA )
                    WritePLYBond(ptr->srcCA,ptr->dstCA,ptr->radius,
                                 ptr->col,index,fp);
            } else WritePLYBond(ptr->src,ptr->dst,ptr->radius,
                                ptr->col,index,fp);
        }
}

 
static void WritePLYBonds( FILE *fp, unsigned int *index )
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register Monitor *mptr;

    ForEachBond
        if( (bptr->flag&CylinderFlag) && !IsDoubleBond(bptr) )
            WritePLYBond(bptr->srcatom,bptr->dstatom,bptr->radius,
                         bptr->col,index,fp);

    if( DrawDoubleBonds ) {
        ForEachBond
            if( bptr->flag & CylinderFlag ) {
                if( bptr->flag & DoubBondFlag ) {
                    WritePLYMultipleBond(bptr,DOUBLE_DELTA,index,fp);
                } else if( bptr->flag & TripBondFlag )
                    WritePLYMultipleBond(bptr,TRIPLE_DELTA,index,fp);
            }
    }

    ForEachBack
        if( bptr->flag & CylinderFlag )
            WritePLYBond(bptr->srcatom,bptr->dstatom,bptr->radius,
                         bptr->col,index,fp);

    WritePLYHBonds(Database->slist,SSBondMode,index,fp);
    WritePLYHBonds(Database->hlist,HBondMode,index,fp);

    if( MonitRadius )
        for( mptr=MonitList; mptr; mptr=mptr->next )
            WritePLYBond(mptr->src,mptr->dst,MonitRadius,
                         mptr->col,index,fp);
}


static void WritePLYTMesh( FILE *fp, unsigned int *index )
{
    register unsigned int base;
    register VrtStruct *v;
    register TriStruct *t;
    register Real x, y, z;
    register int i;

    if (index) {
        base = *index;
        for( i=0; i<TriCount; i++ ) {
            t = &Tri[i];
            WritePLYTriangle(t->u+base,t->v+base,t->w+base,t->col,fp);
        }
        base += VrtCount;
    } else {
        for( i=0; i<VrtCount; i++ ) {
            v = &Vrt[i];
            x = (v->wx+OrigCX)/250.0;
            y = (v->wy+OrigCY)/250.0;
            z = (v->wz+OrigCZ)/250.0;
#ifndef INVERT
            y = -y;
#endif
            WritePLYVertex(x, y, z, fp);
        }
    }
}


int WritePLYFile( char *name )
{
    unsigned int tricount;
    unsigned int vrtcount;
    register FILE *fp;
    register int res;

    if( !Database )
        return True;

    fp = fopen(name,BinaryFlag?"wb":"w");
    if( !fp ) {
        FatalCADError(name);
        return False;
    }

    if( ResolutionFlag != ResDefault ) {
        res = ResolutionFlag;
    } else res = ResMedium;

    InitialiseCAD(res);
    MarkAtomRadii();

    vrtcount = CountPLYVertices();
    tricount = CountPLYTriangles();

    fputs("ply\n",fp);
    if( BinaryFlag ) {
        fputs("format binary_little_endian 1.0\n",fp);
    } else fputs("format ascii 1.0\n",fp);
    fputs("comment created by RasMol v2.6\n",fp);
    fprintf(fp,"element vertex %u\n",vrtcount);
    fputs("property float x\n",fp);
    fputs("property float y\n",fp);
    fputs("property float z\n",fp);
    fprintf(fp,"element face %u\n",tricount);
    fputs("property uchar red\n",fp);
    fputs("property uchar green\n",fp);
    fputs("property uchar blue\n",fp);
    fputs("property list uchar int vertex_index\n",fp);
    fputs("end_header\n",fp);

    WritePLYAtoms(fp,(unsigned int*)0);
    WritePLYBonds(fp,(unsigned int*)0);
    WritePLYTMesh(fp,(unsigned int*)0);

    vrtcount = 0;
    WritePLYAtoms(fp,&vrtcount);
    WritePLYBonds(fp,&vrtcount);
    WritePLYTMesh(fp,&vrtcount);

    fclose(fp);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','TEXT',133);
#endif
    return True;
}


/*===========================================*/
/*  VRML Virtual Reality Modelling Language  */
/*===========================================*/


static int VRMLShade;

static void WriteVRMLTriple( double x, double y, double z, FILE *fp )
{
#ifdef __STDC__
    fprintf(fp,"%.3f %.3f %.3f",x,y,z);
#else
    fprintf(fp,"%.3lf %.3lf %.3lf",x,y,z);
#endif
}

static void WriteVRMLColour( int indent, int shade, FILE *fp )
{
    register int i;

    for( i=0; i<indent; i++ )
        putc(' ',fp);
    fputs("Material { diffuseColor ",fp);
    WriteVRMLTriple(Shade[shade].r/255.0,
                    Shade[shade].g/255.0,
                    Shade[shade].b/255.0,fp);
    fputs(" }\n",fp);
}


static void WriteVRMLSphere( Real *p, double rad, int col, FILE *fp )
{
    register int shade;

    shade = Colour2Shade(col);
    if( shade != VRMLShade ) {
        WriteVRMLColour(2,shade,fp);
        VRMLShade = shade;
    }

    fputs("  Separator {\n",fp);
    fputs("    Translation { translation ",fp);
    fprintf(fp,"%.3f %.3f %.3f }\n",p[0],p[1],p[2]);
    fprintf(fp,"    Sphere { radius %.3f }\n",rad);
    fputs("  }\n",fp);
}


static void WriteVRMLAtoms( FILE *fp )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register double ox,oy,oz;
    register double x,y,z;
    register int i,flag;

    for( i=0; i<LastShade; i++ )
        if( Shade[i].refcount ) {
            flag = False;
            ox = 0.0;  oy = 0.0;  oz = 0.0;
            ForEachAtom
                if( aptr->mbox && (Colour2Shade(aptr->col)==i) ) {
                    if( !flag ) {
                        WriteVRMLColour(2,i,fp);
                        fputs("  Separator {\n",fp);
                        VRMLShade = i;
                        flag = True;
                    }

                    x = (aptr->xorg+OrigCX)/250.0;
                    y = (aptr->yorg+OrigCY)/250.0;
                    z = (aptr->zorg+OrigCZ)/250.0;

                    fputs("    Translation { translation ",fp);
                    WriteVRMLTriple(x-ox,InvertY(y-oy),z-oz,fp);
                    ox = x;  oy = y;  oz = z;
                    fputs(" }\n",fp);

                    fputs("    Sphere { radius ",fp);
#ifdef __STDC__
                    fprintf(fp,"%.3f }\n",(double)aptr->mbox/250.0);
#else
                    fprintf(fp,"%.3f }\n",(double)aptr->mbox/250.0);
#endif
                }

            if( flag )
                fputs("  }\n",fp);
        }
}


static void WriteVRMLCylinder( Real *src, Real *dst, double rad,
                               int col, FILE *fp )
{
    auto Real mid[3];
    auto Real x[3],y[3],z[3];
    register double len;
    register int shade;

    mid[0] = 0.5*(src[0]+dst[0]);
    mid[1] = 0.5*(src[1]+dst[1]);
    mid[2] = 0.5*(src[2]+dst[2]);

    y[0] = dst[0]-src[0];
    y[1] = dst[1]-src[1];
    y[2] = dst[2]-src[2];

    len = Modulus(y);
    y[0] /= len;
    y[1] /= len;
    y[2] /= len;

    if( fabs(y[0]) < fabs(y[1]) ) {
        x[0] =  0.0;
        x[1] =  y[2];
        x[2] = -y[1];
    } else {
        x[0] =  y[2];
        x[1] =  0.0;
        x[2] = -y[0];
    }
    Normalize(x);

    z[0] = x[1]*y[2] - x[2]*y[1];
    z[1] = x[2]*y[0] - x[0]*y[2];
    z[2] = x[0]*y[1] - x[1]*y[0];
    Normalize(z);

    shade = Colour2Shade(col);
    if( shade != VRMLShade ) {
        WriteVRMLColour(2,shade,fp);
        VRMLShade = shade;
    }

    fputs("  Separator {\n",fp);
    fputs("    MatrixTransform { matrix\n",fp);
    fprintf(fp,"      %8.5f %8.5f %8.5f 0\n",x[0],x[1],x[2]);
    fprintf(fp,"      %8.5f %8.5f %8.5f 0\n",y[0],y[1],y[2]);
    fprintf(fp,"      %8.5f %8.5f %8.5f 0\n",z[0],z[1],z[2]);
    fprintf(fp,"      %8.3f %8.3f %8.3f 1\n",mid[0],mid[1],mid[2]);
    fputs("    }\n",fp);
    fputs("    Cylinder { ",fp);
    /* fputs("parts SIDES ",fp); */
    fprintf(fp,"radius %.3f height %.3f }\n",rad,len);
    fputs("  }\n",fp);
}


static void WriteVRMLBond( Atom __far *src, Atom __far *dst,
                           int radius, int col, FILE *fp )
{
    register double rad;
    auto Real s[3];
    auto Real d[3];
    auto Real m[3];

    s[0] = (src->xorg+OrigCX)/250.0;
    s[1] = (src->yorg+OrigCY)/250.0;
    s[2] = (src->zorg+OrigCZ)/250.0;

    d[0] = (dst->xorg+OrigCX)/250.0;
    d[1] = (dst->yorg+OrigCY)/250.0;
    d[2] = (dst->zorg+OrigCZ)/250.0;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    rad = radius/250.0;

    if( col ) {
        WriteVRMLCylinder(s,d,rad,col,fp);
    } else if( src->col == dst->col ) {
        WriteVRMLCylinder(s,d,rad,src->col,fp);
    } else {
        m[0] = 0.5*(s[0]+d[0]);
        m[1] = 0.5*(s[1]+d[1]);
        m[2] = 0.5*(s[2]+d[2]);
        WriteVRMLCylinder(s,m,rad,src->col,fp);
        WriteVRMLCylinder(m,d,rad,dst->col,fp);
    }
}


static void WriteVRMLMultipleBond( Bond __far *ptr, double delta, FILE *fp )
{ 
    register Atom __far *src;
    register Atom __far *dst;
    register double rad;
    auto Real s[3],os[3];
    auto Real d[3],od[3];
    auto Real m[3],v[3];
    register int col;

    src = ptr->srcatom;
    s[0] = (src->xorg+OrigCX)/250.0;
    s[1] = (src->yorg+OrigCY)/250.0;
    s[2] = (src->zorg+OrigCZ)/250.0;

    dst = ptr->dstatom;
    d[0] = (dst->xorg+OrigCX)/250.0;
    d[1] = (dst->yorg+OrigCY)/250.0;
    d[2] = (dst->zorg+OrigCZ)/250.0;

#ifndef INVERT
    s[1] = -s[1];
    d[1] = -d[1];
#endif

    MultipleBondVector(ptr,v);
    rad = ptr->radius/250.0;

    v[0] *= delta*rad;
    v[1] *= delta*rad;
    v[2] *= delta*rad;

    os[0] = s[0]+v[0];  od[0] = d[0]+v[0];
    os[1] = s[1]+v[1];  od[1] = d[1]+v[1];
    os[2] = s[2]+v[2];  od[2] = d[2]+v[2];

    if( ptr->col ) {
        col = ptr->col;
    } else if( src->col == dst->col ) {
        col = src->col;
    } else col = 0;

    if( col ) {
        if( (delta+1.0)*ptr->radius > (double)ptr->srcatom->mbox )
            WriteVRMLSphere(os,rad,col,fp);
        WriteVRMLCylinder(os,od,rad,col,fp);
        if( (delta+1.0)*ptr->radius > (double)ptr->dstatom->mbox )
            WriteVRMLSphere(od,rad,col,fp);
    } else {
        m[0] = 0.5*(os[0]+od[0]);
        m[1] = 0.5*(os[1]+od[1]);
        m[2] = 0.5*(os[2]+od[2]);
        if( (delta+1.0)*ptr->radius > (double)ptr->srcatom->mbox )
            WriteVRMLSphere(os,rad,src->col,fp);
        WriteVRMLCylinder(os,m,rad,src->col,fp);
        WriteVRMLCylinder(m,od,rad,dst->col,fp);
        if( (delta+1.0)*ptr->radius > (double)ptr->dstatom->mbox )
            WriteVRMLSphere(od,rad,dst->col,fp);
    }

    os[0] = s[0]-v[0];  od[0] = d[0]-v[0];
    os[1] = s[1]-v[1];  od[1] = d[1]-v[1];
    os[2] = s[2]-v[2];  od[2] = d[2]-v[2];

    if( col ) {
        if( (delta+1.0)*ptr->radius > (double)ptr->srcatom->mbox )
            WriteVRMLSphere(os,rad,col,fp);
        WriteVRMLCylinder(os,od,rad,col,fp);
        if( (delta+1.0)*ptr->radius > (double)ptr->dstatom->mbox )
            WriteVRMLSphere(od,rad,col,fp);
    } else {
        m[0] = 0.5*(os[0]+od[0]);
        m[1] = 0.5*(os[1]+od[1]);
        m[2] = 0.5*(os[2]+od[2]);
        /* These are reversed to reduce colour changes! */
        if( (delta+1.0)*ptr->radius > (double)ptr->dstatom->mbox )
            WriteVRMLSphere(od,rad,dst->col,fp);
        WriteVRMLCylinder(od,m,rad,dst->col,fp);
        WriteVRMLCylinder(m,os,rad,src->col,fp);
        if( (delta+1.0)*ptr->radius > (double)ptr->srcatom->mbox )
            WriteVRMLSphere(os,rad,src->col,fp);
    }
}


static void WriteVRMLHBonds( HBond __far *list, int mode, FILE *fp )
{
    register HBond __far *ptr;

    for( ptr=list; ptr; ptr=ptr->hnext )
        if( (ptr->flag&CylinderFlag) && (ptr->radius>0) ) {
            if( mode ) {
                if( ptr->srcCA && ptr->dstCA )
                    WriteVRMLBond(ptr->srcCA,ptr->dstCA,ptr->radius,
                                  ptr->col,fp);
            } else WriteVRMLBond(ptr->src,ptr->dst,ptr->radius,ptr->col,fp);
        }
}


static void WriteVRMLBonds( FILE *fp )
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register Monitor *mptr;

    ForEachBond
        if( (bptr->flag&CylinderFlag) && !IsDoubleBond(bptr) )
            WriteVRMLBond(bptr->srcatom,bptr->dstatom,
                          bptr->radius,bptr->col,fp);

    if( DrawDoubleBonds ) {
        ForEachBond
            if( bptr->flag & CylinderFlag ) {
                if( bptr->flag & DoubBondFlag ) {
                    WriteVRMLMultipleBond(bptr,DOUBLE_DELTA,fp);
                } else if( bptr->flag & TripBondFlag )
                    WriteVRMLMultipleBond(bptr,TRIPLE_DELTA,fp);
            }
    }


    ForEachBack
        if( bptr->flag & CylinderFlag )
            WriteVRMLBond(bptr->srcatom,bptr->dstatom,
                          bptr->radius,bptr->col,fp);

    WriteVRMLHBonds(Database->slist,SSBondMode,fp);
    WriteVRMLHBonds(Database->hlist,HBondMode,fp);

    if( MonitRadius )
        for( mptr=MonitList; mptr; mptr=mptr->next )
            WriteVRMLBond(mptr->src,mptr->dst,MonitRadius,mptr->col,fp);
}


static void WriteVRMLLine( int src, int dst, int shade, int *flag, FILE *fp )
{
    if( !*flag ) {
        WriteVRMLColour(4,shade,fp);
        fputs("    IndexedLineSet {\n",fp);
        fputs("      coordIndex [\n",fp);
        *flag = True;
    }
    fprintf(fp,"        %5d, %5d, -1,\n",src-1,dst-1);
}


static void WriteVRMLWireframe( FILE *fp )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Bond __far *bptr;
    register Atom __far *src;
    register Atom __far *dst;
    register double x,y,z;
    register int i,j;
    static int flag;

    ForEachAtom
        aptr->mbox = 0;

    i = 0;
    ForEachBond
        if( bptr->flag & WireFlag ) {
            src = bptr->srcatom;
            dst = bptr->dstatom;
            if( i == 0 ) {
                fputs("  Separator {\n",fp);
                fputs("    Coordinate3 {\n",fp);
                fputs("      point [\n",fp);
            }

            if( !src->mbox ) {
                x = (src->xorg+OrigCX)/250.0;
                y = (src->yorg+OrigCY)/250.0;
                z = (src->zorg+OrigCZ)/250.0;

                fputs("        ",fp);
                WriteVRMLTriple(x,InvertY(y),z,fp);
                fputs(",\n",fp);
                src->mbox = ++i;
            }
            
            if( !dst->mbox ) {
                x = (dst->xorg+OrigCX)/250.0;
                y = (dst->yorg+OrigCY)/250.0;
                z = (dst->zorg+OrigCZ)/250.0;

                fputs("        ",fp);
                WriteVRMLTriple(x,InvertY(y),z,fp);
                fputs(",\n",fp);
                dst->mbox = ++i;
            }

            if( !bptr->col && (src->col!=dst->col) ) {
                x = (0.5*(src->xorg+dst->xorg)+OrigCX)/250.0;
                y = (0.5*(src->yorg+dst->yorg)+OrigCY)/250.0;
                z = (0.5*(src->zorg+dst->zorg)+OrigCZ)/250.0;
                
                fputs("        ",fp);
                WriteVRMLTriple(x,InvertY(y),z,fp);
                fputs(",\n",fp);
                i++;
            }
        }

    /* No wireframe! */
    if( !i )  return;

    fputs("      ]\n",fp);
    fputs("    }\n",fp);
    
    for( j=0; j<LastShade; j++ )
        if( Shade[j].refcount ) {
            i = 1;
            flag = False;
            ForEachBond
                if( bptr->flag & WireFlag ) {
                    src = bptr->srcatom;   
                    dst = bptr->dstatom;

                    if( src->mbox == i ) i++;
                    if( dst->mbox == i ) i++;

                    if( bptr->col ) {
                        if( Colour2Shade(bptr->col) == j )
                            WriteVRMLLine(src->mbox,dst->mbox,j,&flag,fp);
                    } else if( src->col == dst->col ) {
                        if( Colour2Shade(src->col) == j )
                            WriteVRMLLine(src->mbox,dst->mbox,j,&flag,fp);
                    } else /* Two Colour Bond */ {
                        if( Colour2Shade(src->col) == j ) {
                            WriteVRMLLine(src->mbox,i,j,&flag,fp);
                        } else if( Colour2Shade(dst->col) == j )
                            WriteVRMLLine(dst->mbox,i,j,&flag,fp);
                        i++;
                    }
                }

            if( flag ) {
                fputs("      ]\n",fp);
                fputs("    }\n",fp);
            }
        }
    fputs("  }\n",fp);
}


static void WriteVRMLDots( FILE *fp )
{
    auto int hist[LastShade];
    register DotStruct __far *ptr;
    register double x,y,z;
    register int count;
    register int flag;
    register int i,j;

    flag = False;
    for( i=0; i<LastShade; i++ )
        if( Shade[i].refcount ) {
            count = 0;

            for( ptr=DotPtr; ptr; ptr=ptr->next )
                for( j=0; j<ptr->count; j++ )
                    if( Colour2Shade(ptr->col[j]) == i ) {
                        if( !flag ) {
                            fputs("  Separator {\n",fp);
                            fputs("    Coordinate3 {\n",fp);
                            fputs("      point [\n",fp);
                            flag = True;
                        }

                        x = (double)ptr->xpos[j]/250.0;
                        y = (double)ptr->ypos[j]/250.0;
                        z = (double)ptr->zpos[j]/250.0;

                        fputs("        ",fp);
                        WriteVRMLTriple(x,InvertY(y),z,fp);
                        fputs(",\n",fp);
                        count++;
                    }

            hist[i] = count;
        } else hist[i] = 0;

    if( flag ) {
        fputs("      ]\n",fp);
        fputs("    }\n",fp);

        count = 0;
        for( i=0; i<LastShade; i++ )
            if( hist[i] ) {
                WriteVRMLColour(4,i,fp);
                fputs("    PointSet {\n      ",fp);
                fprintf(fp,"startIndex %d numPoints %d\n",count,hist[i]);
                fputs("    }\n",fp);
                count += hist[i];
            }

        fputs("  }\n",fp);
    }
}


static void WriteVRMLTMesh( FILE *fp )
{
    register VrtStruct *v;
    register Real nx, ny, nz;
    register Real x, y, z;
    register int init;
    register int i,j;

    /* Coordinates */
    fputs("  Coordinate3 {\n    point [\n",fp);
    for( i=0; i<VrtCount; i++ ) {
        v = &Vrt[i];
        x = (v->wx+OrigCX)/250.0;
        y = (v->wy+OrigCY)/250.0;
        z = (v->wz+OrigCZ)/250.0;
#ifndef INVERT
        y = -y;
#endif
        fprintf(fp,"      %g %g %g,\n",x,y,z);
    }
    fputs("    ]\n  }\n",fp);

    /* Normals */
    fputs("  Normal {\n    vector [\n",fp);
    for( i=0; i<VrtCount; i++ ) {
        v = &Vrt[i];
        nx = Normal[v->normal][0];
        ny = Normal[v->normal][1];
        nz = Normal[v->normal][2];
#ifndef INVERT
        ny = -ny;
#endif
        fprintf(fp,"      %g %g %g,\n",nx,ny,nz);
    }
    fputs("    ]\n  }\n",fp);

    /* Triangles */
    for( j=0; j<LastShade; j++ ) {
        init = False;
        for( i=0; i<TriCount; i++ ) {
            if( Colour2Shade(Tri[i].col) == j ) {
                if( !init ) {
                    fputs("  Material { diffuseColor ",fp);
                    fprintf(fp,"%.3f %.3f %.3f",Shade[j].r/255.0,
                                                Shade[j].g/255.0,
                                                Shade[j].b/255.0);
                    fputs(" }\n  IndexedFaceSet {\n    coordIndex [\n",fp);
                    init = True;
                }
                fprintf(fp,"      %d, %d, %d, -1,\n",Tri[i].u,
                                                     Tri[i].v,
                                                     Tri[i].w);
            }
        }
        if( init )
            fputs("    ]\n  }\n",fp);
    }
}


int WriteVRMLFile( char *name )
{
    register FILE *fp;

    if( !Database )
        return True;

    fp = fopen(name,"w");
    if( !fp ) {
        FatalCADError(name);
        return False;
    }

    fputs("#VRML V1.0 ascii\n",fp);
    fputs("#Created by RasMol v2.6\n\n",fp);

    fputs("DEF Viewer Info { string \"examiner\" }\n",fp);
    fprintf(fp,"DEF Title Info { string \"%s\" }\n",Info.moleculename);
    fputs("DEF Creator Info { string \"Created by RasMol v2.6\" }\n",fp);
    fputs("DEF BackgroundColor Info { string \"",fp);
    WriteVRMLTriple(BackR/255.0,BackG/255.0,BackB/255.0,fp);
    fputs("\" }\n\n",fp);
    
    fputs("Separator {\n",fp);

#ifdef ORIG
    fputs("  DirectionalLight {\n",fp);
    fputs("    direction -1 -1 -2\n",fp);
    fputs("  }\n\n",fp);
#endif

    VRMLShade = -1;

    MarkAtomRadii();
    WriteVRMLAtoms(fp);
    WriteVRMLBonds(fp);
    WriteVRMLWireframe(fp);
    if( DotPtr )
        WriteVRMLDots(fp);
    if( TriCount )
        WriteVRMLTMesh(fp);

    fputs("}\n",fp);
    fclose(fp);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','TEXT',133);
#endif
    return True;
}


/*=============================*/
/*  Grasp Surface File Format  */
/*=============================*/

static void WriteGraspInt( FILE *fp, int value )
{
    register char *ptr;
    register char ch;

    ptr = (char*)&value;
    if( SwapBytes )
    {   ch = ptr[0];  ptr[0] = ptr[3];  ptr[3] = ch;
        ch = ptr[1];  ptr[1] = ptr[2];  ptr[2] = ch;
    }
    fwrite(ptr,sizeof(int),1,fp);
}

static void WriteGraspFloat( FILE *fp, float value )
{
    register char *ptr;
    register char ch;

    ptr = (char*)&value;
    if( SwapBytes )
    {   ch = ptr[0];  ptr[0] = ptr[3];  ptr[3] = ch;
        ch = ptr[1];  ptr[1] = ptr[2];  ptr[2] = ch;
    }
    fwrite(ptr,sizeof(float),1,fp);
}

static void WriteGraspLine( FILE *fp, const char *buffer, int len )
{
    register const char *ptr = buffer;
    register int i;

    WriteGraspInt(fp,len);
    for( i=0; i<len; i++ )
        if( *ptr ) {
            putc(*ptr,fp);
            ptr++;
        } else putc(' ',fp);
    WriteGraspInt(fp,len);
}

static void WriteGraspHeader( FILE *fp, int vcount, int tcount )
{
    char buffer[80];

    WriteGraspLine(fp,"format=2",80);
    WriteGraspLine(fp,"vertices,normals,triangles",80);
    WriteGraspLine(fp,"",80);

    sprintf(buffer,"%d %d 65 1.0",vcount,tcount);
    WriteGraspLine(fp,buffer,80);

    WriteGraspLine(fp,"0.0 0.0 0.0",80);
}


int WriteGraspFile( char *name )
{
    register TriStruct *t;
    register VrtStruct *vrt;
    register Real nx, ny, nz;
    register Real x, y, z;
    register int u, v, w;
    register FILE *fp;
    register int i;

    if( TriCount == 0 )
        return True;

    fp = fopen(name,"wb");
    if( !fp ) {
        FatalCADError(name);
        return False;
    }

    WriteGraspHeader(fp,VrtCount,TriCount);

    /* Write Vertices */
    WriteGraspInt(fp,12*VrtCount);
    for( i=0; i<VrtCount; i++ ) {
        vrt = &Vrt[i];
        x = (vrt->wx+OrigCX)/250.0;
        y = (vrt->wy+OrigCY)/250.0;
        z = (vrt->wz+OrigCZ)/250.0;
#ifdef INVERT
        y = -y;
#endif
        WriteGraspFloat(fp,x);
        WriteGraspFloat(fp,y);
        WriteGraspFloat(fp,-z);
    }
    WriteGraspInt(fp,12*VrtCount);

    /* Write Normals */
    WriteGraspInt(fp,12*VrtCount);
    for( i=0; i<VrtCount; i++ ) {
        vrt = &Vrt[i];
        nx = vrt->nx;
        ny = vrt->ny;
        nz = vrt->nz;
#ifdef INVERT
        ny = -ny;
#endif
        WriteGraspFloat(fp,nx);
        WriteGraspFloat(fp,ny);
        WriteGraspFloat(fp,-nz);
    }
    WriteGraspInt(fp,12*VrtCount);

    /* Write Triangles */
    WriteGraspInt(fp,12*TriCount);
    for( i=0; i<TriCount; i++ ) {
        t = &Tri[i];
#ifdef INVERT
        u = t->u;
        v = t->v;
        w = t->w;
#else
        u = t->u;
        v = t->w;
        w = t->v;
#endif
        WriteGraspInt(fp,u+1);
        WriteGraspInt(fp,v+1);
        WriteGraspInt(fp,w+1);
    }
    WriteGraspInt(fp,12*TriCount);

    fclose(fp);
#ifdef APPLEMAC
    SetFileInfo(name,'ttxt','TEXT',133);
#endif
    return True;
}
