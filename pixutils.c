/* pixutils.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#include <Types.h>
#include <Errors.h>
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#endif
#ifdef sun386
#include <stdlib.h>
#endif

#include <stdio.h>
#include <math.h>

#define PIXUTILS
#include "pixutils.h"
#include "graphics.h"
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"
#include "font.h"

#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif

/* Sutherland-Cohen Line Clipping Macros */
#define BitAbove    0x01
#define BitBelow    0x02
#define BitRight    0x04
#define BitLeft     0x08
#define BitFront    0x10

#define Reject(x,y)   ((x)&(y))
#define Accept(x,y)   (!((x)|(y)))
#define RootSix       2.44948974278

/* These define light source position */
#define LightDot(x,y,z)  ((x)+InvertY(y)+(z)+(z))
#define LightLength      RootSix


typedef struct {
        Long dx,dz,di;
        Long x,z,i;
    } Edge;

typedef struct {
        short dx,dy,dz;
        short inten;
        Long offset;
    } ArcEntry;


/* Note: DrawCylinderCaps currently employs an
 *       extremely crude hack to avoid stripes
 *       appearing along cylinders.
 */
#define ARCSIZE  8192

static ArcEntry __far *ArcAcPtr;
static ArcEntry __far *ArcDnPtr;
#if defined(IBMPC) || defined(APPLEMAC)
static ArcEntry __far *ArcAc;
static ArcEntry __far *ArcDn;
#else
static ArcEntry ArcAc[ARCSIZE];
static ArcEntry ArcDn[ARCSIZE];
#endif

static char FontDimen[23];
static int ClipStatus;



#define SETPIXEL(dptr,fptr,d,c)    if( (d) > *(dptr) )              \
                                   {   *(dptr) = (d);               \
                                       *(fptr) = (c);               \
                                   }

#define OutCode2(res,x,y)             \
    {   if( (y)<0 )                   \
        {   (res) = BitAbove;         \
        } else if( (y) >= View.ymax ) \
        {   (res) = BitBelow;         \
        } else (res) = 0;             \
                                      \
        if( (x) < 0 )                 \
        {   (res) |= BitLeft;         \
        } else if( (x) >= View.xmax ) \
            (res) |= BitRight;        \
    }

#define OutCode3(res,x,y,z)           \
    {   OutCode2(res,x,y);            \
        if( !ZValid((z)) )            \
            (res) |= BitFront;        \
    }



/*=======================*/
/*  Function Prototypes  */
/*=======================*/

static void DrawArcDn( short __huge*, Pixel __huge*, int, int );
static void DrawArcAc( short __huge*, Pixel __huge*, int, int );
static void ClipArcDn( short __huge*, Pixel __huge*, int, int, int, int );
static void ClipArcAc( short __huge*, Pixel __huge*, int, int, int, int );


void PlotPoint( int x, int y, int z, int col )
{
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register Long offset;

    /* SETPIXEL(dptr,fptr,z,Lut[col]); */

    offset = (Long)y*View.yskip+x;
    dptr = View.dbuf+offset;
    if( z > *dptr )
    {   fptr = View.fbuf+offset;
        *fptr = Lut[col];
        *dptr = z;
    }
}


void ClipPoint( int x, int y, int z, int col )
{
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register Long offset;

    if( XValid(x) && YValid(y) && ZValid(z) )
    {   /* PlotPoint(x,y,z,col); */
        offset = (Long)y*View.yskip+x;
        dptr = View.dbuf+offset;
        if( z > *dptr )
        {   fptr = View.fbuf+offset;
            *fptr = Lut[col];
            *dptr = z;
        }
    }
}


void PlotDeepPoint( int x, int y, int z, int col )
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int inten;

    offset = (Long)y*View.yskip+x;
    dptr = View.dbuf+offset;

    if( z > *dptr )
    {   fptr = View.fbuf+offset;
        inten = (ColourDepth*(z+ImageRadius-ZOffset))/ImageSize;
        if( inten > 0 )
        {   *fptr = Lut[col+(inten&ColourMask)];
        } else *fptr = Lut[col];
        *dptr = z;
    }
}


void ClipDeepPoint( int x, int y, int z, int col )
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int inten;

    if( XValid(x) && YValid(y) && ZValid(z) )
    {   /* PlotDeepPoint(x,y,z,col); */
        offset = (Long)y*View.yskip+x;
        dptr = View.dbuf+offset;

        if( z > *dptr )
        {  fptr = View.fbuf+offset;
           inten = (ColourDepth*(z+ImageRadius-ZOffset))/ImageSize;
           *fptr = Lut[col+inten];
           *dptr = z;
        }
    }
}



/*================================================*/
/*  Macros for Bresenhams Line Drawing Algorithm  */
/*================================================*/

#define CommonStep(s)  z1 += zrate; SETPIXEL(dptr,fptr,z1,c);     \
                       if( (zerr+=dz)>0 ) { zerr-=(s); z1+=iz; }

#define XStep  { if((err+=dy)>0) { fptr+=ystep; dptr+=ystep; err-=dx; } \
                 fptr+=ix; dptr+=ix; x1+=ix; CommonStep(dx); }

#define YStep  { if((err+=dx)>0) { fptr+=ix; dptr+=ix; err-=dy; } \
                 fptr+=ystep; dptr+=ystep; y1+=iy; CommonStep(dy); }
                     

void DrawTwinLine( int x1, int y1, int z1,
                   int x2, int y2, int z2,
                   int col1, int col2 )
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int zrate, zerr;
    register int ystep,err;
    register int ix,iy,iz;
    register int dx,dy,dz;
    register int mid;
    register Pixel c;

    c = Lut[col1];

    offset = (Long)y1*View.yskip + x1;
    fptr = View.fbuf+offset;
    dptr = View.dbuf+offset;

    SETPIXEL(dptr,fptr,z1,c);

    dx = x2-x1;  dy = y2-y1; 
    if( !dx && !dy ) return;
    dz = z2-z1;

    if( dy<0 ) 
    {   ystep = -View.yskip;
        dy = -dy; 
        iy = -1;
    } else
    {   ystep = View.yskip;
        iy = 1;
    }

    if( dx<0 ) 
    {   dx = -dx;
        ix = -1;
    } else ix = 1;

    if( dz<0 ) 
    {   dz = -dz;
        iz = -1;
    } else iz = 1;

    if( dx>dy )
    {   if( dz >= dx )
        {   zrate = dz/dx;
            dz -= dx*zrate;
            if( iz < 0 )
                zrate = -zrate;
        } else zrate = 0;
        err = zerr = -(dx>>1);

        if( col1 != col2 )
        {   mid = (x1+x2)>>1;
            while( x1!=mid ) XStep;
            c = Lut[col2];
        }
        while( x1!=x2 ) XStep;

    } else
    {   if( dz >= dy )
        {   zrate = dz/dy;
            dz -= dy*zrate;
            if( iz < 0 )
                zrate = -zrate;
        } else zrate = 0;
        err = zerr = -(dy>>1);

        if( col1 != col2 )
        {   mid = (y1+y2)>>1;
            while( y1!=mid ) YStep;
            c = Lut[col2];
        }
        while( y1!=y2 ) YStep;
    }
}


void ClipLine( int x1, int y1, int z1,
               int x2, int y2, int z2,
               int col )
{
    register int code1,code2;
    register int delta,rest;
    register int temp;

    OutCode3(code1,x1,y1,z1);
    OutCode3(code2,x2,y2,z2);
    if( Reject(code1,code2) )
        return;

    while( !Accept(code1,code2) )
    {   if( !code1 )
        {   temp=x1; x1=x2; x2=temp;
            temp=y1; y1=y2; y2=temp;
            temp=z1; z1=z2; z2=temp;
            code1 = code2;
            code2 = 0;
        }

        if( code1 & BitAbove )
        {   delta = y2-y1;
            x1 += (int)(((Long)y1*(x1-x2))/delta);  
            z1 += (int)(((Long)y1*(z1-z2))/delta);
            y1 = 0;
        } else if( code1 & BitLeft )
        {   delta = x2-x1;
            y1 += (int)(((Long)x1*(y1-y2))/delta);
            z1 += (int)(((Long)x1*(z1-z2))/delta);
            x1 = 0;
        } else if( code1 & BitRight )
        {   delta = x2-x1;
            temp=View.xmax-1; rest=temp-x1;
            y1 += (int)(((Long)rest*(y2-y1))/delta);
            z1 += (int)(((Long)rest*(z2-z1))/delta);
            x1 = temp;
        } else if( code1 & BitBelow )
        {   delta = y2-y1;
            temp=View.ymax-1; rest=temp-y1;
            x1 += (int)(((Long)rest*(x2-x1))/delta);
            z1 += (int)(((Long)rest*(z2-z1))/delta);
            y1 = temp;
        } else /* SLAB */
        {   delta = z2-z1;
            rest = (SlabValue-1)-z1;
            x1 += (int)(((Long)rest*(x2-x1))/delta);
            y1 += (int)(((Long)rest*(y2-y1))/delta);
            z1 = SlabValue-1;
        }

        OutCode3(code1,x1,y1,z1);
        if( Reject(code1,code2) )
            return;
    }
    DrawTwinLine(x1,y1,z1,x2,y2,z2,col,col);
}


void ClipTwinLine( int x1, int y1, int z1,
                   int x2, int y2, int z2,
                   int col1, int col2 )
{
    register int xmid,ymid,zmid;
    register int code1,code2;

    if( col1!=col2 )
    {   OutCode3(code1,x1,y1,z1);
        OutCode3(code2,x2,y2,z2);
        if( !Reject(code1,code2) )
        {   if( !Accept(code1,code2) )
            {  xmid = (x1+x2)/2;
               ymid = (y1+y2)/2;
               zmid = (z1+z2)/2;
               ClipLine(x1,y1,z1,xmid,ymid,zmid,col1);
               ClipLine(xmid,ymid,zmid,x2,y2,z2,col2);
            } else
               DrawTwinLine(x1,y1,z1,x2,y2,z2,col1,col2);
        }
    } else ClipLine(x1,y1,z1,x2,y2,z2,col1);
}



/*=============================================*/
/*  Macros for 3D Bresenhams Vector Algorithm  */
/*=============================================*/

#define CommonVectStep(s)  z1 += zrate;   c1 += crate;                    \
                           SETPIXEL(dptr,fptr,z1,Lut[col+c1]);            \
                           if( (zerr+=dz)>0 ) { zerr -= (s); z1 += iz; }  \
                           if( (cerr+=dc)>0 ) { cerr -= (s); c1 += iz; }

#define XVectStep  { if((err+=dy)>0) { fptr+=ystep; dptr+=ystep; err-=dx; } \
                     fptr+=ix; dptr+=ix; x1+=ix; CommonVectStep(dx); }

#define YVectStep  { if((err+=dx)>0) { fptr+=ix; dptr+=ix; err-=dy; } \
                     fptr+=ystep; dptr+=ystep; y1+=iy; CommonVectStep(dy); }


void DrawTwinVector( int x1, int y1, int z1,
                     int x2, int y2, int z2,
                     int col1, int col2 )
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int dx,dy,dz,dc;
    register int crate, cerr;
    register int zrate, zerr;
    register int ystep,err;
    register int ix,iy,iz;
    register int col, mid;
    register int c1, c2;

    c1 = (ColourDepth*(z1+ImageRadius-ZOffset))/ImageSize;
    c2 = (ColourDepth*(z2+ImageRadius-ZOffset))/ImageSize;

    offset = (Long)y1*View.yskip + x1;
    fptr = View.fbuf+offset;
    dptr = View.dbuf+offset;

    SETPIXEL(dptr,fptr,z1,Lut[col1+c1]);

    dx = x2 - x1;  dy = y2 - y1;
    dz = z2 - z1;  dc = c2 - c1;
    if( !dx && !dy ) return;

    if( dy<0 ) 
    {   ystep = -View.yskip;
        dy = -dy; 
        iy = -1; 
    } else
    {   ystep = View.yskip;
        iy = 1;
    }

    if( dx<0 ) 
    {   dx = -dx; 
        ix = -1; 
    } else ix = 1;

    iz = (dz<0)? -1 : 1;

    if( dx>dy )
    {   if( dz >= dx )
        {   zrate = dz/dx;
            dz -= dx*zrate;
            if( iz < 0 )
                zrate = -zrate;
        } else zrate = 0;

        if( dc >= dx )
        {   crate = dc/dx;
            dc -= dx*crate;
            if( iz < 0 )
                crate = -crate;
        } else crate = 0;

        err = zerr = cerr = -(dx>>1);
        col = col1;

        if( dz<0 )
        {   dz = -dz;
            dc = -dc;
        }
        
        if( col1 != col2 )
        {   mid = (x1+x2)>>1;
            while( x1!=mid ) XVectStep;
            col = col2;
        }
        while( x1!=x2 ) XVectStep;
    } else
    {   if( dz >= dy )
        {   zrate = dz/dy;
            dz -= dy*zrate;
            if( iz < 0 )
                zrate = -zrate;
        } else zrate = 0;

        if( dc >= dy )
        {   crate = dc/dy;
            dc -= dy*crate;
            if( iz < 0 )
                crate = -crate;
        } else crate = 0;

        err = zerr = cerr = -(dy>>1);
        col = col1;

        if( dz<0 )
        {   dz = -dz;
            dc = -dc;
        }

        if( col1 != col2 )
        {   mid = (y1+y2)>>1;
            while( y1!=mid ) YVectStep;
            col=col2;
        }
        while( y1!=y2 ) YVectStep;
    }
}


static void ClipVector( int x1, int y1, int z1,
                        int x2, int y2, int z2,
                        int col )
{
    register int code1,code2;
    register int delta,rest;
    register int temp;

    OutCode3(code1,x1,y1,z1);
    OutCode3(code2,x2,y2,z2);
    if( Reject(code1,code2) )
        return;

    while( !Accept(code1,code2) )
    {   if( !code1 )
        {   temp=x1; x1=x2; x2=temp;
            temp=y1; y1=y2; y2=temp;
            temp=z1; z1=z2; z2=temp;
            code1 = code2;
            code2 = 0;
        }

        if( code1 & BitAbove )
        {   delta = y2-y1;
            x1 += (int)(((Long)y1*(x1-x2))/delta);  
            z1 += (int)(((Long)y1*(z1-z2))/delta);
            y1 = 0;
        } else if( code1 & BitLeft )
        {   delta = x2-x1;
            y1 += (int)(((Long)x1*(y1-y2))/delta);
            z1 += (int)(((Long)x1*(z1-z2))/delta);
            x1 = 0;
        } else if( code1 & BitRight )
        {   delta = x2-x1;
            temp=View.xmax-1; rest=temp-x1;
            y1 += (int)(((Long)rest*(y2-y1))/delta);
            z1 += (int)(((Long)rest*(z2-z1))/delta);
            x1 = temp;
        } else if( code1 & BitBelow )
        {   delta = y2-y1;
            temp=View.ymax-1; rest=temp-y1;
            x1 += (int)(((Long)rest*(x2-x1))/delta);
            z1 += (int)(((Long)rest*(z2-z1))/delta);
            y1 = temp;
        } else /* SLAB */
        {   delta = z2-z1;
            rest = (SlabValue-1)-z1;
            x1 += (int)(((Long)rest*(x2-x1))/delta);
            y1 += (int)(((Long)rest*(y2-y1))/delta);
            z1 = SlabValue-1;
        }
        OutCode3(code1,x1,y1,z1);
        if( Reject(code1,code2) )
            return;
    }
    DrawTwinVector(x1,y1,z1,x2,y2,z2,col,col);
}


void ClipTwinVector( int x1, int y1, int z1,
                     int x2, int y2, int z2,
                     int col1, int col2 )
{
    register int xmid,ymid,zmid;
    register int code1,code2;

    if( col1!=col2 )
    {   OutCode3(code1,x1,y1,z1);
        OutCode3(code2,x2,y2,z2);
        if( !Reject(code1,code2) )
        {   if( !Accept(code1,code2) )
            {  xmid = (x1+x2)/2;
               ymid = (y1+y2)/2;
               zmid = (z1+z2)/2;
               ClipVector(x1,y1,z1,xmid,ymid,zmid,col1);
               ClipVector(xmid,ymid,zmid,x2,y2,z2,col2);
            } else
               DrawTwinVector(x1,y1,z1,x2,y2,z2,col1,col2);
        }
    } else ClipVector(x1,y1,z1,x2,y2,z2,col1);
}


/*==================================*/
/*  Monochrome Depth-Cued Vectors!  */
/*==================================*/

void DrawTwinVector2( int x1, int y1, int z1,
                      int x2, int y2, int z2,
                      int col1, int col2 )
{
    register int inten;
    register int midz;

    midz = ((z1+z2)/2)+ImageRadius-ZOffset;
    if( midz >= ImageSize )
    {   inten = ColourMask;
    } else if( midz > 0 )
    {   inten = (ColourDepth*midz)/ImageSize;
    } else inten = 0;
    DrawTwinLine(x1,y1,z1,x2,y2,z2,col1+inten,col2+inten);
}


void ClipTwinVector2( int x1, int y1, int z1,
                      int x2, int y2, int z2,
                      int col1, int col2 )
{
    register int inten;
    register int midz;

    midz = ((z1+z2)/2)+ImageRadius-ZOffset;
    if( midz >= ImageSize )
    {   inten = ColourMask;
    } else if( midz > 0 )
    {   inten = (ColourDepth*midz)/ImageSize;
    } else inten = 0;
    ClipTwinLine(x1,y1,z1,x2,y2,z2,col1+inten,col2+inten);
}


void ClipDashVector( int x1, int y1, int z1,
                     int x2, int y2, int z2,
                     int col1, int col2 )
{
    register Long offset;
    register Pixel __huge *fptr;
    register short __huge *dptr;
    register int ix,iy,iz,ic;
    register int dx,dy,dz,dc;
    register int crate, cerr;
    register int zrate, zerr;
    register int ystep,err;
    register int col, mid;
    register int c1, c2;
    register int count;

    if( (x1==x2) && (y1==y2) ) 
        return;

    /* Reject(OutCode3(x1,y1,z1),OutCode3(x2,y2,z2)) */
    if( (x1<0) && (x2<0) ) return;
    if( (y1<0) && (y2<0) ) return;
    if( (x1>=View.xmax) && (x2>=View.xmax) ) return;
    if( (y1>=View.ymax) && (y2>=View.ymax) ) return;
    if( UseSlabPlane && (z1>=SlabValue) && (z2>=SlabValue) )
        return;

    c1 = (ColourDepth*(z1+ImageRadius-ZOffset))/ImageSize;
    c2 = (ColourDepth*(z2+ImageRadius-ZOffset))/ImageSize;

    dx = x2 - x1;  dy = y2 - y1;
    dz = z2 - z1;  dc = c2 - c1;

    offset = (Long)y1*View.yskip + x1;
    fptr = View.fbuf+offset;
    dptr = View.dbuf+offset;
    count = 0;

    ystep = View.yskip;
    ix = iy = iz = ic = 1;
    if( dy<0 ) { dy = -dy; iy = -1; ystep = -ystep; }
    if( dx<0 ) { dx = -dx; ix = -1; }
    if( dz<0 ) { dz = -dz; iz = -1; }
    if( dc<0 ) { dc = -dc; ic = -1; }

    if( dx>dy )
    {   if( x2<x1 )
        {   mid = col1;
            col1 = col2;
            col2 = mid;
        }
        if( dz >= dx )
        {   zrate = dz/dx;
            dz -= dx*zrate;
            if( iz < 0 )
                zrate = -zrate;
        } else zrate = 0;

        if( dc >= dx )
        {   crate = dc/dx;
            dc -= dx*crate;
        } else crate = 0;

        err = zerr = cerr = -(dx>>1);
        mid = (x1+x2)/2;

        while( x1!=x2 )
        {   if( XValid(x1) && YValid(y1) && ZValid(z1) )
            {   if( count<2 )
                {   col = (x1<mid)? col1 : col2;
                    SETPIXEL(dptr,fptr,z1,Lut[col+c1]);
                    count++;
                } else if( count==3 )
                {   count = 0;
                } else count++;
            }

            if( (err+=dy)>0 )
            {   err -= dx;
                fptr+=ystep;
                dptr+=ystep;
                y1+=iy;
            }

            if( (zerr+=dz)>0 )
            {   zerr -= dx;
                z1 += iz;
            }

            if( (cerr+=dc)>0 )
            {   cerr -= dx;
                c1 += ic;
            }

            fptr+=ix; dptr+=ix; x1+=ix;
            z1 += zrate;   c1 += crate;
        }
    } else
    {   if( y1>y2 )
        {   mid = col1;
            col1 = col2;
            col2 = mid;
        }

        if( dz >= dy )
        {   zrate = dz/dy;
            dz -= dy*zrate;
            if( iz < 0 )
                zrate = -zrate;
        } else zrate = 0;

        if( dc >= dy )
        {   crate = dc/dy;
            dc -= dy*crate;
        } else crate = 0;

        err = zerr = cerr = -(dy>>1);
        mid = (y1+y2)/2;

        
        while( y1!=y2 )
        {   if( XValid(x1) && YValid(y1) && ZValid(z1) )
            {   if( count<2 )
                {   col = (y1<mid)? col1 : col2;
                    SETPIXEL(dptr,fptr,z1,Lut[col+c1]);
                    count++;
                } else if( count==3 )
                {   count = 0;
                } else count++;
            }

            if( (err+=dx)>0 )
            {   err-=dy;
                fptr+=ix;
                dptr+=ix;
                x1+=ix;
            }

            if( (zerr+=dz)>0 )
            {   zerr -= dy;
                z1 += iz;
            }

            if( (cerr+=dc)>0 )
            {   cerr -= dy;
                c1 += ic;
            }

            fptr+=ystep; dptr+=ystep; y1+=iy;
            z1 += zrate;   c1 += crate;
        }
    }
}


/* SplineCount is either 1, 2, 3, 4, 5 or 9! */
void StrandRibbon( Knot __far *src, Knot __far *dst, int col1, int col2 )
{
    register int hsx, hsy, hsz;
    register int hdx, hdy, hdz;
    register int qsx, qsy, qsz;
    register int qdx, qdy, qdz;
    register int col;

    if( SplineCount != 4 )
    {   if( SplineCount == 1 ) 
        {   ClipVector( src->px, src->py, src->pz,
                        dst->px, dst->py, dst->pz, col2 );
            return;
        } else if( SplineCount != 2 )
            ClipVector( src->px, src->py, src->pz,
                        dst->px, dst->py, dst->pz, col1 );

        ClipVector( src->px+src->wx, src->py+src->wy, src->pz+src->wz,
                    dst->px+dst->wx, dst->py+dst->wy, dst->pz+dst->wz, col2 );
        ClipVector( src->px-src->wx, src->py-src->wy, src->pz-src->wz,
                    dst->px-dst->wx, dst->py-dst->wy, dst->pz-dst->wz, col2 );
        if( SplineCount<=3 ) return;

        hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;

        ClipVector( src->px+hsx, src->py+hsy, src->pz+hsz,
                    dst->px+hdx, dst->py+hdy, dst->pz+hdz, col1 );
        ClipVector( src->px-hsx, src->py-hsy, src->pz-hsz,
                    dst->px-hdx, dst->py-hdy, dst->pz-hdz, col1 );
        if( SplineCount==5 ) 
            return;
        col = col1;
    } else /* SplineCount == 4 */
    {   hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;
        col = col2;
    }

    qsx = hsx/2;  qsy = hsy/2;  qsz = hsz/2;
    qdx = hdx/2;  qdy = hdy/2;  qdz = hdz/2;

    ClipVector( src->px+hsx+qsx, src->py+hsy+qsy, src->pz+hsz+qsz,
                dst->px+hdx+qdx, dst->py+hdy+qdy, dst->pz+hdz+qdz, col );
    ClipVector( src->px+hsx-qsx, src->py+hsy-qsy, src->pz+hsz-qsz,
                dst->px+hdx-qdx, dst->py+hdy-qdy, dst->pz+hdz-qdz, col1 );
    ClipVector( src->px-hsx+qsx, src->py-hsy+qsy, src->pz-hsz+qsz,
                dst->px-hdx+qdx, dst->py-hdy+qdy, dst->pz-hdz+qdz, col1 );
    ClipVector( src->px-hsx-qsx, src->py-hsy-qsy, src->pz-hsz-qsz,
                dst->px-hdx-qdx, dst->py-hdy-qdy, dst->pz-hdz-qdz, col );
}


void DashRibbon( Knot __far *src, Knot __far *dst, int col1, int col2 )
{
    register int hsx, hsy, hsz;
    register int hdx, hdy, hdz;
    register int qsx, qsy, qsz;
    register int qdx, qdy, qdz;
    register int col;

    if( SplineCount != 4 )
    {   if( SplineCount == 1 ) 
        {   ClipDashVector( src->px, src->py, src->pz,
                            dst->px, dst->py, dst->pz, col2, col2 );
            return;
        } else if( SplineCount != 2 )
            ClipDashVector( src->px, src->py, src->pz,
                            dst->px, dst->py, dst->pz, col1, col1 );

        ClipDashVector(src->px+src->wx,src->py+src->wy,src->pz+src->wz,
                       dst->px+dst->wx,dst->py+dst->wy,dst->pz+dst->wz,
                       col2,col2);
        ClipDashVector(src->px-src->wx,src->py-src->wy,src->pz-src->wz,
                       dst->px-dst->wx,dst->py-dst->wy,dst->pz-dst->wz,
                       col2,col2);
        if( SplineCount<=3 ) return;

        hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;

        ClipDashVector( src->px+hsx, src->py+hsy, src->pz+hsz,
                        dst->px+hdx, dst->py+hdy, dst->pz+hdz, col1, col1 );
        ClipDashVector( src->px-hsx, src->py-hsy, src->pz-hsz,
                        dst->px-hdx, dst->py-hdy, dst->pz-hdz, col1, col1 );
        if( SplineCount==5 ) 
            return;
        col = col1;
    } else /* SplineCount == 4 */
    {   hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;
        col = col2;
    }

    qsx = hsx/2;  qsy = hsy/2;  qsz = hsz/2;
    qdx = hdx/2;  qdy = hdy/2;  qdz = hdz/2;

    ClipDashVector(src->px+hsx+qsx,src->py+hsy+qsy,src->pz+hsz+qsz,
                   dst->px+hdx+qdx,dst->py+hdy+qdy,dst->pz+hdz+qdz,col,col);
    ClipDashVector(src->px+hsx-qsx,src->py+hsy-qsy,src->pz+hsz-qsz,
                   dst->px+hdx-qdx,dst->py+hdy-qdy,dst->pz+hdz-qdz,col1,col1);
    ClipDashVector(src->px-hsx+qsx,src->py-hsy+qsy,src->pz-hsz+qsz,
                   dst->px-hdx+qdx,dst->py-hdy+qdy,dst->pz-hdz+qdz,col1,col1);
    ClipDashVector(src->px-hsx-qsx,src->py-hsy-qsy,src->pz-hsz-qsz,
                   dst->px-hdx-qdx,dst->py-hdy-qdy,dst->pz-hdz-qdz,col,col);
}


void DrawPolygonOutline( Poly *p )
{
    register int i;

    for( i=0; i<p->count-1; i++ )
        DrawTwinLine( p->v[i].x, p->v[i].y, p->v[i].z, 
                      p->v[i+1].x, p->v[i+1].y, p->v[i+1].z,
                      p->v[i].inten, p->v[i].inten);
    DrawTwinLine( p->v[i].x, p->v[i].y, p->v[i].z,
                  p->v[0].x, p->v[0].y, p->v[0].z,
                  p->v[i].inten, p->v[i].inten);
}


void ClipPolygonOutline( Poly *p )
{
    register int i;

    for( i=0; i<p->count-1; i++ )
        ClipLine( p->v[i].x, p->v[i].y, p->v[i].z, 
                  p->v[i+1].x, p->v[i+1].y, p->v[i+1].z,
                  p->v[i].inten);
    ClipLine( p->v[i].x, p->v[i].y, p->v[i].z,
              p->v[0].x, p->v[0].y, p->v[0].z,
              p->v[i].inten);
}


/*====================================*/
/*  Flat (Constant) Shaded Triangles  */
/*====================================*/

static void DrawFlatTriFlatTop( Vert *li, Vert *ri, Vert *bot, int col )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int x,xmin,xmax;
    register int y,ymin,dy;
    register int offset;
    register int z,dz;

    if( li->x >= ri->x )
        return;

    dz = ((ri->z-li->z)<<16)/(ri->x-li->x);

    y = li->y;
    ymin = bot->y;
    dy = ymin - y;

    lft.dx = ((bot->x - li->x)<<16)/dy;
    lft.dz = ((bot->z - li->z)<<16)/dy;
    rgt.dx = ((bot->x - ri->x)<<16)/dy;
    lft.x = li->x<<16;
    lft.z = li->z<<16;
    rgt.x = ri->x<<16;

    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        /* dz = (rgt.z-lft.z)/(xmax-xmin); */
        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ )
        {   if( (int)(z>>16) > *dptr )
            {   fbase[x] = Lut[col];
                *dptr = (int)(z>>16);
            }
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void DrawFlatTriFlatBot( Vert *top, Vert *li, Vert *ri, int col )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int x,xmin,xmax;
    register int y,ymin,dy;
    register int offset;
    register int z,dz;

    if( li->x >= ri->x )
        return;
    if( li->y == top->y+1 )
        return;

    dz = ((ri->z-li->z)<<16)/(ri->x-li->x);

    ymin = li->y;
    dy = ymin - top->y;
    lft.dx = ((li->x - top->x)<<16)/dy;
    lft.dz = ((li->z - top->z)<<16)/dy;
    rgt.dx = ((ri->x - top->x)<<16)/dy;
    rgt.x = lft.x = top->x<<16;
    lft.z = top->z<<16;

    /* Top vertex never drawn! */
    lft.x += lft.dx;  rgt.x += rgt.dx;
    lft.z += lft.dz;

    y = top->y+1;
    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        /* dz = (rgt.z-lft.z)/(xmax-xmin); */
        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[col];
                *dptr = (int)(z>>16);
            }
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void DrawFlatTriLeft( Vert *top, Vert *li, Vert *bot, int col )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int offset;
    register int dx,dy,ry,ly;
    register int xmin,xmax;
    register int x,y,ymin;
    register int z,dz;

    /* Determine scanline gradients */
    ly = li->y - top->y;
    ry = bot->y - top->y;

    dx = (li->x-top->x)*ry - (bot->x-top->x)*ly;
    if( dx )
    {   z = (li->z-top->z)*ry - (bot->z-top->z)*ly;
        dz = (z<<16)/dx;
    } else dz = 0;

    lft.dx = ((li->x - top->x)<<16)/ly;
    rgt.dx = ((bot->x - top->x)<<16)/ry;
    rgt.dz = ((bot->z - top->z)<<16)/ry;

    rgt.x = lft.x = top->x<<16;
    rgt.z = top->z<<16;

    /* Top vertex never drawn! */
    lft.x += lft.dx;
    rgt.x += rgt.dx;
    rgt.z += rgt.dz;

    y = top->y+1;
    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    ymin = li->y;
    while( y<ymin )
    {   xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        z = rgt.z;

        dptr = dbase+(xmax-1);
        for( x=xmax-1; x>=xmin; x-- )
        {   if( (int)(z>>16) > *dptr )
            {   fbase[x] = Lut[col];
                *dptr = (int)(z>>16);
            }
            z -= dz;
            dptr--;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        rgt.z += rgt.dz;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }

    dy = bot->y - li->y;
    lft.dx = ((bot->x - li->x)<<16)/dy;
    lft.x = li->x<<16;

    ymin = bot->y;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        z = rgt.z;

        dptr = dbase+(xmax-1);
        for( x=xmax-1; x>=xmin; x-- ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[col];
                *dptr = (int)(z>>16);
            }
            z -= dz;
            dptr--;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        rgt.z += rgt.dz;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void DrawFlatTriRight( Vert *top, Vert *ri, Vert *bot, int col )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int offset;
    register int dx,dy,ry,ly;
    register int xmin,xmax;
    register int x,y,ymin;
    register int z,dz;

    /* Determine scanline gradients */
    ly = bot->y - top->y;
    ry = ri->y - top->y;

    dx = (bot->x-top->x)*ry - (ri->x-top->x)*ly;
    if( dx )
    {   z = (bot->z-top->z)*ry - (ri->z-top->z)*ly;
        dz = (z<<16)/dx;
    } else dz = 0;

    lft.dz = ((bot->z - top->z)<<16)/ly;
    lft.dx = ((bot->x - top->x)<<16)/ly;
    rgt.dx = ((ri->x - top->x)<<16)/ry;

    lft.z = top->z<<16;
    rgt.x = lft.x = top->x<<16;

    /* Top vertex never drawn! */
    lft.x += lft.dx;
    rgt.x += rgt.dx;
    lft.z += lft.dz;

    y = top->y+1;
    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    ymin = ri->y;
    while( y<ymin )
    {   xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ )
        {   if( (int)(z>>16) > *dptr )
            {   fbase[x] = Lut[col];
                *dptr = (int)(z>>16);
            }
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }

    dy = bot->y - ri->y;
    rgt.dx = ((bot->x - ri->x)<<16)/dy;
    rgt.x = ri->x<<16;

    ymin = bot->y;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[col];
                *dptr = (int)(z>>16);
            }
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


void DrawFlatTriangle( Vert *p, Vert *q, Vert *r, int col )
{
    if( p->y < q->y )
    {   if( p->y < r->y )
        {   if( q->y < r->y )
            {   DrawFlatTriRight(p,q,r,col);
            } else if( q->y > r->y )
            {   DrawFlatTriLeft(p,r,q,col);
            } else DrawFlatTriFlatBot(p,r,q,col);
        } else  if( p->y > r->y )
        {   DrawFlatTriRight(r,p,q,col);
        } else DrawFlatTriFlatTop(r,p,q,col);

    } else if( p->y > q->y )
    {   if( q->y < r->y )
        {   if( p->y < r->y )
            {   DrawFlatTriLeft(q,p,r,col);
            } else if( p->y > r->y )
            {   DrawFlatTriRight(q,r,p,col);
            } else DrawFlatTriFlatBot(q,p,r,col);
         } else if( q->y > r->y )
         {  DrawFlatTriLeft(r,q,p,col);
         } else DrawFlatTriFlatTop(q,r,p,col);

    } else /* p->y == q->y */
         if( p->y < r->y )
         {   DrawFlatTriFlatTop(p,q,r,col);
         } else if( p->y > r->y )
             DrawFlatTriFlatBot(r,q,p,col);
}


/*=====================================*/
/*  Gouraud (Smooth) Shaded Triangles  */
/*=====================================*/

static void DrawTriFlatTop( Vert *li, Vert *ri, Vert *bot )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int x,xmin,xmax;
    register int y,ymin,dy;
    register int z,inten;
    register int offset;
    register int dz,di;

    if( li->x >= ri->x )
        return;

    di = ((ri->inten-li->inten)<<16)/(ri->x-li->x);
    dz = ((ri->z-li->z)<<16)/(ri->x-li->x);

    y = li->y;
    ymin = bot->y;
    dy = ymin - y;

    lft.di = ((bot->inten - li->inten)<<16)/dy;
    lft.dx = ((bot->x - li->x)<<16)/dy;
    lft.dz = ((bot->z - li->z)<<16)/dy;
    rgt.dx = ((bot->x - ri->x)<<16)/dy;
    lft.i = li->inten<<16;
    lft.x = li->x<<16;
    lft.z = li->z<<16;
    rgt.x = ri->x<<16;

    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        inten = lft.i;
        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void DrawTriFlatBot( Vert *top, Vert *li, Vert *ri )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int x,xmin,xmax;
    register int y,ymin,dy;
    register int z,inten;
    register int offset;
    register int dz,di;

    if( li->x >= ri->x )
        return;
    if( li->y == top->y+1 )
        return;

    di = ((ri->inten-li->inten)<<16)/(ri->x-li->x);
    dz = ((ri->z-li->z)<<16)/(ri->x-li->x);

    ymin = li->y;
    dy = ymin - top->y;
    lft.di = ((li->inten - top->inten)<<16)/dy;
    lft.dx = ((li->x - top->x)<<16)/dy;
    lft.dz = ((li->z - top->z)<<16)/dy;
    rgt.dx = ((ri->x - top->x)<<16)/dy;

    rgt.x = lft.x = top->x<<16;
    lft.i = top->inten<<16;
    lft.z = top->z<<16;

    /* Top vertex never drawn! */
    lft.x += lft.dx;  rgt.x += rgt.dx;
    lft.z += lft.dz;  lft.i += lft.di;

    y = top->y+1;
    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        inten = lft.i;
        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void DrawTriLeft( Vert *top, Vert *li, Vert *bot )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int dx,dy,ry,ly;
    register int xmin,xmax;
    register int x,y,ymin;
    register int z,inten;
    register int offset;
    register int dz,di;

    /* Determine scanline gradients */
    ly = li->y - top->y;
    ry = bot->y - top->y;

    dx = (li->x-top->x)*ry - (bot->x-top->x)*ly;
    if( dx )
    {   inten = (li->inten-top->inten)*ry - (bot->inten-top->inten)*ly;
        z = (li->z-top->z)*ry - (bot->z-top->z)*ly;
        di = (inten<<16)/dx;
        dz = (z<<16)/dx;
    } else dz = di = 0;

    rgt.di = ((bot->inten - top->inten)<<16)/ry;
    rgt.dz = ((bot->z - top->z)<<16)/ry;
    rgt.dx = ((bot->x - top->x)<<16)/ry;
    lft.dx = ((li->x - top->x)<<16)/ly;

    rgt.x = lft.x = top->x<<16;
    rgt.i = top->inten<<16;
    rgt.z = top->z<<16;

    /* Top vertex never drawn! */
    lft.x += lft.dx;
    rgt.x += rgt.dx;
    rgt.z += rgt.dz;
    rgt.i += rgt.di;

    y = top->y+1;
    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    ymin = li->y;
    while( y<ymin )
    {   xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        inten = rgt.i;
        z = rgt.z;

        dptr = dbase+(xmax-1);
        for( x=xmax-1; x>=xmin; x-- )
        {   if( (int)(z>>16) > *dptr )
            {   fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten -= di;
            z -= dz;
            dptr--;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        rgt.z += rgt.dz;
        rgt.i += rgt.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }

    dy = bot->y - li->y;
    lft.dx = ((bot->x - li->x)<<16)/dy;
    lft.x = li->x<<16;

    ymin = bot->y;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        inten = rgt.i;
        z = rgt.z;

        dptr = dbase+(xmax-1);
        for( x=xmax-1; x>=xmin; x-- ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten -= di;
            z -= dz;
            dptr--;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        rgt.z += rgt.dz;
        rgt.i += rgt.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void DrawTriRight( Vert *top, Vert *ri, Vert *bot )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int dx,dy,ry,ly;
    register int xmin,xmax;
    register int x,y,ymin;
    register int z,inten;
    register int offset;
    register int dz,di;

    /* Determine scanline gradients */
    ly = bot->y - top->y;
    ry = ri->y - top->y;

    dx = (bot->x-top->x)*ry - (ri->x-top->x)*ly;
    if( dx )
    {   inten = (bot->inten-top->inten)*ry - (ri->inten-top->inten)*ly;
        z = (bot->z-top->z)*ry - (ri->z-top->z)*ly;
        di = (inten<<16)/dx;
        dz = (z<<16)/dx;
    } else dz = di = 0;

    lft.di = ((bot->inten - top->inten)<<16)/ly;
    lft.dz = ((bot->z - top->z)<<16)/ly;
    lft.dx = ((bot->x - top->x)<<16)/ly;
    rgt.dx = ((ri->x - top->x)<<16)/ry;

    rgt.x = lft.x = top->x<<16;
    lft.i = top->inten<<16;
    lft.z = top->z<<16;

    /* Top vertex never drawn! */
    lft.x += lft.dx;
    rgt.x += rgt.dx;
    lft.z += lft.dz;
    lft.i += lft.di;

    y = top->y+1;
    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    ymin = ri->y;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        inten = lft.i;
        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }

    dy = bot->y - ri->y;
    rgt.dx = ((bot->x - ri->x)<<16)/dy;
    rgt.x = ri->x<<16;

    ymin = bot->y;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        inten = lft.i;
        z = lft.z;

        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


void DrawTriangle( Vert *p, Vert *q, Vert *r )
{
    if( p->y < q->y )
    {   if( p->y < r->y )
        {   if( q->y < r->y )
            {   DrawTriRight(p,q,r);
            } else if( q->y > r->y )
            {   DrawTriLeft(p,r,q);
            } else DrawTriFlatBot(p,r,q);
        } else  if( p->y > r->y )
        {   DrawTriRight(r,p,q);
        } else DrawTriFlatTop(r,p,q);

    } else if( p->y > q->y )
    {   if( q->y < r->y )
        {   if( p->y < r->y )
            {   DrawTriLeft(q,p,r);
            } else if( p->y > r->y )
            {   DrawTriRight(q,r,p);
            } else DrawTriFlatBot(q,p,r);
         } else if( q->y > r->y )
         {  DrawTriLeft(r,q,p);
         } else DrawTriFlatTop(q,r,p);

    } else /* p->y == q->y */
         if( p->y < r->y )
         {   DrawTriFlatTop(p,q,r);
         } else if( p->y > r->y )
             DrawTriFlatBot(r,q,p);
}


static void ClipTriFlatTop( Vert *li, Vert *ri, Vert *bot )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int x,xmin,xmax;
    register int y,ymin,dy;
    register int z,inten;
    register int offset;
    register int dz,di;

    if( li->x >= ri->x )
        return;

    di = ((ri->inten-li->inten)<<16)/(ri->x-li->x);
    dz = ((ri->z-li->z)<<16)/(ri->x-li->x);

    y = li->y;
    ymin = bot->y;
    dy = ymin - y;

    lft.di = ((bot->inten - li->inten)<<16)/dy;
    lft.dx = ((bot->x - li->x)<<16)/dy;
    lft.dz = ((bot->z - li->z)<<16)/dy;
    rgt.dx = ((bot->x - ri->x)<<16)/dy;
    lft.i = li->inten<<16;
    lft.x = li->x<<16;
    lft.z = li->z<<16;
    rgt.x = ri->x<<16;

    if( y < 0 )
    {
      rgt.x -= y*rgt.dx;
      lft.x -= y*lft.dx;
      lft.i -= y*lft.di;
      lft.z -= y*lft.dz;
      y = 0;
    }

    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    if( ymin > View.ymax )
      ymin = View.ymax;

    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        if (xmin < 0) {
            inten = lft.i - xmin*di;
            z = lft.z - xmin*dz;
            xmin = 0;
        } else {
            inten = lft.i;
            z = lft.z;
        }

        if (xmax > View.xmax)
          xmax = View.xmax;
        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void ClipTriFlatBot( Vert *top, Vert *li, Vert *ri )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int x,xmin,xmax;
    register int y,ymin,dy;
    register int z,inten;
    register int offset;
    register int dz,di;

    if( li->x >= ri->x )
        return;
    if( li->y == top->y+1 )
        return;

    di = ((ri->inten-li->inten)<<16)/(ri->x-li->x);
    dz = ((ri->z-li->z)<<16)/(ri->x-li->x);

    ymin = li->y;
    dy = ymin - top->y;
    lft.di = ((li->inten - top->inten)<<16)/dy;
    lft.dx = ((li->x - top->x)<<16)/dy;
    lft.dz = ((li->z - top->z)<<16)/dy;
    rgt.dx = ((ri->x - top->x)<<16)/dy;

    rgt.x = lft.x = top->x<<16;
    lft.i = top->inten<<16;
    lft.z = top->z<<16;

    y = top->y;
    if (y < 0) {
        rgt.x -= y*rgt.dx;
        lft.x -= y*lft.dx;
        lft.z -= y*lft.dz;
        lft.i -= y*lft.di;
        y = 0;
    } else {
        /* Top vertex never drawn! */
        lft.x += lft.dx;  rgt.x += rgt.dx;
        lft.z += lft.dz;  lft.i += lft.di;
        y++;
    }

    offset = y*View.yskip;
    fbase = (Pixel*)View.fbuf+offset;
    dbase = View.dbuf+offset;

    if (ymin > View.ymax)
        ymin = View.ymax;

    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        if (xmin < 0) {
            inten = lft.i - xmin*di;
            z = lft.z - xmin*dz;
            xmin = 0;
        } else {
            inten = lft.i;
            z = lft.z;
        }

        if (xmax > View.xmax)
          xmax = View.xmax;
        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void ClipTriLeft( Vert *top, Vert *li, Vert *bot )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int dx,dy,ry,ly;
    register int xmin,xmax;
    register int x,y,ymin;
    register int z,inten;
    register int offset;
    register int dz,di;

    /* Determine scanline gradients */
    ly = li->y - top->y;
    ry = bot->y - top->y;

    dx = (li->x-top->x)*ry - (bot->x-top->x)*ly;
    if( dx )
    {   inten = (li->inten-top->inten)*ry - (bot->inten-top->inten)*ly;
        z = (li->z-top->z)*ry - (bot->z-top->z)*ly;
        di = (inten<<16)/dx;
        dz = (z<<16)/dx;
    } else dz = di = 0;

    rgt.di = ((bot->inten - top->inten)<<16)/ry;
    rgt.dz = ((bot->z - top->z)<<16)/ry;
    rgt.dx = ((bot->x - top->x)<<16)/ry;

    rgt.x = lft.x = top->x<<16;
    rgt.i = top->inten<<16;
    rgt.z = top->z<<16;

    if( li->y < 0 ) {
        rgt.x += ly*rgt.dx;
        rgt.i += ly*rgt.di;
        rgt.z += ly*rgt.dz;
        y = li->y;

        /* Avoid compiler warnings */
        fbase = (Pixel*)0;
        dbase = (short*)0;
    } else {
        lft.dx = ((li->x - top->x)<<16)/ly;
        lft.x = rgt.x;
        y = top->y;
        if( y < 0 ) {
            lft.x -= y*lft.dx;
            rgt.x -= y*rgt.dx;
            rgt.i -= y*rgt.di;
            rgt.z -= y*rgt.dz;
            y = 0;

            fbase = (Pixel*)View.fbuf;
            dbase = View.dbuf;
        } else {
            /* Top vertex never drawn! */
            lft.x += lft.dx;
            rgt.x += rgt.dx;
            rgt.z += rgt.dz;
            rgt.i += rgt.di;
            y++;

            offset = y*View.yskip;
            fbase = (Pixel*)View.fbuf+offset;
            dbase = View.dbuf+offset;
        }
    }

    ymin = li->y;
    if( ymin > View.ymax )
        ymin = View.ymax;
    while( y<ymin )
    {   xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        if( xmax > View.xmax ) {
            offset = xmax - View.xmax;
            inten = rgt.i - offset*di;
            z = rgt.z - offset*dz;
            xmax = View.xmax;
        } else {
            inten = rgt.i;
            z = rgt.z;
        }

        if( xmin < 0 )
            xmin = 0;
        dptr = dbase+(xmax-1);
        for( x=xmax-1; x>=xmin; x-- )
        {   if( (int)(z>>16) > *dptr )
            {   fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten -= di;
            z -= dz;
            dptr--;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        rgt.z += rgt.dz;
        rgt.i += rgt.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }

    dy = bot->y - li->y;
    lft.dx = ((bot->x - li->x)<<16)/dy;
    lft.x = li->x<<16;

    if( y < 0 ) {
        lft.x -= y*lft.dx;
        rgt.x -= y*rgt.dx;
        rgt.i -= y*rgt.di;
        rgt.z -= y*rgt.dz;
        y = 0;

        fbase = (Pixel*)View.fbuf;
        dbase = View.dbuf;
    } else {
        offset = y*View.yskip;
        fbase = (Pixel*)View.fbuf+offset;
        dbase = View.dbuf+offset;
    }

    ymin = bot->y;
    if( ymin > View.ymax )
        ymin = View.ymax;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        if( xmax > View.xmax ) {
            offset = xmax - View.xmax;
            inten = rgt.i - offset*di;
            z = rgt.z - offset*dz;
            xmax = View.xmax;
        } else {
            inten = rgt.i;
            z = rgt.z;
        }

        if( xmin < 0 )
            xmin = 0;
        dptr = dbase+(xmax-1);
        for( x=xmax-1; x>=xmin; x-- ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten -= di;
            z -= dz;
            dptr--;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        rgt.z += rgt.dz;
        rgt.i += rgt.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


static void ClipTriRight( Vert *top, Vert *ri, Vert *bot )
{
    static Edge lft, rgt;
    register Pixel *fbase;
    register short *dbase;
    register short *dptr;
    register int dx,dy,ry,ly;
    register int xmin,xmax;
    register int x,y,ymin;
    register int z,inten;
    register int offset;
    register int dz,di;

    /* Determine scanline gradients */
    ly = bot->y - top->y;
    ry = ri->y - top->y;

    dx = (bot->x-top->x)*ry - (ri->x-top->x)*ly;
    if( dx )
    {   inten = (bot->inten-top->inten)*ry - (ri->inten-top->inten)*ly;
        z = (bot->z-top->z)*ry - (ri->z-top->z)*ly;
        di = (inten<<16)/dx;
        dz = (z<<16)/dx;
    } else dz = di = 0;

    lft.di = ((bot->inten - top->inten)<<16)/ly;
    lft.dz = ((bot->z - top->z)<<16)/ly;
    lft.dx = ((bot->x - top->x)<<16)/ly;

    lft.x = top->x<<16;
    lft.i = top->inten<<16;
    lft.z = top->z<<16;

    if( ri->y < 0 ) {
        lft.x += ry*lft.dx;
        lft.i += ry*lft.di;
        lft.z += ry*lft.dz;
        y = ri->y;

        /* Avoid compiler warnings */
        fbase = (Pixel*)0;
        dbase = View.dbuf;
    } else {
        rgt.dx = ((ri->x - top->x)<<16)/ry;
        rgt.x = lft.x;
        y = top->y;
        if( y < 0 ) {
            rgt.x -= y*rgt.dx;
            lft.x -= y*lft.dx;
            lft.i -= y*lft.di;
            lft.z -= y*lft.dz;
            y = 0;

            fbase = (Pixel*)View.fbuf;
            dbase = View.dbuf;
        } else {
            /* Top vertex never drawn! */
            rgt.x += rgt.dx;
            lft.x += lft.dx;
            lft.z += lft.dz;
            lft.i += lft.di;
            y++;

            offset = y*View.yskip;
            fbase = (Pixel*)View.fbuf+offset;
            dbase = View.dbuf+offset;
        }
    }

    ymin = ri->y;
    if( ymin > View.ymax )
        ymin = View.ymax;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        if( xmin < 0 ) {
            inten = lft.i - xmin*di;
            z = lft.z - xmin*dz;
            xmin = 0;
        } else {
            inten = lft.i;
            z = lft.z;
        }

        if( xmax > View.xmax )
            xmax = View.xmax;
        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }

    dy = bot->y - ri->y;
    rgt.dx = ((bot->x - ri->x)<<16)/dy;
    rgt.x = ri->x<<16;

    if( y < 0 ) {
        rgt.x -= y*rgt.dx;
        lft.x -= y*lft.dx;
        lft.i -= y*lft.di;
        lft.z -= y*lft.dz;
        y = 0;

        fbase = (Pixel*)View.fbuf;
        dbase = View.dbuf;
    } else {
        offset = y*View.yskip;
        fbase = (Pixel*)View.fbuf+offset;
        dbase = View.dbuf+offset;
    }

    ymin = bot->y;
    if( ymin > View.ymax )
        ymin = View.ymax;
    while( y < ymin ) {
        xmax = (int)(rgt.x>>16);
        xmin = (int)(lft.x>>16);

        if( xmin < 0 ) {
            inten = lft.i - xmin*di;
            z = lft.z - xmin*dz;
            xmin = 0;
        } else {
            inten = lft.i;
            z = lft.z;
        }

        if( xmax > View.xmax )
            xmax = View.xmax;
        dptr = dbase+xmin;
        for( x=xmin; x<xmax; x++ ) {
            if( (int)(z>>16) > *dptr ) {
                fbase[x] = Lut[(int)(inten>>16)];
                *dptr = (int)(z>>16);
            }
            inten += di;
            z += dz;
            dptr++;
        }

        lft.x += lft.dx;
        rgt.x += rgt.dx;
        lft.z += lft.dz;
        lft.i += lft.di;
        dbase += View.yskip;
        fbase += View.yskip;
        y++;
    }
}


void ClipTriangle( Vert *p, Vert *q, Vert *r )
{
    register int ocp, ocq, ocr;

    if( UseSlabPlane )
        if( p->z >= SlabValue ||
            q->z >= SlabValue ||
            r->z >= SlabValue )
            return;

    OutCode2(ocp,p->x,p->y);
    OutCode2(ocq,q->x,q->y);
    OutCode2(ocr,r->x,r->y);

    if( (ocp|ocq|ocr) == 0 )
    {   DrawTriangle(p,q,r);
        return;
    }

    if( (ocp&ocq&ocr) != 0 )
        return;

    if( p->y < q->y )
    {   if( p->y < r->y )
        {   if( q->y < r->y )
            {   ClipTriRight(p,q,r);
            } else if( q->y > r->y )
            {   ClipTriLeft(p,r,q);
            } else ClipTriFlatBot(p,r,q);
        } else  if( p->y > r->y )
        {   ClipTriRight(r,p,q);
        } else ClipTriFlatTop(r,p,q);

    } else if( p->y > q->y )
    {   if( q->y < r->y )
        {   if( p->y < r->y )
            {   ClipTriLeft(q,p,r);
            } else if( p->y > r->y )
            {   ClipTriRight(q,r,p);
            } else ClipTriFlatBot(q,p,r);
         } else if( q->y > r->y )
         {  ClipTriLeft(r,q,p);
         } else ClipTriFlatTop(q,r,p);

    } else /* p->y == q->y */
         if( p->y < r->y )
         {   ClipTriFlatTop(p,q,r);
         } else if( p->y > r->y )
             ClipTriFlatBot(r,q,p);
}


/*============================*/
/*  Convex Polygon Rendering  */
/*============================*/

void DrawFlatPolygon( Poly *p )
{
    static Edge lft, rgt;
    register Edge *pmin, *pmax;
    register Pixel __huge *fbase;
    register short __huge *dbase;
    register short __huge *dptr;
    register Long offset;

    register Long z,dz;
    register int ri,li,ry,ly;
    register int xmin,xmax;
    register int dy,ymin;
    register int top,rem;
    register int x,y,i;

    /* Find top vertex */
    top = 0;  
    ymin = p->v[0].y;
    for( i=1; i<p->count; i++ )
       if( p->v[i].y < ymin )
       {   ymin = p->v[i].y;
           top = i;
       }

    rem = p->count;
    ly = ry = y = ymin;
    li = ri = top;

    offset = (Long)y*View.yskip;
    fbase = View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( rem )
    {   while( ly<=y && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ly;
                lft.dx = (((Long)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((Long)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.x = ((Long)p->v[li].x)<<16;
                lft.z = ((Long)p->v[li].z)<<16;
            }
            ly = p->v[i].y;
            rem--;  li = i;
        }

        while( ry<=y && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ry;
                rgt.dx = (((Long)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((Long)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.x = ((Long)p->v[ri].x)<<16;
                rgt.z = ((Long)p->v[ri].z)<<16;
            }
            ry = p->v[i].y;
            rem--; ri = i;
        }


        ymin = MinFun(ly,ry);
        
        while( y<ymin )
        {   if( lft.x < rgt.x )
            {   pmin = &lft;
                pmax = &rgt;
            } else
            {   pmin = &rgt;
                pmax = &lft;
            }

            xmax = (int)(pmax->x>>16)+1;
            xmin = (int)(pmin->x>>16);

            dz = (Long)((pmax->z-pmin->z)/(xmax-xmin));
            z = pmin->z;

            dptr = dbase+xmin;
            for( x=xmin; x<xmax; x++ )
            {   if( (int)(z>>16) > *dptr )
                {   fbase[x] = Lut[p->v[0].inten];
                    *dptr = (int)(z>>16);
                }
                z += dz;
                dptr++;
            }

            lft.x += lft.dx;  rgt.x += rgt.dx;
            lft.z += lft.dz;  rgt.z += rgt.dz;
            dbase += View.yskip;
            fbase += View.yskip;
            y++;
        }
    }
}


void ClipFlatPolygon( Poly *p )
{
    static Edge lft, rgt;
    register Edge *pmin, *pmax;
    register Pixel __huge *fbase;
    register short __huge *dbase;
    register short __huge *dptr;
    register Long offset;

    register Long z,dz;
    register int ri,li,ry,ly;
    register int xmin,xmax;
    register int dy,ymin;
    register int top,rem;
    register int x,y,i;

    /* Reject Clip Polygon */
    if( UseSlabPlane )
        for( i=0; i<p->count; i++ )
            if( p->v[i].z >= SlabValue )
                return;

    /* Find top vertex */
    top = 0;  
    ymin = p->v[0].y;
    for( i=1; i<p->count; i++ )
       if( p->v[i].y < ymin )
       {   ymin = p->v[i].y;
           top = i;
       }

    rem = p->count;
    ly = ry = y = ymin;
    li = ri = top;

    if( y<0 )
    {   rem--;

        while( ly<=0 && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > 0 )
            {   dy = p->v[i].y - ly;
                lft.dx = (((Long)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((Long)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.x = ((Long)p->v[li].x)<<16;
                lft.z = ((Long)p->v[li].z)<<16;
            } else rem--;
            ly = p->v[i].y;
            li = i;
        }

        while( ry<=0 && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > 0 )
            {   dy = p->v[i].y - ry;
                rgt.dx = (((Long)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((Long)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.x = ((Long)p->v[ri].x)<<16;
                rgt.z = ((Long)p->v[ri].z)<<16;
            } else rem--;
            ry = p->v[i].y;
            ri = i;
        }

        fbase = View.fbuf;
        dbase = View.dbuf;
        y = 0;
    } else /* y >= 0 */
    {   offset = (Long)y*View.yskip;
        fbase = View.fbuf+offset;
        dbase = View.dbuf+offset;
    }

    while( rem )
    {   while( ly<=y && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ly;
                lft.dx = (((Long)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((Long)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.x = ((Long)p->v[li].x)<<16;
                lft.z = ((Long)p->v[li].z)<<16;
            }
            ly = p->v[i].y;
            rem--;  li = i;
        }

        while( ry<=y && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ry;
                rgt.dx = (((Long)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((Long)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.x = ((Long)p->v[ri].x)<<16;
                rgt.z = ((Long)p->v[ri].z)<<16;
            }
            ry = p->v[i].y;
            rem--; ri = i;
        }

        ymin = MinFun(ly,ry);
        if( ymin>View.ymax )
        {   ymin = View.ymax;
            rem = 0;
        }
        
        while( y<ymin )
        {   if( lft.x < rgt.x )
            {   pmin = &lft;
                pmax = &rgt;
            } else
            {   pmin = &rgt;
                pmax = &lft;
            }

            xmax = (int)(pmax->x>>16)+1;
            xmin = (int)(pmin->x>>16);

            if( (xmin<View.xmax) && (xmax>=0) )
            {   dz = (Long)((pmax->z-pmin->z)/(xmax-xmin));
                if( xmin<0 )
                {   z = pmin->z - xmin*dz;
                    xmin = 0;
                } else /* xmin >= 0 */
                    z = pmin->z;

                if( xmax>=View.xmax )
                    xmax = View.xmax;

                dptr = dbase+xmin;
                for( x=xmin; x<xmax; x++ )
                {   if( (int)(z>>16) > *dptr )
                    {   fbase[x] = Lut[p->v[0].inten];
                        *dptr = (int)(z>>16);
                    }
                    z += dz;
                    dptr++;
                }
            }

            lft.x += lft.dx;  rgt.x += rgt.dx;
            lft.z += lft.dz;  rgt.z += rgt.dz;
            dbase += View.yskip;
            fbase += View.yskip;
            y++;
        }
    }
}


void DrawPolygon( Poly *p )
{
    static Edge lft, rgt;
    register Edge *pmin, *pmax;
    register Pixel __huge *fbase;
    register short __huge *dbase;
    register short __huge *dptr;
    register Long offset;

    register Long dz,di;
    register Long z,inten;
    register int ri,li,ry,ly;
    register int xmin,xmax;
    register int dy,ymin;
    register int top,rem;
    register int x,y,i;

    /* Find top vertex */
    top = 0;  
    ymin = p->v[0].y;
    for( i=1; i<p->count; i++ )
       if( p->v[i].y < ymin )
       {   ymin = p->v[i].y;
           top = i;
       }

    rem = p->count;
    ly = ry = y = ymin;
    li = ri = top;

    offset = (Long)y*View.yskip;
    fbase = View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( rem )
    {   while( ly<=y && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ly;
                lft.di = (((Long)(p->v[i].inten - p->v[li].inten))<<16)/dy;
                lft.dx = (((Long)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((Long)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.i = ((Long)p->v[li].inten)<<16;
                lft.x = ((Long)p->v[li].x)<<16;
                lft.z = ((Long)p->v[li].z)<<16;
            }
            ly = p->v[i].y;
            rem--;  li = i;
        }

        while( ry<=y && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ry;
                rgt.di = (((Long)(p->v[i].inten - p->v[ri].inten))<<16)/dy;
                rgt.dx = (((Long)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((Long)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.i = ((Long)p->v[ri].inten)<<16;
                rgt.x = ((Long)p->v[ri].x)<<16;
                rgt.z = ((Long)p->v[ri].z)<<16;
            }
            ry = p->v[i].y;
            rem--; ri = i;
        }

        ymin = MinFun(ly,ry);
        
        while( y<ymin )
        {   if( lft.x < rgt.x )
            {   pmin = &lft;
                pmax = &rgt;
            } else
            {   pmin = &rgt;
                pmax = &lft;
            }

            xmax = (int)(pmax->x>>16)+1;
            xmin = (int)(pmin->x>>16);

            di = (Long)((pmax->i-pmin->i)/(xmax-xmin));
            dz = (Long)((pmax->z-pmin->z)/(xmax-xmin));
            inten = pmin->i;  
            z = pmin->z;

            dptr = dbase+xmin;
            for( x=xmin; x<xmax; x++ )
            {   if( (int)(z>>16) > *dptr )
                {   fbase[x] = Lut[(int)(inten>>16)];
                    *dptr = (int)(z>>16);
                }
                inten += di;
                z += dz;
                dptr++;
            }

            lft.x += lft.dx;  rgt.x += rgt.dx;
            lft.z += lft.dz;  rgt.z += rgt.dz;
            lft.i += lft.di;  rgt.i += rgt.di;
            dbase += View.yskip;
            fbase += View.yskip;
            y++;
        }
    }
}


void ClipPolygon( Poly *p )
{
    static Edge lft, rgt;
    register Edge *pmin, *pmax;
    register Pixel __huge *fbase;
    register short __huge *dbase;
    register short __huge *dptr;
    register Long offset;

    register Long dz,di;
    register Long z,inten;
    register int ri,li,ry,ly;
    register int xmin,xmax;
    register int dy,ymin;
    register int top,rem;
    register int x,y,i;

    /* Reject Clip Polygon */
    if( UseSlabPlane )
        for( i=0; i<p->count; i++ )
            if( p->v[i].z >= SlabValue )
                return;

    /* Find top vertex */
    top = 0;  
    ymin = p->v[0].y;
    for( i=1; i<p->count; i++ )
       if( p->v[i].y < ymin )
       {   ymin = p->v[i].y;
           top = i;
       }

    rem = p->count;
    ly = ry = y = ymin;
    li = ri = top;

    if( y<0 )
    {   rem--;

        while( ly<=0 && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > 0 )
            {   dy = p->v[i].y - ly;
                lft.di = (((Long)(p->v[i].inten - p->v[li].inten))<<16)/dy;
                lft.dx = (((Long)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((Long)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.i = (((Long)p->v[li].inten)<<16) - (Long)ly*lft.di;
                lft.x = (((Long)p->v[li].x)<<16) - (Long)ly*lft.dx;
                lft.z = (((Long)p->v[li].z)<<16) - (Long)ly*lft.dz;
            } else rem--;
            ly = p->v[i].y;
            li = i;
        }

        while( ry<=0 && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > 0 )
            {   dy = p->v[i].y - ry;
                rgt.di = (((Long)(p->v[i].inten - p->v[ri].inten))<<16)/dy;
                rgt.dx = (((Long)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((Long)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.i = (((Long)p->v[ri].inten)<<16) - (Long)ry*rgt.di;
                rgt.x = (((Long)p->v[ri].x)<<16) - (Long)ry*rgt.dx;
                rgt.z = (((Long)p->v[ri].z)<<16) - (Long)ry*rgt.dz;
            } else rem--;
            ry = p->v[i].y;
            ri = i;
        }

        fbase = View.fbuf;
        dbase = View.dbuf;
        y = 0;
    } else /* y >= 0 */
    {   offset = (Long)y*View.yskip;
        fbase = View.fbuf+offset;
        dbase = View.dbuf+offset;
    }

    while( rem )
    {   while( ly<=y && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ly;
                lft.di = (((Long)(p->v[i].inten - p->v[li].inten))<<16)/dy;
                lft.dx = (((Long)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((Long)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.i = ((Long)p->v[li].inten)<<16;
                lft.x = ((Long)p->v[li].x)<<16;
                lft.z = ((Long)p->v[li].z)<<16;
            }
            ly = p->v[i].y;
            rem--;  li = i;
        }

        while( ry<=y && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ry;
                rgt.di = (((Long)(p->v[i].inten - p->v[ri].inten))<<16)/dy;
                rgt.dx = (((Long)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((Long)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.i = ((Long)p->v[ri].inten)<<16;
                rgt.x = ((Long)p->v[ri].x)<<16;
                rgt.z = ((Long)p->v[ri].z)<<16;
            }
            ry = p->v[i].y;
            rem--; ri = i;
        }


        ymin = MinFun(ly,ry);
        if( ymin>View.ymax )
        {   ymin = View.ymax;
            rem = 0;
        }
        
        while( y<ymin )
        {   if( lft.x < rgt.x )
            {   pmin = &lft;
                pmax = &rgt;
            } else
            {   pmin = &rgt;
                pmax = &lft;
            }

            xmax = (int)(pmax->x>>16)+1;
            xmin = (int)(pmin->x>>16);

            if( (xmin<View.xmax) && (xmax>=0) )
            {   di = (Long)((pmax->i-pmin->i)/(xmax-xmin));
                dz = (Long)((pmax->z-pmin->z)/(xmax-xmin));
                if( xmin<0 )
                {   inten = pmin->i - xmin*di;
                    z = pmin->z - xmin*dz;
                    xmin = 0;
                } else /* xmin >= 0 */
                {   inten = pmin->i;  
                    z = pmin->z;
                }

                if( xmax>=View.xmax )
                    xmax = View.xmax;

                dptr = dbase+xmin;
                for( x=xmin; x<xmax; x++ )
                {   if( (int)(z>>16) > *dptr )
                    {   fbase[x] = Lut[(int)(inten>>16)];
                        *dptr = (int)(z>>16);
                    }
                    inten += di;
                    z += dz;
                    dptr++;
                }
            }

            lft.x += lft.dx;  rgt.x += rgt.dx;
            lft.z += lft.dz;  rgt.z += rgt.dz;
            lft.i += lft.di;  rgt.i += rgt.di;
            dbase += View.yskip;
            fbase += View.yskip;
            y++;
        }
    }
}


/*==========================*/
/*  Render Ribbon Segments  */
/*==========================*/

void SolidRibbon( Knot __far *src, Knot __far *dst, int col )
{
    register int dx, dy;
    Vert p, q, r;

    p.x = src->px-src->wx;  
    p.y = src->py-src->wy;  
    p.z = src->pz-src->wz;
    p.inten = src->vinten+col;

    q.x = dst->px+dst->wx;  
    q.y = dst->py+dst->wy;  
    q.z = dst->pz+dst->wz;
    q.inten = dst->vinten+col;

    dx = q.x - p.x;
    dy = q.y - p.y;

    r.x = dst->px-dst->wx;
    r.y = dst->py-dst->wy;  
    r.z = dst->pz-dst->wz;
    r.inten = dst->vinten+col;

    if( dst->wx*dy > dst->wy*dx )
      ClipTriangle( &p, &q, &r );
   else
      ClipTriangle( &p, &r, &q );

    r.x = src->px+src->wx;  
    r.y = src->py+src->wy;  
    r.z = src->pz+src->wz;
    r.inten = src->vinten+col;

    if( src->wx*dy < src->wy*dx )
      ClipTriangle( &p, &q, &r );
   else
      ClipTriangle( &p, &r, &q );
}


void SolidRibbon2( Knot __far *src, Knot __far *dst, int col1, int col2 )
{
    register int dx,dy;
    Vert p, q, r;

    p.x = src->px-src->wx;  
    p.y = src->py-src->wy;  
    p.z = src->pz-src->wz;

    q.x = dst->px+dst->wx;  
    q.y = dst->py+dst->wy;  
    q.z = dst->pz+dst->wz;

    dx = q.x - p.x;
    dy = q.y - p.y;

    r.x = dst->px-dst->wx;
    r.y = dst->py-dst->wy;  
    r.z = dst->pz-dst->wz;

    if( dst->wx*dy > dst->wy*dx )
    {
      p.inten = src->vinten+col2;
      q.inten = dst->vinten+col2;
      r.inten = dst->vinten+col2;
      ClipTriangle( &p, &q, &r );
    } else {
      p.inten = src->vinten+col1;
      q.inten = dst->vinten+col1;
      r.inten = dst->vinten+col1;
      ClipTriangle( &p, &r, &q );
    }

    r.x = src->px+src->wx;  
    r.y = src->py+src->wy;  
    r.z = src->pz+src->wz;

    if( src->wx*dy < src->wy*dx )
    {
      p.inten = src->vinten+col1;
      q.inten = dst->vinten+col1;
      r.inten = src->vinten+col1;
      ClipTriangle( &p, &q, &r );
    }
    else
    {
      p.inten = src->vinten+col2;
      q.inten = dst->vinten+col2;
      r.inten = src->vinten+col2;
      ClipTriangle( &p, &r, &q );
    }
}


void RectRibbon( Knot __far *src, Knot __far *dst, int col,
                 int scap, int dcap )
{
    Vert s00, s01, s10, s11;
    Vert d00, d01, d10, d11;
    int svcol, shcol;
    int dvcol, dhcol;
    int size, inten;
    Real rinten;

    if (src) {
        s00.x = src->px-src->wx-src->dx;
        s00.y = src->py-src->wy-src->dy;
        s00.z = src->pz-src->wz-src->dz;

        s01.x = src->px-src->wx+src->dx;
        s01.y = src->py-src->wy+src->dy;
        s01.z = src->pz-src->wz+src->dz;

        s10.x = src->px+src->wx-src->dx;
        s10.y = src->py+src->wy-src->dy;
        s10.z = src->pz+src->wz-src->dz;

        s11.x = src->px+src->wx+src->dx;
        s11.y = src->py+src->wy+src->dy;
        s11.z = src->pz+src->wz+src->dz;

        shcol = src->hinten+col;
        svcol = src->vinten+col;
    } else shcol = svcol = 0;

    if (dst) {
        d00.x = dst->px-dst->wx-dst->dx;
        d00.y = dst->py-dst->wy-dst->dy;
        d00.z = dst->pz-dst->wz-dst->dz;
    
        d01.x = dst->px-dst->wx+dst->dx;
        d01.y = dst->py-dst->wy+dst->dy;
        d01.z = dst->pz-dst->wz+dst->dz;
    
        d10.x = dst->px+dst->wx-dst->dx;
        d10.y = dst->py+dst->wy-dst->dy;
        d10.z = dst->pz+dst->wz-dst->dz;

        d11.x = dst->px+dst->wx+dst->dx;
        d11.y = dst->py+dst->wy+dst->dy;
        d11.z = dst->pz+dst->wz+dst->dz;
    
        dhcol = dst->hinten+col;
        dvcol = dst->vinten+col;
    } else dhcol = dvcol = 0;
    
    if (src && dst)
    {
        s00.inten = shcol;  d00.inten = dhcol;
        s01.inten = shcol;  d01.inten = dhcol;
        ClipTriangle( &s00, &d00, &s01 );
        ClipTriangle( &d00, &d01, &s01 );

        s01.inten = svcol;  d01.inten = dvcol;
        s11.inten = svcol;  d11.inten = dvcol;
        ClipTriangle( &s01, &d01, &s11 );
        ClipTriangle( &d01, &d11, &s11 );

        s11.inten = shcol;  d11.inten = dhcol;
        s10.inten = shcol;  d10.inten = dhcol;
        ClipTriangle( &s11, &d11, &s10 );
        ClipTriangle( &d11, &d10, &s10 );

        s10.inten = svcol;  d10.inten = dvcol;
        s00.inten = svcol;  d00.inten = dvcol;
        ClipTriangle( &s10, &d10, &s00 );
        ClipTriangle( &d10, &d00, &s00 );
    }

    if (scap)
    {
      size = isqrt( (Long)src->tx*src->tx +
                    (Long)src->ty*src->ty +
                    (Long)src->tz*src->tz ) + 1;
#ifdef ORIGINAL
      rinten = LightDot(-src->tx,InvertY(src->ty),-src->tz);
      rinten /= size*LightLength;
#else
      rinten = (Real)-src->tz/size;
#endif

      if( src->tz > 0 ) rinten = -rinten;

      if( rinten > 0.0 )
      {   inten = (char)(ColourMask*rinten) + col;
      } else inten = col;

      s00.inten = inten;  s01.inten = inten;
      s10.inten = inten;  s11.inten = inten;
      ClipTriangle( &s00, &s01, &s10 );
      ClipTriangle( &s10, &s01, &s11 );
    }

    if (dcap)
    {
      size = isqrt( (Long)dst->tx*dst->tx +
                    (Long)dst->ty*dst->ty +
                    (Long)dst->tz*dst->tz ) + 1;
#ifdef ORIGINAL
      rinten = LightDot(dst->tx,-InvertY(dst->ty),dst->tz);
      rinten /= size*LightLength;
#else
      rinten = (Real)dst->tz/size;
#endif

      if( dst->tz < 0 ) rinten = -rinten;

      if( rinten > 0.0 )
      {   inten = (char)(ColourMask*rinten) + col;
      } else inten = col;

      d00.inten = inten;  d01.inten = inten;
      d10.inten = inten;  d11.inten = inten;
      ClipTriangle( &d00, &d10, &d01 );
      ClipTriangle( &d01, &d10, &d11 );
    }
}


static int TestSphere( int x, int y, int z, int rad )
{
    register int temp;

    ClipStatus = 0;

    if( UseSlabPlane )
    {   if( z-rad>=SlabValue )
            return( False );

        if( z+rad>=SlabValue )
        {   if( SlabMode )
            {   ClipStatus |= BitFront;
            } else return( False );
        } else if( SlabMode==SlabSection )
            return( False );
    }

    temp = x+rad;
    if( temp<0 ) return( False );
    if( temp>=View.xmax ) ClipStatus |= BitRight;

    temp = x-rad;
    if( temp>=View.xmax ) return( False );
    if( temp<0 ) ClipStatus |= BitLeft;

    temp = y+rad;
    if( temp<0 ) return( False );
    if( temp>=View.ymax ) ClipStatus |= BitBelow;

    temp = y-rad;
    if( temp>=View.ymax ) return( False );
    if( temp<0 ) ClipStatus |= BitAbove;

    return True;
}



/*===========================*/
/*  Sphere Rendering Macros  */
/*===========================*/

#define UpdateAcross(dz)    \
        depth = (dz)+z;                    \
        if( depth > *dptr )                \
        {   *dptr = depth;                 \
            fptr = fold+dx;                \
            inten = LightDot(dx,dy,dz);    \
            if( inten>0 )                  \
            {      inten = (int)((inten*ColConst[rad])>>ColBits); \
                   *fptr = Lut[col+inten]; \
            } else *fptr = Lut[col];       \
        }                                  \
        dptr++;  dx++;

#define UpdateLine  \
        dx = -wide;                   \
        dptr = dold-wide;             \
        tptr = LookUp[wide]+wide;     \
        while( dx<0 ) { UpdateAcross(*tptr); tptr--; }       \
        do { UpdateAcross(*tptr); tptr++; } while(dx<=wide); \
        dold += View.yskip;  fold += View.yskip;             \
        dy++;


void DrawSphere( int x, int y, int z, int rad, int col )
{
    register Pixel __huge *fptr, __huge *fold;
    register short __huge *dptr, __huge *dold;
    register unsigned char __far *tptr;

    register Long offset;
    register int depth,wide,inten;
    register int dx,dy;

    /* Avoid Lookup Table Overflow! */
    if( rad > MAXRAD ) rad = MAXRAD;

    offset = (Long)(y-rad)*View.yskip + x;
    fold=View.fbuf+offset;  
    dold=View.dbuf+offset;

    dy = -rad;
    while( dy<0 ) 
    {   wide = LookUp[rad][-dy]; 
        UpdateLine; 
    }

    do { 
        wide = LookUp[rad][dy];  
        UpdateLine; 
    } while( dy<=rad );
}


void ClipSphere( int x, int y, int z, int rad, int col )
{
    register Pixel __huge *fptr, __huge *fold;
    register short __huge *dptr, __huge *dold;

    register int lastx,lasty,dx,dy,dz;
    register int depth,wide,inten,side;
    register int crad,cwide,temp;
    register Long offset;

    /* Avoid Lookup Table Overflow! */
    if( rad > MAXRAD ) rad = MAXRAD;

    /* Visibility Tests */
    if( !TestSphere(x,y,z,rad) )
        return;

    if( !ClipStatus )
    {   DrawSphere(x,y,z,rad,col);
        return;
    }

    if( ClipStatus&BitAbove )
    {   dy = -y;
        fold = View.fbuf + x;
        dold = View.dbuf + x;
    } else
    {   dy = -rad;
        offset = (Long)(y+dy)*View.yskip+x;
        fold = View.fbuf + offset;
        dold = View.dbuf + offset;
    }

    if( ClipStatus&BitBelow )
    {   lasty = (View.ymax-1)-y;
    } else lasty = rad;

    side = (View.xmax-1)-x;
    /* No Slab Plane Clipping */
    if( !(ClipStatus&BitFront) )
    {   while( dy<=lasty )
        {   wide = LookUp[rad][AbsFun(dy)];
            lastx = MinFun(wide,side);
            dx = - MinFun(wide,x);
            dptr = dold + dx;

            while( dx<=lastx )
            {   dz = LookUp[wide][AbsFun(dx)];
                UpdateAcross(dz);
            }
            dold += View.yskip;
            fold += View.yskip;
            dy++;
        }
        return;
    }

    dz = SlabValue-z;
    crad = LookUp[rad][AbsFun(dz)];

    if( (z>SlabValue) || (SlabMode==SlabSection) )
    {   if( crad<lasty ) lasty = crad;
        if( -crad>dy ) 
        {   dy = -crad;
            offset = (Long)(y+dy)*View.yskip+x;
            fold = View.fbuf + offset;
            dold = View.dbuf + offset;
        }
    }

    while( dy<=lasty )
    {   temp = AbsFun(dy);
        wide = LookUp[rad][temp];
        lastx = MinFun(wide,side);
        dx = - MinFun(x,wide);
        dptr = dold + dx;

        if( temp<=crad )
        {   cwide = LookUp[crad][temp];
            while( dx<=lastx )
            {   temp = AbsFun(dx);
                if( temp<=cwide )
                {    /* Slab Plane Clipping Modes */
                    switch( SlabMode )
                    {   case( SlabFinal ):
                                fold[dx] = Lut[col+SlabInten];
                                *dptr = SliceValue;
                                break;

                        case( SlabHollow ):
                                dz = LookUp[wide][temp];
                                depth = z - dz;
                                if( depth>*dptr )
                                {   *dptr = depth;
                                    inten = LightDot(-dx,-dy,dz);

                                    if( inten>0 )
                                    {   inten=(int)( (inten*ColConst[rad])
                                                     >>(ColBits+1));
                                        fold[dx] = Lut[col+inten];
                                    } else fold[dx] = Lut[col];
                                }
                                break;

                        case( SlabSection ):
                        case( SlabClose ):
                                dz = SlabValue-z;
                                depth = dx*dx+dy*dy+dz*dz+SliceValue;
                                if( (*dptr<SliceValue) || (depth<*dptr) )
                                {   fold[dx] = Lut[col+SlabInten];
                                    *dptr = depth;
                                }
                                break;
                    }
                    dptr++;  dx++;
                } else if( (z<SlabValue) && (SlabMode!=SlabSection) )
                {    dz = LookUp[wide][temp];
                     UpdateAcross(dz);
                } else
                {   dptr++;  dx++;
                }
            }
        } else /* Slabless ScanLine */
            while( dx<=lastx )
            {   dz = LookUp[wide][AbsFun(dx)];
                UpdateAcross(dz);
            }

        dold += View.yskip;
        fold += View.yskip;
        dy++;
    }
}


static void DrawArcAc( short __huge *dbase,
                       Pixel __huge *fbase,
                       int z, int c )
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;

    for( ptr=ArcAc; ptr<ArcAcPtr; ptr++ )
    {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
        SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
    }
}


static void DrawArcDn( short __huge *dbase,
                       Pixel __huge *fbase,
                       int z, int c )
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;

    for( ptr=ArcDn; ptr<ArcDnPtr; ptr++ )
    {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
        SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
    }
}


static void DrawCylinderCaps( int x1, int y1, int z1,
                              int x2, int y2, int z2,
                              int c1, int c2, int rad )
{
    register short __huge *dold, __huge *dptr;
    register Pixel __huge *fold;
#ifdef UNUSED
    register int ax,ay,ix,iy;
    register int zrate,lz;
#endif
    register Long offset,temp,end;
    register int inten,absx;
    register int wide,depth;
    register int dx,dy,dz;
    register int lx,ly;

    lx = x2-x1;
    ly = y2-y1;

#ifdef UNUSED
    lz = z2-z1;
    if( ly>0 ) { ay = ly; iy = 1; }
    else { ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);
#endif

    end = (Long)ly*View.yskip+lx;
    temp = (Long)y1*View.yskip+x1;
    fold = View.fbuf+temp;
    dold = View.dbuf+temp;

    ArcAcPtr = ArcAc;
    ArcDnPtr = ArcDn;

    temp = (Long)-(rad*View.yskip);
    for( dy= -rad; dy<=rad; dy++ )
    {   wide = LookUp[rad][AbsFun(dy)];

        for( dx= -wide; dx<=wide; dx++ )
        {   absx = AbsFun(dx);
            dz = LookUp[wide][absx];
            inten = LightDot(dx,dy,dz);
            if( inten>0 )
            {   inten = (int)((inten*ColConst[rad])>>ColBits);
            } else inten = 0;
            offset = temp+dx;

            if( XValid(x1+dx) && YValid(y1+dy) )
            {   dptr = dold+offset; depth = z1+dz;
                SETPIXEL(dptr,fold+offset,depth,Lut[c1+inten]);
            }

            if( XValid(x2+dx) && YValid(y2+dy) )
            {   dptr = dold+(offset+end); depth = z2+dz;
                SETPIXEL(dptr,fold+(offset+end),depth,Lut[c2+inten]);
            }

#ifdef UNUSED
            k1 = AbsFun(dx+ix); 
            k2 = AbsFun(dx-ix);

            if( ((k1>wide)||(dz>=LookUp[wide][k1]-zrate)) &&
                ((k2>wide)||(dz>LookUp[wide][k2]+zrate)) )
#else
            if( ArcAcPtr < ArcAc+ARCSIZE )
#endif
            {   ArcAcPtr->offset = offset; ArcAcPtr->inten = inten;
                ArcAcPtr->dx=dx; ArcAcPtr->dy=dy; ArcAcPtr->dz=dz;
                ArcAcPtr++;
            }

#ifdef UNUSED
            k1 = AbsFun(dy+iy);
            k2 = AbsFun(dy-iy);

            high = LookUp[rad][absx];
            if( ((k1>high)||(dz>=LookUp[LookUp[rad][k1]][absx]-zrate)) &&
                ((k2>high)||(dz>LookUp[LookUp[rad][k2]][absx]+zrate)) )
#else
            if( ArcDnPtr < ArcDn+ARCSIZE )
#endif
            {   ArcDnPtr->offset = offset; ArcDnPtr->inten = inten;
                ArcDnPtr->dx=dx; ArcDnPtr->dy=dy; ArcDnPtr->dz=dz;
                ArcDnPtr++;
            }
        }
        temp += View.yskip;
    }
}


void DrawCylinder( int x1, int y1, int z1,
                   int x2, int y2, int z2,
                   int c1, int c2, int rad )
{
    register short __huge *dbase;
    register Pixel __huge *fbase;

    register int zrate,zerr,ystep,err;
    register int ix,iy,ax,ay;
    register int lx,ly,lz;
    register int mid,tmp;
    register Long temp;

    /* Avoid Lookup Table Overflow! */
    if( rad > MAXRAD ) rad = MAXRAD;

    /* Trivial Case */
    if( (x1==x2) && (y1==y2) )
    {   if( z1>z2 )
        {      DrawSphere(x1,y1,z1,rad,c1);
        } else DrawSphere(x2,y2,z2,rad,c2);
        return;
    }

    if( z1<z2 )
    {   tmp=x1; x1=x2; x2=tmp;
        tmp=y1; y1=y2; y2=tmp;
        tmp=z1; z1=z2; z2=tmp;
        tmp=c1; c1=c2; c2=tmp;
    }

    DrawCylinderCaps(x1,y1,z1,x2,y2,z2,c1,c2,rad);

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ystep = View.yskip; ay = ly; iy = 1; }
    else {   ystep = -View.yskip; ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);

    temp = (Long)y1*View.yskip+x1;
    fbase = View.fbuf+temp;
    dbase = View.dbuf+temp;

    if( ax>ay )
    {   lz -= ax*zrate;
        zerr = err = -(ax>>1);

        if( c1 != c2 )
        {   mid = (x1+x2)>>1;
            while( x1!=mid )
            {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
                fbase+=ix; dbase+=ix; x1+=ix;
                if( (err+=ay)>0 )
                {   fbase+=ystep; dbase+=ystep; err-=ax;
                       DrawArcDn(dbase,fbase,z1,c1);
                } else DrawArcAc(dbase,fbase,z1,c1);
            }
        }

        while( x1!=x2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
            fbase+=ix; dbase+=ix; x1+=ix;
            if( (err+=ay)>0 )
            {   fbase+=ystep; dbase+=ystep; err-=ax;
                   DrawArcDn(dbase,fbase,z1,c2);
            } else DrawArcAc(dbase,fbase,z1,c2);
        }
    } else /*ay>=ax*/
    {   lz -= ay*zrate;
        zerr = err = -(ay>>1);

        if( c1 != c2 )
        {   mid = (y1+y2)>>1;
            while( y1!=mid )
            {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
                fbase+=ystep; dbase+=ystep; y1+=iy;
                if( (err+=ax)>0 )
                {   fbase+=ix; dbase+=ix; err-=ay; 
                       DrawArcAc(dbase,fbase,z1,c1);
                } else DrawArcDn(dbase,fbase,z1,c1);
            }
        }

        while( y1!=y2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
            fbase+=ystep; dbase+=ystep; y1+=iy;
            if( (err+=ax)>0 )
            {   fbase+=ix; dbase+=ix; err-=ay; 
                   DrawArcAc(dbase,fbase,z1,c2);
            } else DrawArcDn(dbase,fbase,z1,c2);
        }
    }
}


static int TestCylinder( int x1, int y1, int z1,
                         int x2, int y2, int z2,
                         int rad )
{
    register int tmp1, tmp2;

    if( UseSlabPlane )
        if( (z1+rad>SlabValue) || (z2+rad>SlabValue) )
            return(False);

    ClipStatus = False;

    tmp1 = x1+rad;  tmp2 = x2+rad;
    if( (tmp1<0) && (tmp2<0) )
        return( False );
    if( (tmp1>=View.xmax) || (tmp2>=View.xmax) )
        ClipStatus = True;

    tmp1 = x1-rad;  tmp2 = x2-rad;
    if( (tmp1>=View.xmax) && (tmp2>=View.xmax) )
        return( False );
    if( (tmp1<0) || (tmp2<0) )
        ClipStatus = True;

    tmp1 = y1+rad;  tmp2 = y2+rad;
    if( (tmp1<0) && (tmp2<0) )
        return( False );
    if( (tmp1>=View.ymax) || (tmp2>=View.ymax) )
        ClipStatus = True;

    tmp1 = y1-rad;  tmp2 = y2-rad;
    if( (tmp1>=View.ymax) && (tmp2>=View.ymax) )
        return( False );
    if( (tmp1<0) || (tmp2<0) )
        ClipStatus = True;

    return True;
}


static void ClipArcAc( short __huge *dbase, 
                       Pixel __huge *fbase,
                       int x, int y, int z, int c )
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;
    register int temp;

    ptr = ArcAc;
    while( (temp=y+ptr->dy) < 0 )
        if( ++ptr == ArcAcPtr )
            return;

    while( (temp<View.ymax) && (ptr<ArcAcPtr) )
    {   temp = x+ptr->dx;
        if( XValid(temp) )
        {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
            SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
        }
        ptr++;
        temp = y+ptr->dy;
    }
}


static void ClipArcDn( short __huge *dbase,
                       Pixel __huge *fbase,
                       int x, int y, int z, int c )
{
    register ArcEntry __far *ptr;
    register short __huge *dptr;
    register short depth;
    register int temp;

    ptr = ArcDn;
    while( (temp=y+ptr->dy) < 0 )
        if( ++ptr == ArcDnPtr )
            return;

    while( (temp<View.ymax) && (ptr<ArcDnPtr) )
    {   temp = x+ptr->dx;
        if( XValid(temp) )
        {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
            SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
        }
        ptr++;
        temp = y+ptr->dy;
    }
}


void ClipCylinder( int x1, int y1, int z1,
                   int x2, int y2, int z2,
                   int c1, int c2, int rad )
{
    register short __huge *dbase;
    register Pixel __huge *fbase;

    register int zrate,zerr,ystep,err;
    register int ix,iy,ax,ay;
    register int lx,ly,lz;
    register int mid,tmp;
    register Long temp;

    /* Avoid Lookup Table Overflow! */
    if( rad > MAXRAD ) rad = MAXRAD;

    /* Visibility Tests */
    if( !TestCylinder(x1,y1,z1,x2,y2,z2,rad) )
        return;

    if( !ClipStatus )
    {   DrawCylinder(x1,y1,z1,x2,y2,z2,c1,c2,rad);
        return;
    }

    /* Trivial Case */
    if( (x1==x2) && (y1==y2) )
    {   if( z1>z2 )
        {      ClipSphere(x1,y1,z1,rad,c1);
        } else ClipSphere(x2,y2,z2,rad,c2);
        return;
    }

    if( z1<z2 )
    {   tmp=x1; x1=x2; x2=tmp;
        tmp=y1; y1=y2; y2=tmp;
        tmp=z1; z1=z2; z2=tmp;
        tmp=c1; c1=c2; c2=tmp;
    }

    DrawCylinderCaps(x1,y1,z1,x2,y2,z2,c1,c2,rad);

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ystep = View.yskip; ay = ly; iy = 1; }
    else {   ystep = -View.yskip; ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);

    temp = (Long)y1*View.yskip+x1;
    fbase = View.fbuf+temp;
    dbase = View.dbuf+temp;

    if( ax>ay )
    {   if( x2<x1 )
        {   tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        lz -= ax*zrate;
        zerr = err = -(ax>>1);
        mid = (x1+x2)/2;

        while( x1!=x2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
            fbase+=ix; dbase+=ix; x1+=ix;
            if( (err+=ay)>0 )
            {   fbase += ystep;  err -= ax;
                dbase += ystep;  y1 += iy;
                   ClipArcDn(dbase,fbase,x1,y1,z1,(x1<mid?c1:c2));
            } else ClipArcAc(dbase,fbase,x1,y1,z1,(x1<mid?c1:c2));
        }
    } else /*ay>=ax*/
    {   if( y2<y1 )
        {   tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        lz -= ay*zrate;
        zerr = err = -(ay>>1);
        mid = (y1+y2)/2;

        while( y1!=y2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
            fbase+=ystep; dbase+=ystep; y1+=iy;
            if( (err+=ax)>0 )
            {   fbase += ix;  err -= ay;
                dbase += ix;  x1 += ix; 
                   ClipArcAc(dbase,fbase,x1,y1,z1,(y1<mid?c1:c2));
            } else ClipArcDn(dbase,fbase,x1,y1,z1,(y1<mid?c1:c2));
        }
    }
}



/*================================*/
/*  Character Rendering Routines  */
/*================================*/

void SetFontSize( int size )
{
    register int count;
    register int i;

    count = 0;
    for( i=0; i<23; i++ )
    {   FontDimen[i] = count>>4;
        count += size;
    }
    FontSize = size;
}


static void ClipCharacter( int x, int y, int z, int glyph, int col )
{
    register char *ptr;
    register int sx,sy;
    register int ex,ey;

    /* Avoid compiler warnings! */
    ex = 0;  ey = 0;

    ptr = VectFont[glyph];
    while( *ptr )
    {   /* Uppercase test */
        if( ptr[0] < 'a' )
        {   sx = x + FontDimen[ptr[0]-'A'];
            sy = y + InvertY(FontDimen[ptr[1]-'a']);
            ptr += 2;
        } else
        {   sx = ex;
            sy = ey;
        }

        ex = x + FontDimen[ptr[0]-'a'];
        ey = y + InvertY(FontDimen[ptr[1]-'a']);
        if( (ex!=sx) || (ey!=sy) )
        {   ClipLine(sx,sy,z,ex,ey,z,col);
        } else ClipPoint(ex,ey,z,col);
        ptr += 2;
    }
}


void DisplayString( int x, int y, int z, char *label, int col )
{
    register int clip,high,max;
    register char *ptr;
    register int sx,sy;
    register int ex,ey;

    /* Avoid compiler warnings! */
    ex = 0;  ey = 0;

    high = (FontSize*3)>>1;
#ifdef INVERT
    if( ((y+high)<0) || (y>=View.ymax) ) return;
    clip = (y<0) || (y+high>=View.ymax);
#else
    if( (y<0) || ((y-high)>=View.ymax) ) return;
    clip = (y-high<0) || (y>=View.ymax);
#endif

    if( x < 0 )
    {   while( *label && (x<=-FontSize) )
        {   x += FontSize;
            label++;
        }

        if( *label )
        {   ClipCharacter(x,y,z,(*label-32),col);
            x += FontSize;
            label++;
        } else return;
    }

    if( !clip )
    {   max = View.xmax-FontSize;
        while( *label && (x<max) )
        {  ptr = VectFont[*label-32];
           while( *ptr )
           {   /* Uppercase test */
               if( ptr[0] < 'a' )
               {   sx = x + FontDimen[ptr[0]-'A'];
                   sy = y + InvertY(FontDimen[ptr[1]-'a']);
                   ptr += 2;
               } else
               {   sx = ex;
                   sy = ey;
               }

               ex = x + FontDimen[ptr[0]-'a'];
               ey = y + InvertY(FontDimen[ptr[1]-'a']);
               if( (ex!=sx) || (ey!=sy) )
               {   DrawTwinLine(sx,sy,z,ex,ey,z,col,col);
               } else PlotPoint(ex,ey,z,col);
               ptr += 2;
           }
           x += FontSize;
           label++;
        }

        if( *label )
            ClipCharacter(x,y,z,(*label-32),col);
    } else /* Always Clip! */
        while( *label && (x<View.xmax) )
        {   ClipCharacter(x,y,z,(*label-32),col);
            x += FontSize;
            label++;
        }
}


void InitialisePixUtils( void )
{
#if defined(IBMPC) || defined(APPLEMAC)
    ArcAc = (ArcEntry __far*)_fmalloc(ARCSIZE*sizeof(ArcEntry));
    ArcDn = (ArcEntry __far*)_fmalloc(ARCSIZE*sizeof(ArcEntry));
#endif
    SplineCount = 5;
    SetFontSize(8);
}

