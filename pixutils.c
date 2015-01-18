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
#define ARCSIZE  2048

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

#define OutCode(res,x,y,z)            \
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
                                      \
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


#ifdef UNUSED
static int OutCode( int x, int y, int z )
{
    register int result;

    if( y < 0 )
    {   result = BitAbove;
    } else if( y >= View.ymax )
    {   result = BitBelow;
    } else result = 0;

    if( x < 0 )
    {   result |= BitLeft;
    } else if( x >= View.xmax )
        result |= BitRight;

    if( !ZValid(z) )
        result |= BitFront;
    return result;
}
#endif


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
    {  fptr = View.fbuf+offset;
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

    OutCode(code1,x1,y1,z1);
    OutCode(code2,x2,y2,z2);
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

        OutCode(code1,x1,y1,z1);
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
    {   OutCode(code1,x1,y1,z1);
        OutCode(code2,x2,y2,z2);
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

    OutCode(code1,x1,y1,z1);
    OutCode(code2,x2,y2,z2);
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
        OutCode(code1,x1,y1,z1);
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
    {   OutCode(code1,x1,y1,z1);
        OutCode(code2,x2,y2,z2);
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

    /* Reject(OutCode(x1,y1,z1),OutCode(x2,y2,z2)) */
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


#ifdef UNUSED  /* Unused Function */
static void OutLinePolygon( Poly *p )
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
#endif


#ifdef UNUSED
static void DrawPolygon( Poly *p )
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
#endif


#ifdef UNUSED
static void DrawFlatPolygon( Poly *p )
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
#endif


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


#ifdef UNUSED
static void ClipFlatPolygon( Poly *p )
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
#endif


void SolidRibbon( Knot __far *src, Knot __far *dst, int col )
{
    static Poly p;

    p.v[0].x = src->px+src->wx;  
    p.v[0].y = src->py+src->wy;  
    p.v[0].z = src->pz+src->wz;
    p.v[0].inten = src->vinten+col;

    p.v[1].x = dst->px+dst->wx;  
    p.v[1].y = dst->py+dst->wy;  
    p.v[1].z = dst->pz+dst->wz;
    p.v[1].inten = dst->vinten+col;

    p.v[2].x = dst->px-dst->wx;
    p.v[2].y = dst->py-dst->wy;  
    p.v[2].z = dst->pz-dst->wz;
    p.v[2].inten = dst->vinten+col;

    p.v[3].x = src->px-src->wx;  
    p.v[3].y = src->py-src->wy;  
    p.v[3].z = src->pz-src->wz;
    p.v[3].inten = src->vinten+col;

    p.count = 4;
    /* OutLinePolygon( &p ); */
    ClipPolygon( &p );
}


void SolidRibbon2( Knot __far *src, Knot __far *dst, int col1, int col2 )
{
    register int dx,dy;
    register int col;
    static Poly p;

    p.count = 3;
    p.v[0].x = src->px+src->wx;  
    p.v[0].y = src->py+src->wy;  
    p.v[0].z = src->pz+src->wz;
    p.v[1].x = dst->px-dst->wx;  
    p.v[1].y = dst->py-dst->wy;  
    p.v[1].z = dst->pz-dst->wz;

    dx = p.v[1].x - p.v[0].x;
    dy = p.v[1].y - p.v[0].y;

    p.v[2].x = dst->px+dst->wx;
    p.v[2].y = dst->py+dst->wy;  
    p.v[2].z = dst->pz+dst->wz;

#ifdef INVERT
    col = ( dst->wx*dy > dst->wy*dx )? col2 : col1;
#else
    col = ( dst->wx*dy < dst->wy*dx )? col2 : col1;
#endif
    p.v[0].inten = src->vinten+col;
    p.v[1].inten = dst->vinten+col;
    p.v[2].inten = dst->vinten+col;

    /* OutLinePolygon( &p ); */
    ClipPolygon( &p );

    p.v[2].x = src->px-src->wx;  
    p.v[2].y = src->py-src->wy;  
    p.v[2].z = src->pz-src->wz;

#ifdef INVERT
    col = ( src->wx*dy > src->wy*dx )? col2 : col1;
#else
    col = ( src->wx*dy < src->wy*dx )? col2 : col1;
#endif

    p.v[0].inten = src->vinten+col;
    p.v[1].inten = dst->vinten+col;
    p.v[2].inten = src->vinten+col;
    /* OutLinePolygon( &p ); */
    ClipPolygon( &p );
}


void RectRibbon( Knot __far *src, Knot __far *dst, int col )
{
    static Poly p;

    p.count = 4;

    p.v[0].inten = src->vinten+col;
    p.v[1].inten = dst->vinten+col;
    p.v[2].inten = dst->vinten+col;
    p.v[3].inten = src->vinten+col;

    /* Top Surface */
    p.v[0].x = src->px+src->wx+src->dx;  
    p.v[0].y = src->py+src->wy+src->dy;  
    p.v[0].z = src->pz+src->wz+src->dz;

    p.v[1].x = dst->px+dst->wx+dst->dx;  
    p.v[1].y = dst->py+dst->wy+dst->dy;  
    p.v[1].z = dst->pz+dst->wz+dst->dz;

    p.v[2].x = dst->px-dst->wx+dst->dx;
    p.v[2].y = dst->py-dst->wy+dst->dy;  
    p.v[2].z = dst->pz-dst->wz+dst->dz;

    p.v[3].x = src->px-src->wx+src->dx;  
    p.v[3].y = src->py-src->wy+src->dy;  
    p.v[3].z = src->pz-src->wz+src->dz;
    ClipPolygon( &p );

    /* Bottom Surface */
    p.v[0].x = src->px+src->wx-src->dx;  
    p.v[0].y = src->py+src->wy-src->dy;  
    p.v[0].z = src->pz+src->wz-src->dz;

    p.v[1].x = dst->px+dst->wx-dst->dx;  
    p.v[1].y = dst->py+dst->wy-dst->dy;  
    p.v[1].z = dst->pz+dst->wz-dst->dz;

    p.v[2].x = dst->px-dst->wx-dst->dx;
    p.v[2].y = dst->py-dst->wy-dst->dy;  
    p.v[2].z = dst->pz-dst->wz-dst->dz;

    p.v[3].x = src->px-src->wx-src->dx;  
    p.v[3].y = src->py-src->wy-src->dy;  
    p.v[3].z = src->pz-src->wz-src->dz;
    ClipPolygon( &p );

    p.v[0].inten = src->hinten+col;
    p.v[1].inten = dst->hinten+col;
    p.v[2].inten = dst->hinten+col;
    p.v[3].inten = src->hinten+col;

    /* Left Surface */
    p.v[0].x = src->px+src->wx+src->dx;  
    p.v[0].y = src->py+src->wy+src->dy;  
    p.v[0].z = src->pz+src->wz+src->dz;

    p.v[1].x = dst->px+dst->wx+dst->dx;  
    p.v[1].y = dst->py+dst->wy+dst->dy;  
    p.v[1].z = dst->pz+dst->wz+dst->dz;

    p.v[2].x = dst->px+dst->wx-dst->dx;
    p.v[2].y = dst->py+dst->wy-dst->dy;  
    p.v[2].z = dst->pz+dst->wz-dst->dz;

    p.v[3].x = src->px+src->wx-src->dx;  
    p.v[3].y = src->py+src->wy-src->dy;  
    p.v[3].z = src->pz+src->wz-src->dz;
    ClipPolygon( &p );

    /* Right Surface */
    p.v[0].x = src->px-src->wx+src->dx;  
    p.v[0].y = src->py-src->wy+src->dy;  
    p.v[0].z = src->pz-src->wz+src->dz;

    p.v[1].x = dst->px-dst->wx+dst->dx;  
    p.v[1].y = dst->py-dst->wy+dst->dy;  
    p.v[1].z = dst->pz-dst->wz+dst->dz;

    p.v[2].x = dst->px-dst->wx-dst->dx;
    p.v[2].y = dst->py-dst->wy-dst->dy;  
    p.v[2].z = dst->pz-dst->wz-dst->dz;

    p.v[3].x = src->px-src->wx-src->dx;  
    p.v[3].y = src->py-src->wy-src->dy;  
    p.v[3].z = src->pz-src->wz-src->dz;
    ClipPolygon( &p );
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
    register Byte __far *tptr;

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

