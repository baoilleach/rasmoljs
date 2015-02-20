/* transfor.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
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
#include <stdio.h>
#include <math.h>

#define TRANSFORM
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "cmndline.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "graphics.h"


typedef struct {
                short col;
                short shade;
                unsigned char r;
                unsigned char g;
                unsigned char b;
              } ShadeRef;


#define CPKMAX  16
static ShadeRef CPKShade[] = {
     { 0, 0, 200, 200, 200 },       /*  0 Light Grey   */
     { 0, 0, 143, 143, 255 },       /*  1 Sky Blue     */
     { 0, 0, 240,   0,   0 },       /*  2 Red          */
     { 0, 0, 255, 200,  50 },       /*  3 Yellow       */
     { 0, 0, 255, 255, 255 },       /*  4 White        */
     { 0, 0, 255, 192, 203 },       /*  5 Pink         */
     { 0, 0, 218, 165,  32 },       /*  6 Golden Rod   */
     { 0, 0,   0,   0, 255 },       /*  7 Blue         */
     { 0, 0, 255, 165,   0 },       /*  8 Orange       */
     { 0, 0, 128, 128, 144 },       /*  9 Dark Grey    */
     { 0, 0, 165,  42,  42 },       /* 10 Brown        */
     { 0, 0, 160,  32, 240 },       /* 11 Purple       */
     { 0, 0, 255,  20, 147 },       /* 12 Deep Pink    */
     { 0, 0,   0, 255,   0 },       /* 13 Green        */
     { 0, 0, 178,  34,  34 },       /* 14 Fire Brick   */
     { 0, 0,  34, 139,  34 } };     /* 15 Forest Green */


static ShadeRef Shapely[] = {
     { 0, 0, 140, 255, 140 },    /* ALA */
     { 0, 0, 255, 255, 255 },    /* GLY */
     { 0, 0,  69,  94,  69 },    /* LEU */
     { 0, 0, 255, 112,  66 },    /* SER */
     { 0, 0, 255, 140, 255 },    /* VAL */
     { 0, 0, 184,  76,   0 },    /* THR */
     { 0, 0,  71,  71, 184 },    /* LYS */
     { 0, 0, 160,   0,  66 },    /* ASP */
     { 0, 0,   0,  76,   0 },    /* ILE */
     { 0, 0, 255, 124, 112 },    /* ASN */
     { 0, 0, 102,   0,   0 },    /* GLU */
     { 0, 0,  82,  82,  82 },    /* PRO */
     { 0, 0,   0,   0, 124 },    /* ARG */
     { 0, 0,  83,  76,  66 },    /* PHE */
     { 0, 0, 255,  76,  76 },    /* GLN */
     { 0, 0, 140, 112,  76 },    /* TYR */
     { 0, 0, 112, 112, 255 },    /* HIS */
     { 0, 0, 255, 255, 112 },    /* CYS */
     { 0, 0, 184, 160,  66 },    /* MET */
     { 0, 0,  79,  70,   0 },    /* TRP */

     { 0, 0, 255,   0, 255 },    /* ASX */
     { 0, 0, 255,   0, 255 },    /* GLX */
     { 0, 0, 255,   0, 255 },    /* PCA */
     { 0, 0, 255,   0, 255 },    /* HYP */

     { 0, 0, 160, 160, 255 },    /*   A */
     { 0, 0, 255, 140,  75 },    /*   C */
     { 0, 0, 255, 112, 112 },    /*   G */
     { 0, 0, 160, 255, 160 },    /*   T */

     { 0, 0, 184, 184, 184 },    /* 28 -> BackBone */
     { 0, 0,  94,   0,  94 },    /* 29 -> Special  */
     { 0, 0, 255,   0, 255 } };  /* 30 -> Default  */

     
static ShadeRef AminoShade[] = {
     { 0, 0, 230,  10,  10 },    /*  0  ASP, GLU      */
     { 0, 0,  20,  90, 255 },    /*  1  LYS, ARG      */
     { 0, 0, 130, 130, 210 },    /*  2  HIS           */
     { 0, 0, 250, 150,   0 },    /*  3  SER, THR      */
     { 0, 0,   0, 220, 220 },    /*  4  ASN, GLN      */
     { 0, 0, 230, 230,   0 },    /*  5  CYS, MET      */
     { 0, 0, 200, 200, 200 },    /*  6  ALA           */
     { 0, 0, 235, 235, 235 },    /*  7  GLY           */
     { 0, 0,  15, 130,  15 },    /*  8  LEU, VAL, ILE */
     { 0, 0,  50,  50, 170 },    /*  9  PHE, TYR      */
     { 0, 0, 180,  90, 180 },    /* 10  TRP           */
     { 0, 0, 220, 150, 130 },    /* 11  PRO, PCA, HYP */
     { 0, 0, 190, 160, 110 } };  /* 12  Others        */

static int AminoIndex[] = {
      6, /*ALA*/   7, /*GLY*/   8, /*LEU*/   3,  /*SER*/
      8, /*VAL*/   3, /*THR*/   1, /*LYS*/   0,  /*ASP*/
      8, /*ILE*/   4, /*ASN*/   0, /*GLU*/  11,  /*PRO*/
      1, /*ARG*/   9, /*PHE*/   4, /*GLN*/   9,  /*TYR*/
      2, /*HIS*/   5, /*CYS*/   5, /*MET*/  10,  /*TRP*/
      4, /*ASX*/   4, /*GLX*/  11, /*PCA*/  11   /*HYP*/
                          };

static ShadeRef HBondShade[] = {
     { 0, 0, 255, 255, 255 },    /* 0  Offset =  2   */
     { 0, 0, 255,   0, 255 },    /* 1  Offset =  3   */
     { 0, 0, 255,   0,   0 },    /* 2  Offset =  4   */
     { 0, 0, 255, 165,   0 },    /* 3  Offset =  5   */
     { 0, 0,   0, 255, 255 },    /* 4  Offset = -3   */
     { 0, 0,   0, 255,   0 },    /* 5  Offset = -4   */
     { 0, 0, 255, 255,   0 } };  /* 6  Others        */


static ShadeRef StructShade[] = {
     { 0, 0, 255, 255, 255 },    /* 0  Default     */
     { 0, 0, 255,   0, 128 },    /* 1  Alpha Helix */
     { 0, 0, 255, 200,   0 },    /* 2  Beta Sheet  */
     { 0, 0,  96, 128, 255 } };  /* 3  Turn        */

static ShadeRef PotentialShade[] = {
     { 0, 0, 255,   0,   0 },    /* 0  Red     25 < V       */
     { 0, 0, 255, 165,   0 },    /* 1  Orange  10 < V <  25 */
     { 0, 0, 255, 255,   0 },    /* 2  Yellow   3 < V <  10 */
     { 0, 0,   0, 255,   0 },    /* 3  Green    0 < V <   3 */
     { 0, 0,   0, 255, 255 },    /* 4  Cyan    -3 < V <   0 */
     { 0, 0,   0,   0, 255 },    /* 5  Blue   -10 < V <  -3 */
     { 0, 0, 160,  32, 240 },    /* 6  Purple -25 < V < -10 */
     { 0, 0, 255, 255, 255 } };  /* 7  White        V < -25 */


/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(ptr=group->alist;ptr;ptr=ptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext) 
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)

#define MatchChar(a,b)   (((a)=='#')||((a)==(b)))


static ShadeRef ScaleRef[LastShade];
static int MaskColour[MAXMASK];
static int MaskShade[MAXMASK];
static int ScaleCount;

static Real LastRX,LastRY,LastRZ;
static Real DScale,IScale;
static Real Zoom;



void DetermineClipping( void )
{
    register int temp;
    register int max;

    max = 0;
    if( DrawAtoms && (MaxAtomRadius>max) )  max = MaxAtomRadius;
    if( DrawBonds && (MaxBondRadius>max) )  max = MaxBondRadius;
       
    temp = ImageRadius + max;
    if( (YOffset>=temp) && (XOffset>=temp) && (YOffset+temp<YRange) )
    {   if( UseStereo )
        {   UseScreenClip = (XOffset+temp) >= (XRange>>1);
        } else UseScreenClip = (XOffset+temp) >= XRange;
    } else UseScreenClip = True;
}


void SetRadiusValue( int rad )
{
    register int irad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    irad = (int)(Scale*rad);
    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( ptr->flag & SelectFlag )
        {   if( irad>MaxAtomRadius )
                MaxAtomRadius = irad;
            ptr->flag |= SphereFlag;
            ptr->radius = rad;
            ptr->irad = irad;
            change = True;
        } else if( ptr->flag & SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void SetRadiusTemperature( void )
{
    register int rad,irad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( (ptr->flag&SelectFlag) && (ptr->temp>0) )
        {   rad = (5*ptr->temp)>>1;
            if( rad>750 ) rad = 750;

            irad = (int)(Scale*rad);
            if( irad>MaxAtomRadius )
                MaxAtomRadius = irad;
            ptr->flag |= SphereFlag;
            ptr->radius = rad;
            ptr->irad = irad;
            change = True;
        } else if( ptr->flag & SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void SetVanWaalRadius( void )
{
    register int rad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   rad = ElemVDWRadius(ptr->elemno);
            ptr->irad = (int)(Scale*rad);
            ptr->radius = rad;
            change = True;

            ptr->flag |=SphereFlag;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        } else if( ptr->flag&SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void DisableSpacefill( void )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database || !DrawAtoms )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    
    ForEachAtom
        if( !(ptr->flag&SelectFlag) )
        {   if( ptr->flag&SphereFlag )
            {   if( ptr->irad>MaxAtomRadius )
                    MaxAtomRadius = ptr->irad;
                DrawAtoms = True;
            }
        } else if( ptr->flag&SphereFlag )
            ptr->flag &= ~SphereFlag;

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


void EnableWireframe( int mask, int rad )
{
    register Bond __far *bptr;
    register int flag, irad;

    if( !Database )
        return;

    DrawBonds = False;
    MaxBondRadius = 0;
    irad = (int)(Scale*rad);

    ForEachBond
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   DrawBonds = True;
            bptr->flag &= ~DrawBondFlag;
            bptr->flag |= mask;
            if( mask == CylinderFlag )
            {   if( irad>MaxBondRadius )
                    MaxBondRadius = irad;
                bptr->radius = rad;
                bptr->irad = irad;
            }
        } else if( bptr->flag&DrawBondFlag )
        {    DrawBonds = True;
             if( bptr->flag&CylinderFlag )
                 if( bptr->irad>MaxBondRadius )
                     MaxBondRadius = bptr->irad;
        }
    }
    DetermineClipping();
}


void DisableWireframe( void )
{
    register Bond __far *bptr;
    register int flag;

    if( !Database || !DrawBonds )
        return;

    DrawBonds = False;
    MaxBondRadius = 0;

    ForEachBond
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   bptr->flag &= ~DrawBondFlag;
        } else if( bptr->flag&DrawBondFlag )
        {   DrawBonds = True;
            if( bptr->flag&CylinderFlag )
                if( bptr->irad>MaxBondRadius )
                    MaxBondRadius = bptr->irad;
        }
    }
    DetermineClipping();
}


void EnableBackbone( int mask, int rad )
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register int flag,irad;

    if( !Database )
        return;

    irad = (int)(Scale*rad);

    ForEachBack
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   bptr->flag &= ~DrawBondFlag;
            bptr->flag |= mask;
            if( mask == CylinderFlag )
            {   bptr->radius = rad;
                bptr->irad = irad;
            }
        } 
    }
}


void DisableBackbone( void )
{
    register Chain __far *chain;
    register Bond __far *bptr;

    if( !Database )
        return;

    if( ZoneBoth )
    {   ForEachBack
            if( (bptr->dstatom->flag&bptr->srcatom->flag) & SelectFlag )
                bptr->flag &= ~DrawBondFlag;
    } else ForEachBack
        if( (bptr->dstatom->flag|bptr->srcatom->flag) & SelectFlag )
            bptr->flag &= ~DrawBondFlag;
}


void SetHBondStatus( int hbonds, int enable, int rad )
{
    register HBond __far *list;
    register HBond __far *ptr;
    register Atom __far *src;
    register Atom __far *dst;
    register int flag, irad;

    if( !Database )
        return;

    if( hbonds )
    {   if( enable && (Info.hbondcount<0) )
            CalcHydrogenBonds();
        list = Database->hlist;
    } else 
    {   if( enable && (Info.ssbondcount<0) )
            FindDisulphideBridges();
        list = Database->slist;
    }

    irad = (int)(Scale*rad);
    for( ptr=list; ptr; ptr=ptr->hnext )
    {   src = ptr->src;  dst = ptr->dst;

        flag = ZoneBoth? src->flag&dst->flag : src->flag|dst->flag;
        if( flag & SelectFlag ) 
        {   ptr->flag &= ~DrawBondFlag;
            if( enable )
            {   if( rad )
                {   ptr->flag |= CylinderFlag;
                    ptr->radius = rad;
                    ptr->irad = irad;
                } else ptr->flag |= WireFlag;
            }
        }
    }
}


void SetRibbonStatus( int enable, int flag, int width )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    /* Ribbons already disabled! */
    if( !enable && !DrawRibbon )
        return;

    if( Info.helixcount < 0 )
        DetermineStructure(False);

    DrawRibbon = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
            if( enable )
            {   if( group->flag & DrawKnotFlag )
                    DrawRibbon = True;
                
                for( ptr=group->alist; ptr; ptr=ptr->anext )
                    if( IsAlphaCarbon(ptr->refno) )
                    {   if( ptr->flag&SelectFlag )
                        {   group->flag &= ~DrawKnotFlag;
                            group->flag |= flag;
                            if( !width )
                            {   if( group->struc & (HelixFlag|SheetFlag) )
                                {      group->width = 380;
                                } else group->width = 100;
                            } else group->width = width;
                            DrawRibbon = True;
                        }
                        break;

                    } else if( IsSugarPhosphate(ptr->refno) )
                    {   if( ptr->flag&SelectFlag )
                        {   group->width = width? width : 720;
                            group->flag &= ~DrawKnotFlag;
                            group->flag |= flag;
                            DrawRibbon = True;
                        }
                        break;
                    }


            } else  /* Disable Ribbon */
                if( group->flag & DrawKnotFlag )
                {   for( ptr=group->alist; ptr; ptr=ptr->anext )
                        if( IsAlphaCarbon(ptr->refno) ||
                            IsSugarPhosphate(ptr->refno) )
                        {   if( ptr->flag&SelectFlag )
                                group->flag &= ~DrawKnotFlag;
                            break;
                        }
                    if( group->flag & DrawKnotFlag ) 
                        DrawRibbon = True;
                }
}


void SetRibbonCartoons( void )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    if( Info.helixcount < 0 )
        DetermineStructure(False);

    /* DrawBetaArrows = True; */
    /* CartoonHeight = 120;   */

    DrawRibbon = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
        {   if( group->flag & DrawKnotFlag )
                DrawRibbon = True;
                
            for( ptr=group->alist; ptr; ptr=ptr->anext )
                if( IsAlphaCarbon(ptr->refno) )
                {   if( ptr->flag&SelectFlag )
                    {   group->flag &= ~DrawKnotFlag;
                        if( group->struc & (HelixFlag|SheetFlag) )
                        {   group->flag |= CartoonFlag;
                            group->width = 380;
                        } else 
                        {   group->flag |= TraceFlag;
                            group->width = 100;
                        }
                        DrawRibbon = True;
                    }
                    break;

                } else if( IsSugarPhosphate(ptr->refno) )
                {   if( ptr->flag&SelectFlag )
                    {   group->flag &= ~DrawKnotFlag;
                        group->flag |= RibbonFlag;
                        group->width = 720;
                        DrawRibbon = True;
                    }
                    break;
                }
            }
}


void SetTraceTemperature( void )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int init,flag;
    register int min,max;
    register Real coeff;

    if( !Database )
        return;

    /* Avoid Compiler Warnings! */
    min = max = 0;

    flag = 0;
    init = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
            for( ptr=group->alist; ptr; ptr=ptr->anext )
                if( IsAlphaCarbon(ptr->refno) || IsSugarPhosphate(ptr->refno) )
                {   flag |= ptr->flag;
                    if( init )
                    {   if( ptr->temp<min ) 
                        {   min = ptr->temp;
                        } else if( ptr->temp>max )
                            max = ptr->temp;
                    } else
                    {   min = max = ptr->temp;
                        init = True;
                    }
                    break;
                }

    /* No groups selected! */
    if( !(flag&SelectFlag) )
        return;

    if( Info.helixcount < 0 )
        DetermineStructure(False);

    if( max != min )
    {   coeff = 200.0/(max-min);
    } else coeff = 0.0;

    DrawRibbon = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
        {   if( group->flag & DrawKnotFlag )
                DrawRibbon = True;
                
            for( ptr=group->alist; ptr; ptr=ptr->anext )
                if( IsAlphaCarbon(ptr->refno) || IsSugarPhosphate(ptr->refno) )
                {   if( ptr->flag&SelectFlag )
                    {   group->width = (int)(coeff*(ptr->temp-min))+50;
                        group->flag &= ~DrawKnotFlag;
                        group->flag |= TraceFlag;
                        DrawRibbon = True;
                    }
                    break;
                }
        }
}



/*===========================*/
/* Atom Selection Functions! */
/*===========================*/

static void DisplaySelectCount( void )
{
    char buffer[40];

    if( FileDepth == -1 )
    {   InvalidateCmndLine();
        if( SelectCount==0 )
        {   WriteString("No atoms selected!\n");
        } else if( SelectCount>1 )
        {   sprintf(buffer,"%ld atoms selected!\n",(long)SelectCount);
            WriteString(buffer);
        } else WriteString("1 atom selected!\n");
    }

    if( DisplayMode )
        ReDrawFlag |= RFRefresh;
    AdviseUpdate(AdvSelectCount);
}


void SelectZone( int mask )
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    SelectCount = 0;
    ForEachAtom
        if( ptr->flag & mask )
        {   ptr->flag |= SelectFlag;
            SelectCount++;
        } else ptr->flag &= ~SelectFlag;
    DisplaySelectCount();

    if( ZoneBoth )
    {   ForEachBond
           if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
    } else
        ForEachBond
           if( (bptr->srcatom->flag|bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;

}


void RestrictZone( int mask )
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int flag;

    if( !Database )
        return;

    DrawAtoms = False;   MaxAtomRadius = 0;
    DrawBonds = False;   MaxBondRadius = 0;
    
    SelectCount = 0;
    ForEachAtom
        if( ptr->flag & mask )
        {   ptr->flag |= SelectFlag;
            SelectCount++;

            if( ptr->flag & SphereFlag )
            {   DrawAtoms = True;
                if( ptr->irad>MaxAtomRadius )
                    MaxAtomRadius = ptr->irad;
            }
        } else 
        {   ptr->flag &= ~(SelectFlag|SphereFlag);
            if( ptr->label )
            {   DeleteLabel( (Label*)ptr->label );
                ptr->label = (void*)0;
            }
        }
    DisplaySelectCount();
    
    ForEachBond
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( flag & SelectFlag )
        {   bptr->flag |= SelectFlag;
            if( bptr->flag&DrawBondFlag )
            {   DrawBonds = True;
                if( bptr->flag & CylinderFlag )
                    if( bptr->irad>MaxBondRadius )
                        MaxBondRadius = bptr->irad;
            } 
        } else bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    ForEachBack
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( !(flag&SelectFlag) )
            bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    if( DrawRibbon )
    {   DrawRibbon = False;
        for( chain=Database->clist; chain; chain=chain->cnext )
            for( group=chain->glist; group; group=group->gnext )
                if( group->flag & DrawKnotFlag )
                {   for( ptr=group->alist; ptr; ptr=ptr->anext )
                        if( IsAlphaCarbon(ptr->refno) ||
                            IsSugarPhosphate(ptr->refno) )
                        {   if( !(ptr->flag&SelectFlag) )
                                group->flag &= ~DrawKnotFlag;
                            break;
                        }
                    if( group->flag & DrawKnotFlag ) 
                        DrawRibbon = True;
                }
    }

    RestrictMonitors();
    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


void SelectZoneExpr( Expr *expr )
{
    register Bond __far *bptr;

    if( !Database )
        return;

    SelectCount = 0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
        for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
            for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
                if( EvaluateExpr(expr) )
                {   QAtom->flag |= SelectFlag;
                    SelectCount++;
                } else QAtom->flag &= ~SelectFlag;
    DisplaySelectCount();

    if( ZoneBoth )
    {   ForEachBond
           if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
    } else
        ForEachBond
           if( (bptr->srcatom->flag|bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
}


void RestrictZoneExpr( Expr *expr )
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int flag;

    if( !Database )
        return;

    DrawAtoms = False;   MaxAtomRadius = 0;
    DrawBonds = False;   MaxBondRadius = 0;

    SelectCount = 0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
        for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
            for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
                if( EvaluateExpr(expr) )
                {   QAtom->flag |= SelectFlag;
                    SelectCount++;

                    if( QAtom->flag & SphereFlag )
                    {   DrawAtoms = True;
                        if( QAtom->irad>MaxAtomRadius )
                            MaxAtomRadius = QAtom->irad;
                    }
                }  else 
                {   QAtom->flag &= ~(SelectFlag|SphereFlag);
                    if( QAtom->label )
                    {   DeleteLabel( (Label*)QAtom->label );
                        QAtom->label = (void*)0;
                    }
                }
    DisplaySelectCount();

    ForEachBond
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( flag & SelectFlag )
        {   bptr->flag |= SelectFlag;
            if( bptr->flag & CylinderFlag )
            {   DrawBonds = True;
                if( bptr->irad>MaxBondRadius )
                    MaxBondRadius = bptr->irad;
            } else if( bptr->flag&WireFlag )
                DrawBonds = True;
        } else bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    ForEachBack
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( !(flag&SelectFlag) )
            bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    if( DrawRibbon )
    {   DrawRibbon = False;
        for( chain=Database->clist; chain; chain=chain->cnext )
            for( group=chain->glist; group; group=group->gnext )
                if( group->flag & DrawKnotFlag )
                {   for( ptr=group->alist; ptr; ptr=ptr->anext )
                        if( IsAlphaCarbon(ptr->refno) ||
                            IsSugarPhosphate(ptr->refno) )
                        {   if( !(ptr->flag&SelectFlag) )
                                group->flag &= ~DrawKnotFlag;
                            break;
                        }
                    if( group->flag & DrawKnotFlag )
                        DrawRibbon = True;
                }
    }

    RestrictMonitors();
    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


int DefineShade( int r, int g, int b )
{
    register int d,dr,dg,db;
    register int dist,best;
    register int i;

    /* Already defined! */
    for( i=0; i<LastShade; i++ )
        if( Shade[i].refcount )
            if( (Shade[i].r==r)&&(Shade[i].g==g)&&(Shade[i].b==b) )
                return i;

    /* Allocate request */
    for( i=0; i<LastShade; i++ )
         if( !Shade[i].refcount )
         {   Shade[i].r = r;
             Shade[i].g = g;
             Shade[i].b = b;
             Shade[i].refcount = 0;
             return i;
         }

    InvalidateCmndLine();
    WriteString("Warning: Unable to allocate shade!\n");

    /* To avoid lint warning! */
    best = dist = 0;

    /* Nearest match */
    for( i=0; i<LastShade; i++ )
    {   dr = Shade[i].r - r;
        dg = Shade[i].g - g;
        db = Shade[i].b - b;
        d = dr*dr + dg*dg + db*db;
        if( !i || (d<dist) )
        {   dist = d;
            best = i;
        }
    }
    return best;
}


static void ScaleColourMap( int count )
{
    register int i, r, g, b;
    register int fract;

    ScaleCount=0;
    for( i=0; i<LastShade; i++ )
        if( !Shade[i].refcount )
            ScaleCount++;

    /* If there are no shades free! */
    if( !ScaleCount ) ScaleCount = LastShade;

    if( count && (count<ScaleCount) )
        ScaleCount = count;

    if( ScaleCount == 1 )
    {   ScaleRef[i].r = 0;
        ScaleRef[i].g = 0;
        ScaleRef[i].b = 255;
        ScaleRef[i].shade = 0;
        ScaleRef[i].col = 0;
        return;
    }
    
    for( i=0; i<ScaleCount; i++ )
    {   fract = (int)((1023*i)/(ScaleCount-1));
        if( fract < 256 )
        {   r = 0;  g = fract;  b = 255;
        } else if( fract < 512 )
        {   r = 0;  g = 255;  b = 511-fract;
        } else if( fract < 768 )
        {   r = fract-512;  g = 255;  b = 0;
        } else /* fract < 1024 */                             
        {   r = 255;  g = 1023-fract;  b = 0;
        }
        ScaleRef[i].r = r;
        ScaleRef[i].g = g;
        ScaleRef[i].b = b;
        ScaleRef[i].shade = 0;
        ScaleRef[i].col = 0;
    }
}


static void SetLutEntry( int i, int r, int g, int b )
{
    ULut[i] = True;
    RLut[i] = r;
    GLut[i] = g;
    BLut[i] = b;

#ifdef EIGHTBIT
    Lut[i] = i;
#else
    Lut[i] = ( (Card)((r<<8)|g)<<8 ) | b;
#endif
}


static Real Power( Real x, int y )
{
    register Real result;

    result = x;
    while( y > 1 )
    {   if( y & 1 )
            result *= x;
        result *= result; 
        y >>= 1;
    }
    return result;
}


void DefineColourMap( void )
{
    register Real diffuse,fade;
    register Real temp,inten;
    register int col,r,g,b;
    register int i,j,k;

    /* Avoid compiler warning! */
    k = 0;

    for( i=0; i<LutSize; i++ )
        ULut[i] = False;

    if( !DisplayMode )
    {   SetLutEntry(BackCol,BackR,BackG,BackB);
        SetLutEntry(LabelCol,LabR,LabG,LabB);
        SetLutEntry(BoxCol,BoxR,BoxG,BoxB);
    } else SetLutEntry(BackCol,80,80,80);

    diffuse = 1.0 - Ambient;
    if( DisplayMode )
    {   for( i=0; i<ColourDepth; i++ )
        {   temp = (Real)i/ColourMask;
            inten = diffuse*temp + Ambient;

            /* Unselected [40,40,255] */
            /* Selected   [255,160,0]  */
            r = (int)(255*inten);
            g = (int)(160*inten);
            b = (int)(40*inten);

            /* Avoid Borland Compiler Warning! */
            /* Shade2Colour(0) == FirstCol     */
            SetLutEntry( FirstCol+i, b, b, r );
            SetLutEntry( Shade2Colour(1)+i, r, g, 0 );
        }
    } else
        for( i=0; i<ColourDepth; i++ )
        {   temp = (Real)i/ColourMask;
            inten = diffuse*temp + Ambient;
            fade = 1.0-inten;

            if( FakeSpecular )
            {   temp = Power(temp,SpecPower);
                k = (int)(255*temp);
                temp = 1.0 - temp;
                inten *= temp;
                fade *= temp;
            }

            for( j=0; j<LastShade; j++ )
                if( Shade[j].refcount )
                {   col = Shade2Colour(j);
                    if( UseBackFade )
                    {   temp = 1.0-inten;
                        r = (int)(Shade[j].r*inten + fade*BackR); 
                        g = (int)(Shade[j].g*inten + fade*BackG);
                        b = (int)(Shade[j].b*inten + fade*BackB);
                    } else
                    {   r = (int)(Shade[j].r*inten); 
                        g = (int)(Shade[j].g*inten);
                        b = (int)(Shade[j].b*inten);
                    }

                    if( FakeSpecular )
                    {   r += k;
                        g += k;
                        b += k;
                    }
                    SetLutEntry( col+i, r, g, b );
                }
        }

    if( Interactive )
        AllocateColourMap();
}


void ResetColourMap( void )
{
    register int i;

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
        ULut[i] = False;
#endif

    SpecPower = 8;
    FakeSpecular = False;
    Ambient = DefaultAmbient;
    UseBackFade = False;

    BackR = BackG = BackB = 0;
    BoxR = BoxG = BoxB = 255;
    LabR = LabG = LabB = 255;

    for( i=0; i<LastShade; i++ )
        Shade[i].refcount = 0;
    ScaleCount = 0;
}


void ColourBondNone( void )
{
    register Bond __far *bptr;

    if( Database )
        ForEachBond
            if( (bptr->flag&SelectFlag) && bptr->col )
            {   Shade[Colour2Shade(bptr->col)].refcount--;
                bptr->col = 0;
            }
}


void ColourBondAttrib( int r, int g, int b )
{
    register Bond __far *bptr;
    register int shade,col;

    if( Database )
    {   ForEachBond
            if( (bptr->flag&SelectFlag) && bptr->col )
                Shade[Colour2Shade(bptr->col)].refcount--;

        shade = DefineShade(r,g,b);
        col = Shade2Colour(shade);

        ForEachBond
            if( bptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                bptr->col = col;
            }
    }
}


void ColourBackNone( void )
{
    register Chain __far *chain;
    register Bond __far *bptr;
    register int flag;

    if( Database )
        ForEachBack
        {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                           : bptr->dstatom->flag | bptr->srcatom->flag;

            if( flag&SelectFlag )
            {   bptr->flag |= SelectFlag;
                if( bptr->col )
                {   Shade[Colour2Shade(bptr->col)].refcount--;
                    bptr->col = 0;
                }
            } else bptr->flag &= ~SelectFlag;
        }
}


void ColourBackAttrib( int r, int g, int b )
{
    register int shade,col;
    register Chain __far *chain;
    register Bond __far *bptr;

    if( Database )
    {   ColourBackNone();
        shade = DefineShade(r,g,b);
        col = Shade2Colour(shade);

        ForEachBack
            if( bptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                bptr->col = col;
            }
    }
}


void ColourHBondNone( int hbonds )
{
    register HBond __far *list;
    register HBond __far *ptr;
    register Atom __far *src;
    register Atom __far *dst;

    if( !Database )
        return;

    list = hbonds? Database->hlist : Database->slist;

    if( ZoneBoth )
    {   for( ptr=list; ptr; ptr=ptr->hnext )
        {   src = ptr->src;  dst = ptr->dst;

            if( (src->flag&dst->flag) & SelectFlag )
            {   ptr->flag |= SelectFlag;
                if( ptr->col )
                {   Shade[Colour2Shade(ptr->col)].refcount--;
                    ptr->col = 0;
                }
            } else ptr->flag &= ~SelectFlag;
        }
    } else
        for( ptr=list; ptr; ptr=ptr->hnext )
        {   src = ptr->src;  dst = ptr->dst;

            if( (src->flag|dst->flag) & SelectFlag )
            {   ptr->flag |= SelectFlag;
                if( ptr->col )
                {   Shade[Colour2Shade(ptr->col)].refcount--;
                    ptr->col = 0;
                }
            } else ptr->flag &= ~SelectFlag;
        }
}


void ColourHBondType( void )
{
    register HBond __far *ptr;
    register ShadeRef *ref;
    register int i;

    if( !Database ) return;
    for( i=0; i<7; i++ )
        HBondShade[i].col = 0;

    if( Info.hbondcount < 0 )
    {   CalcHydrogenBonds();
    } else ColourHBondNone( True );

    for( ptr=Database->hlist; ptr; ptr=ptr->hnext )
        if( ptr->flag & SelectFlag )
        {   switch( ptr->offset )
            {   case(  2 ):  ref = HBondShade;     break;
                case(  3 ):  ref = HBondShade+1;   break;
                case(  4 ):  ref = HBondShade+2;   break;
                case(  5 ):  ref = HBondShade+3;   break;
                case( -3 ):  ref = HBondShade+4;   break;
                case( -4 ):  ref = HBondShade+5;   break;
                default:     ref = HBondShade+6;   break;
            }

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = (unsigned char)ref->col;
        }
}


void ColourHBondAttrib( int hbonds, int r, int g, int b )
{
    register HBond __far *list;
    register HBond __far *ptr;
    register int col,shade;

    if( !Database )
        return;

    if( hbonds )
    {   if( Info.hbondcount < 0 )
        {   CalcHydrogenBonds();
        } else ColourHBondNone(True);
    } else
        if( Info.ssbondcount < 0 )
        {   FindDisulphideBridges();
        } else ColourHBondNone(False);


    shade = DefineShade(r,g,b);
    col = Shade2Colour(shade);

    list = hbonds? Database->hlist : Database->slist;
    for( ptr=list; ptr; ptr=ptr->hnext )
        if( ptr->flag & SelectFlag )
        {   Shade[shade].refcount++;
            ptr->col = col;
        }
}


void ColourRibbonNone( int flag )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    if( !Database )
        return;

    if( Info.helixcount < 0 )
        return;

    for( chain=Database->clist; chain; chain=chain->cnext )
        for( group=chain->glist; group; group=group->gnext )
            for( aptr=group->alist; aptr; aptr=aptr->anext )
                if( (aptr->flag&SelectFlag) && 
                    (IsAlphaCarbon(aptr->refno)||
                     IsSugarPhosphate(aptr->refno)) )
                {   if( (flag&RibColInside) && group->col1 )
                    {   Shade[Colour2Shade(group->col1)].refcount--;
                        group->col1 = 0;
                    }
                    if( (flag&RibColOutside) && group->col2 )
                    {   Shade[Colour2Shade(group->col2)].refcount--;
                        group->col2 = 0;
                    }
                    break;
                }
}


void ColourRibbonAttrib( int flag, int r, int g, int b )
{
    register int shade, col;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    if( Database )
    {   if( Info.helixcount >= 0 )
        {   ColourRibbonNone( flag );
        } else DetermineStructure(False);

        shade = DefineShade(r,g,b);
        col = Shade2Colour(shade);

        for( chain=Database->clist; chain; chain=chain->cnext )
            for( group=chain->glist; group; group=group->gnext )
                for( aptr=group->alist; aptr; aptr=aptr->anext )
                    if( (aptr->flag&SelectFlag) && 
                        (IsAlphaCarbon(aptr->refno)||
                         IsSugarPhosphate(aptr->refno)) )
                    {   if( flag & RibColInside )
                        {   Shade[shade].refcount++;
                            group->col1 = col;
                        }
                        if( flag & RibColOutside )
                        {   Shade[shade].refcount++;
                            group->col2 = col;
                        }
                        break;
                    }
    }
}


void ColourMonitNone( void )
{
    register Monitor *ptr;
    register int flag;

    if( Database )
        for( ptr=MonitList; ptr; ptr=ptr->next )
            if( ptr->col )
            {   flag = ZoneBoth? ptr->src->flag & ptr->dst->flag
                               : ptr->src->flag | ptr->dst->flag;
                if( flag & SelectFlag )
                {   Shade[Colour2Shade(ptr->col)].refcount--;
                    ptr->col = 0;
                }
            }
}


void ColourMonitAttrib( int r, int g, int b )
{
    register Monitor *ptr;
    register int shade,col;
    register int flag;

    if( !Database )
        return;

    ColourMonitNone();
    shade = DefineShade(r,g,b);
    col = Shade2Colour(shade);

    for( ptr=MonitList; ptr; ptr=ptr->next )
    {   flag = ZoneBoth? ptr->src->flag & ptr->dst->flag 
                       : ptr->src->flag | ptr->dst->flag;
        if( flag & SelectFlag )
        {   Shade[shade].refcount++;
            ptr->col = col;
        }
    }
}


void ColourDotsAttrib( int r, int g, int b )
{
    register DotStruct __far *ptr;
    register int i,shade,col;

    if( Database )
    {   for( ptr=DotPtr; ptr; ptr=ptr->next )
            for( i=0; i<ptr->count; i++ )
            {    shade = Colour2Shade(ptr->col[i]);
                 Shade[shade].refcount--;
            }

        shade = DefineShade(r,g,b);
        col = Shade2Colour(shade);
        for( ptr=DotPtr; ptr; ptr=ptr->next )
            for( i=0; i<ptr->count; i++ )
            {   Shade[shade].refcount++;
                ptr->col[i] = col;
            }
    }
}


/* Coulomb's Law */
#define CoulombScale  ((Long)(1<<12))
static int CalculatePotential( Long x, Long y, Long z )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register Long dx,dy,dz;
    register Long result;
    register Card dist;
    register Card max;


    /* Calculated charges have b-values < 0.0     */
    /* if( MinFun(MinMainTemp,MinHetaTemp) >= 0 ) */
    /*     CalculateCharges();                    */

    /* 8.0 Angstrom Cut Off */
    max = (Long)2000*2000;

    result = 0;
    ForEachAtom
    {   dx = ptr->xorg-x;
        if( (dist=dx*dx) < max )
        {   dy = ptr->yorg - y;
            if( (dist+=dy*dy) < max )
            {   dz = ptr->zorg - z;
                if( (dist+=dz*dz) < max )
                    result += (CoulombScale*ptr->temp) / (int)isqrt(dist);
            }
        }
    }
    /* Dielectric Constant = 10.0 */
    /* (332.0*250.0)/(10.0*100.0) */
    result = (result*83)/CoulombScale;
    return (int)result;
}


void ColourDotsPotential( void )
{
    register DotStruct __far *ptr;
    register int i,shade,result;
    register ShadeRef *ref;

    if( Database )
    {   for( i=0; i<8; i++ )
            PotentialShade[i].col = 0;

        /* Colour Dots None! */
        for( ptr=DotPtr; ptr; ptr=ptr->next )
            for( i=0; i<ptr->count; i++ )
            {    shade = Colour2Shade(ptr->col[i]);
                 Shade[shade].refcount--;
            }

        for( ptr=DotPtr; ptr; ptr=ptr->next )
            for( i=0; i<ptr->count; i++ )
            {   result = CalculatePotential( ptr->xpos[i],
                                             ptr->ypos[i],
                                             ptr->zpos[i] );

                /* Determine Colour Bucket */
                if( result >= 0 )
                {   if( result > 10 )
                    {      if( result > 24 )
                           {      ref = PotentialShade + 0;
                           } else ref = PotentialShade + 1;
                    } else if( result > 3 )
                           {      ref = PotentialShade + 2;
                           } else ref = PotentialShade + 3;
                } else 
                    if( result > -10 )
                    {      if( result > -3 )
                           {      ref = PotentialShade + 4;
                           } else ref = PotentialShade + 5;
                    } else if( result > -24 )
                           {      ref = PotentialShade + 6;
                           } else ref = PotentialShade + 7;

                if( !ref->col )
                {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                    ref->col = Shade2Colour(ref->shade);
                }
                Shade[ref->shade].refcount++;
                ptr->col[i] = ref->col;
            }
    }
}


static void ResetColourAttrib( void )
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    ForEachAtom
        if( (ptr->flag&SelectFlag) && ptr->col )
            Shade[Colour2Shade(ptr->col)].refcount--;
}


void MonoColourAttrib( int r, int g, int b )
{
    register int shade,col;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( Database )
    {   ResetColourAttrib();
        shade = DefineShade(r,g,b);
        col = Shade2Colour(shade);

        ForEachAtom
            if( ptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                ptr->col = col;
            }
    }
}


void ScaleColourAttrib( int attr )
{
    register ShadeRef *ref;
    register int count, attrno, factor;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register Long temp;
    register int i;

    if( !Database ) return;

    switch( attr )
    {   case(ChainAttr):   attrno = Info.chaincount;   
                           factor = 1;
                           break;

        case(GroupAttr):   factor = MinMainRes;
                           attrno = MaxMainRes;
                           if( HetaGroups && HetaGroupCount )
                           {   if( MinHetaRes < factor )
                                   factor = MinHetaRes;
                               if( MaxHetaRes > attrno )
                                   attrno = MaxHetaRes;
                           } 
                           attrno -= (factor-1);
                           break;

        case(ChargeAttr):
        case(TempAttr):    factor = MinMainTemp;
                           attrno = MaxMainTemp;
                           if( HetaGroups && HetaGroupCount )
                           {   if( MinHetaTemp < factor )
                                   factor = MinHetaTemp;
                               if( MaxHetaTemp > attrno )
                                   attrno = MaxHetaTemp;
                           }
                           attrno -= (factor-1);
                           break;

        default:           return;
    }

    if( attrno<2 )
    {   MonoColourAttrib(255,255,255);
        return;
    }

    ResetColourAttrib();
    ScaleColourMap(attrno);

    switch( attr )
    {    case(ChainAttr):
                 count = 0;
                 for( chain=Database->clist; chain; chain=chain->cnext )
                 {   ref = &(ScaleRef[(count*ScaleCount)/attrno]);
                     if( !(ref->col && Shade[ref->shade].refcount) )
                     {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                         ref->col = Shade2Colour(ref->shade);
                     }
                     for( group=chain->glist; group; group=group->gnext )
                         for( ptr=group->alist; ptr; ptr=ptr->anext )
                             if( ptr->flag&SelectFlag )
                             {   Shade[ref->shade].refcount++;
                                 ptr->col = ref->col;
                             }
                     count++;
                 }
                 break;


         case(GroupAttr):
                 for( chain=Database->clist; chain; chain=chain->cnext )
                     for( group=chain->glist; group; group=group->gnext )
                     {   temp = (Long)ScaleCount*(group->serno-factor);
                         i = (int)(temp/attrno);

                         if( i >= ScaleCount )
                         {   ref = ScaleRef + (ScaleCount-1);
                         } else if( i >= 0 )
                         {   ref = ScaleRef + i;
                         } else ref = ScaleRef;

                         if( !(ref->col && Shade[ref->shade].refcount) )
                         {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                             ref->col = Shade2Colour(ref->shade);
                         }

                         for( ptr=group->alist; ptr; ptr=ptr->anext )
                             if( ptr->flag&SelectFlag )
                             {   Shade[ref->shade].refcount++;
                                 ptr->col = ref->col;
                             }
                     }
                 break;


         case(TempAttr):
                 ForEachAtom
                     if( ptr->flag&SelectFlag )
                     {   i = (int)(((Long)ScaleCount*(ptr->temp-factor))
                                    /attrno);

                         if( i >= ScaleCount )
                         {   ref = ScaleRef + (ScaleCount-1);
                         } else if( i >= 0 )
                         {   ref = ScaleRef + i;
                         } else ref = ScaleRef;

                         if( !(ref->col && Shade[ref->shade].refcount) )
                         {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                             ref->col = Shade2Colour(ref->shade);
                         }
                         Shade[ref->shade].refcount++;
                         ptr->col = ref->col;
                     }
                 break;

        case(ChargeAttr):
                ForEachAtom
                     if( ptr->flag&SelectFlag )
                     {   i = (int)(((Long)ScaleCount*(ptr->temp-factor))
                                    /attrno);

                         if( i <= 0 )
                         {   ref = ScaleRef + (ScaleCount-1);
                         } else if( i < ScaleCount )
                         {   ref = ScaleRef + ((ScaleCount-1)-i);
                         } else ref = ScaleRef;

                         if( !(ref->col && Shade[ref->shade].refcount) )
                         {   ref->shade = DefineShade(ref->r,ref->g,ref->b);
                             ref->col = Shade2Colour(ref->shade);
                         }
                         Shade[ref->shade].refcount++;
                         ptr->col = ref->col;
                     }
                 break;
    }
}



/*====================================*/
/*  Raster3D Color Record Processing  */
/*====================================*/

static int MatchNumber( int len, int value, char *mask )
{
    register int result;
    register int digit;
    register int ch;
    register int i;

    result = True;
    for( i=0; i<len; i++ )
    {   digit = (value%10) + '0';
        ch = mask[len-i];
        if( ch == ' ' )
        {   if( value ) result = False;
        } else if( !MatchChar(ch,digit) )
            result = False;
        value /= 10;
    }
    return result;
}


void UserMaskAttrib( int fields )
{
    register MaskDesc *mptr;
    register char *temp, *name;
    register int shade, change;
    register int i, rad, match;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database ) return;

    if( !MaskCount )
    {   InvalidateCmndLine();
        WriteString("Warning: No user supplied colour records!\n");
        return;
    }

    /* Avoid Compiler Warning! */
    mptr = (MaskDesc*)0;
    match = False;

    change = False;
    ResetColourAttrib();
    if( fields & MaskColourFlag )
        for( i=0; i<MaskCount; i++ )
            MaskShade[i] = -1;

    if( fields & MaskRadiusFlag )
    {   MaxAtomRadius = 0;
        DrawAtoms = False;
    }

    ForEachAtom
    if( ptr->flag & SelectFlag )
    {   for( i=0; i<MaskCount; i++ )
        {   mptr = UserMask+i;
            temp = mptr->mask;
            match = True;

            if( !MatchChar(temp[13],chain->ident) ) match=False;
            if( !MatchChar(temp[9],ptr->altl) )     match=False;

            /* Atom Name */
            if( match )
            {   name = ElemDesc[ptr->refno];
                if( !MatchChar(temp[5],name[0]) ) match=False;
                if( !MatchChar(temp[6],name[1]) ) match=False;
                if( !MatchChar(temp[7],name[2]) ) match=False;
                if( !MatchChar(temp[8],name[3]) ) match=False;
            }

            /* Group Name */
            if( match )
            {   name = Residue[group->refno];
                if( !MatchChar(temp[10],name[0]) ) match=False;
                if( !MatchChar(temp[11],name[1]) ) match=False;
                if( !MatchChar(temp[12],name[2]) ) match=False;
            }

            if( match && (mptr->flags&SerNoFlag) )
                match = MatchNumber(4,ptr->serno,&temp[0]);
            if( match && (mptr->flags&ResNoFlag) )
                match = MatchNumber(3,group->serno,&temp[14]);
            if( match ) break;
        }

        if( fields & MaskColourFlag )
        {   if( match )
            {   if( MaskShade[i] == -1 )
                {   MaskShade[i] = DefineShade(mptr->r,mptr->g,mptr->b);
                    MaskColour[i] = Shade2Colour(MaskShade[i]);
                }
                Shade[MaskShade[i]].refcount++;
                ptr->col = MaskColour[i];
            } else
            {   shade = DefineShade(255,255,255);
                ptr->col = Shade2Colour(shade);
                Shade[shade].refcount++;
            }
        }

        if( fields & MaskRadiusFlag )
        {   rad = match? mptr->radius : 375;
            ptr->irad = (int)(Scale*rad);
            ptr->flag |= SphereFlag;
            ptr->radius = rad;

            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
            change = True;
        }
    } else if( ptr->flag & SphereFlag )
    {   DrawAtoms = True;
        if( ptr->irad>MaxAtomRadius )
            MaxAtomRadius = ptr->irad;
    }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void CPKColourAttrib( void )
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<CPKMAX; i++ )
        CPKShade[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   ref = CPKShade + Element[ptr->elemno].cpkcol;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void AminoColourAttrib( void )
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<13; i++ )
        AminoShade[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   if( IsAmino(group->refno) )
            {   ref = AminoShade + AminoIndex[group->refno];
            } else ref = AminoShade+12;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void ShapelyColourAttrib( void )
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<30; i++ )
        Shapely[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   if( IsAminoNucleo(group->refno) )
            {   ref = Shapely + group->refno;
            } else ref = Shapely+30;

/*  Original Colour Scheme
 *
 *  ref = &(Shapely[26]);
 *  if( IsNucleo(group->refno) )
 *  {   ref = Shapely + group->refno;
 *  } else if( IsShapelyBackbone(ptr->refno) )
 *  {   ref = &(Shapely[24]);
 *  } else if( IsShapelySpecial(ptr->refno) )
 *  {   ref = &(Shapely[25]);
 *  } else if( IsAmino(group->refno) )
 *      ref = Shapely + group->refno;
 */

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void StructColourAttrib( void )
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database )
        return;

    if( Info.helixcount < 0 )
        DetermineStructure(False);

    for( i=0; i<4; i++ )
        StructShade[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   if( group->struc & HelixFlag )
            {   ref = StructShade+1;
            } else if( group->struc & SheetFlag )
            {   ref = StructShade+2;
            } else if( group->struc & TurnFlag )
            {   ref = StructShade+3;
            } else ref = StructShade;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


int IsCPKColour( Atom __far *ptr )
{
    register ShadeRef *cpk;
    register ShadeDesc *col;

    cpk = CPKShade + Element[ptr->elemno].cpkcol;
    col = Shade + Colour2Shade(ptr->col);
    return( (col->r==cpk->r) && 
            (col->g==cpk->g) && 
            (col->b==cpk->b) );
}


int IsVDWRadius( Atom __far *ptr )
{
    register int rad;

    if( !(ptr->flag&SphereFlag) )
        return False;

    rad = ElemVDWRadius( ptr->elemno );
    return ptr->radius == rad;
}


void DefaultRepresentation( void )
{
    if( Database )
    {   ReDrawFlag |= RFRefresh | RFColour;
        if( Info.bondcount < 1 )
        {   EnableBackbone(CylinderFlag,80);
        } else EnableWireframe(WireFlag,0);
        CPKColourAttrib();
    }
}


void InitialTransform( void )
{
    register Card dist,max;
    register double fdist,fmax;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register Card ax, ay, az;
    register Long dx, dy, dz;

    dx = MaxX-MinX;   OrigCX = (dx>>1)+MinX;
    dy = MaxY-MinY;   OrigCY = (dy>>1)+MinY;
    dz = MaxZ-MinZ;   OrigCZ = (dz>>1)+MinZ;

    MaxX -= OrigCX;   MinX -= OrigCX;
    MaxY -= OrigCY;   MinY -= OrigCY;
    MaxZ -= OrigCZ;   MinZ -= OrigCZ;

    SideLen = MaxFun(dx,dy);
    if( dz>SideLen ) SideLen = dz;
    SideLen += 1500;  Offset = SideLen>>1;
    XOffset = WRange;  YOffset = HRange;
    ZOffset = 10000;

    ForEachAtom
    {   ptr->xorg -= OrigCX;
        ptr->yorg -= OrigCY;
        ptr->zorg -= OrigCZ;
    }

    if( Offset > 37836 )
    {   fmax = 0.0;
        ForEachAtom
        {   ax = (Card)AbsFun(ptr->xorg);
            ay = (Card)AbsFun(ptr->yorg);
            az = (Card)AbsFun(ptr->zorg);
            fdist = (double)ax*ax + 
                    (double)ay*ay + 
                    (double)az*az;
            if( fdist > fmax )
                fmax = fdist;
        }
    } else
    {   max = 1;
        ForEachAtom
        {   ax = (Card)AbsFun(ptr->xorg);
            ay = (Card)AbsFun(ptr->yorg);
            az = (Card)AbsFun(ptr->zorg);
            dist = ax*ax + ay*ay + az*az;
            if( dist > max )
                max = dist;
        }
        fmax = (double)max;
    }

    WorldRadius = ((Card)sqrt(fmax))+750;
    WorldSize = WorldRadius<<1;
    DScale = 1.0/WorldSize;

    /* Code should match ReSizeScreen() */
    /* MaxZoom*DScale*Range*750 == 252  */
    MaxZoom = 0.324*WorldSize/Range;
    if( MaxZoom < 1.0 )
    {   DScale *= MaxZoom;
        MaxZoom = 1.0;
    }
    ZoomRange = Range;
    MaxZoom -= 1.0;
}


void ReviseInvMatrix( void )
{
    /* The inverse of a rotation matrix
     * is its transpose, and the inverse
     * of Scale is 1.0/Scale [IScale]!
     */
    InvX[0] = IScale*RotX[0];
    InvX[1] = IScale*RotY[0];
    InvX[2] = IScale*RotZ[0];

    InvY[0] = IScale*RotX[1];
    InvY[1] = IScale*RotY[1];
    InvY[2] = IScale*RotZ[1];

    InvZ[0] = IScale*RotX[2];
    InvZ[1] = IScale*RotY[2];
    InvZ[2] = IScale*RotZ[2];
    ShadowTransform();
}


void TranslateToCentre( Long x, Long y, Long z )
{
    register int ix,iy;

    ix = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2]);
    iy = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2]);
    if( PickMode == PickOrign )
    {   DialValue[4] = (Real)-ix/(Zoom*XRange);
        DialValue[5] = (Real)-iy/(Zoom*YRange);
    } else
    {   DialValue[4] = (Real)-ix/XRange;
        DialValue[5] = (Real)-iy/YRange;
    }
    ReDrawFlag |= RFTransX | RFTransY;
}


void PrepareTransform( void )
{
    register Real theta, temp;
    register Real cost, sint;
    register Real x, y, z;
    register Real ncost;

    if( (ReDrawFlag&RFRotateX) && (DialValue[0]!=LastRX) )
    {   theta = PI*(DialValue[0]-LastRX);
        cost = cos(theta);  sint = sin(theta);
        LastRX = DialValue[0];

        y=RotY[0]; z=RotZ[0];
        RotY[0]=cost*y+sint*z; 
        RotZ[0]=cost*z-sint*y;

        y=RotY[1]; z=RotZ[1];
        RotY[1]=cost*y+sint*z;
        RotZ[1]=cost*z-sint*y;

        y=RotY[2]; z=RotZ[2];
        RotY[2]=cost*y+sint*z;
        RotZ[2]=cost*z-sint*y;

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotY[0] + sint*RotZ[0]);
            temp += CenY*(ncost*RotY[1] + sint*RotZ[1]);
            temp += CenZ*(ncost*RotY[2] + sint*RotZ[2]);
            if( PickMode == PickOrign )
            {   temp = DialValue[5] - (Scale*temp)/(Zoom*YRange);
            } else temp = DialValue[5] - (Scale*temp)/YRange;

            if( temp < -1.0 )
            {   DialValue[5] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[5] = 1.0;
            } else DialValue[5] = temp;
        }
    }

    if( (ReDrawFlag&RFRotateY) && (DialValue[1]!=LastRY) )
    {   theta = PI*(DialValue[1]-LastRY);
        cost = cos(theta);  sint = sin(theta);
        LastRY = DialValue[1];

        x=RotX[0]; z=RotZ[0];
        RotX[0]=cost*x+sint*z;
        RotZ[0]=cost*z-sint*x;

        x=RotX[1]; z=RotZ[1];
        RotX[1]=cost*x+sint*z;
        RotZ[1]=cost*z-sint*x;

        x=RotX[2]; z=RotZ[2];
        RotX[2]=cost*x+sint*z;
        RotZ[2]=cost*z-sint*x;

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotX[0] + sint*RotZ[0]);
            temp += CenY*(ncost*RotX[1] + sint*RotZ[1]);
            temp += CenZ*(ncost*RotX[2] + sint*RotZ[2]);
            if( PickMode == PickOrign )
            {   temp = DialValue[4] - (Scale*temp)/(Zoom*XRange);
            } else temp = DialValue[4] - (Scale*temp)/XRange;

            if( temp < -1.0 )
            {   DialValue[4] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[4] = 1.0;
            } else DialValue[4] = temp;
        }
    }

    if( (ReDrawFlag&RFRotateZ) && (DialValue[2]!=LastRZ) )
    {   theta = PI*(DialValue[2]-LastRZ);
        cost = cos(theta);  sint = sin(theta);
        LastRZ = DialValue[2];

        x=RotX[0]; y=RotY[0];
        RotX[0]=cost*x-sint*y;
        RotY[0]=cost*y+sint*x;

        x=RotX[1]; y=RotY[1];
        RotX[1]=cost*x-sint*y;
        RotY[1]=cost*y+sint*x;

        x=RotX[2]; y=RotY[2];
        RotX[2]=cost*x-sint*y;
        RotY[2]=cost*y+sint*x;

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotX[0] - sint*RotY[0]);
            temp += CenY*(ncost*RotX[1] - sint*RotY[1]);
            temp += CenZ*(ncost*RotX[2] - sint*RotY[2]);

            if( PickMode == PickOrign )
            {   temp = DialValue[4] - (Scale*temp)/(Zoom*XRange);
            } else temp = DialValue[4] - (Scale*temp)/XRange;

            if( temp < -1.0 )
            {   DialValue[4] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[4] = 1.0;
            } else DialValue[4] = temp;

            ncost = 1.0-cost;
            temp =  CenX*(ncost*RotY[0] + sint*RotX[0]);
            temp += CenY*(ncost*RotY[1] + sint*RotX[1]);
            temp += CenZ*(ncost*RotY[2] + sint*RotX[2]);

            if( PickMode == PickOrign )
            {   temp = DialValue[5] - (Scale*temp)/(Zoom*YRange);
            } else temp = DialValue[5] - (Scale*temp)/YRange;

            if( temp < -1.0 )
            {   DialValue[5] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[5] = 1.0;
            } else DialValue[5] = temp;
        }
    }
}


void ApplyTransform( void )
{
    register int temp;
    register Real x, y, z;
    register int oldx,oldy;
    register Chain __far *chain;
    register Group __far *group;
    register HBond __far *hptr;
    register Bond __far *bptr;
    register Atom __far *ptr;


    if( ReDrawFlag & RFMagnify )
    {   if( DialValue[3] <= 0.0 )
        {   Zoom = DialValue[3]+1.0;
            if( Zoom<0.1 ) Zoom=0.1;
        } else Zoom = (DialValue[3]*MaxZoom) + 1.0;

        Scale = Zoom*DScale*Range;
        ImageSize = (int)(Scale*WorldSize);
        if( ImageSize < 2 )
        {   ImageRadius = 1;
            ImageSize = 2;
        } else 
            ImageRadius = ImageSize>>1;
        IScale = 1.0/Scale;

        MaxAtomRadius = 0;
        MaxBondRadius = 0;
    }

    if( ReDrawFlag & RFRotate )
    {   PrepareTransform();
        if( UseShadow )
            ShadowTransform();
    }

    if( ReDrawFlag & (RFRotate|RFMagnify) )
    {   MatX[0] = Scale*RotX[0]; 
        MatX[1] = Scale*RotX[1];
        MatX[2] = Scale*RotX[2];

        MatY[0] = Scale*RotY[0];
        MatY[1] = Scale*RotY[1];
        MatY[2] = Scale*RotY[2];

        MatZ[0] = Scale*RotZ[0];
        MatZ[1] = Scale*RotZ[1];
        MatZ[2] = Scale*RotZ[2];

        if( UseShadow )
        {   InvX[0] = IScale*RotX[0]; 
            InvX[1] = IScale*RotY[0];
            InvX[2] = IScale*RotZ[0];

            InvY[0] = IScale*RotX[1];
            InvY[1] = IScale*RotY[1];
            InvY[2] = IScale*RotZ[1];

            InvZ[0] = IScale*RotX[2];
            InvZ[1] = IScale*RotY[2];
            InvZ[2] = IScale*RotZ[2];
        }
    }

    oldx = XOffset;
    oldy = YOffset;

    /* Zoom dependent Translation! */
    /* XOffset = WRange + (int)(DialValue[4]*ImageSize); */
    /* YOffset = HRange + (int)(DialValue[5]*ImageSize); */

    if( PickMode == PickOrign )
    {   /* Zoom about Center of Screen */
        XOffset = WRange + (int)(Zoom*DialValue[4]*XRange);
        YOffset = HRange + (int)(Zoom*DialValue[5]*YRange);
    } else
    {   XOffset = WRange + (int)(DialValue[4]*XRange);
        YOffset = HRange + (int)(DialValue[5]*YRange);
    }

    if( UseStereo ) XOffset /= 2;

    switch( ReDrawFlag )
    {   case(RFTransX):
                if( XOffset != oldx ) 
                {   temp = XOffset - oldx;
                    ForEachAtom ptr->x += temp;
                }
                break;

        case(RFTransY):
                if( YOffset != oldy ) 
                {   temp = YOffset - oldy;
                    ForEachAtom ptr->y += temp;
                }
                break;

        case(RFRotateX):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
            }
            break;

        case(RFRotateY):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
            }
            break;

        case(RFRotateZ):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
            }
            break;

        default:
            /* This condition scales atomic radii! */
            if( DrawAtoms && (ReDrawFlag&RFMagnify) )
            {   ForEachAtom 
                {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                    ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                    ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                    ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
                    if( ptr->flag&SphereFlag )
                    {   ptr->irad = (int)(Scale*ptr->radius);
                        if( ptr->irad>MaxAtomRadius )
                            MaxAtomRadius = ptr->irad;
                    }
                }
            } else
                ForEachAtom 
                {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                    ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                    ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                    ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
                }

            if( ReDrawFlag & RFMagnify )
            {   if( DrawBonds )
                    ForEachBond
                        if( bptr->flag&CylinderFlag )
                        {   bptr->irad = (int)(Scale*bptr->radius);
                            if( bptr->irad>MaxBondRadius )
                            MaxBondRadius = bptr->irad;
                        }

                for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
                    if( hptr->flag&CylinderFlag )
                        hptr->irad = (int)(Scale*hptr->radius);

                for( hptr=Database->slist; hptr; hptr=hptr->hnext )
                    if( hptr->flag&CylinderFlag )
                        hptr->irad = (int)(Scale*hptr->radius);

                ForEachBack
                    if( bptr->flag&CylinderFlag )
                        bptr->irad = (int)(Scale*bptr->radius);
            }
    }

    DetermineClipping();
    if( UseScreenClip || ReDrawFlag!=RFRotateY )
        BucketFlag = False;
}


void ResetTransform( void )
{
    RotX[0] = 1.0;  RotX[1] = 0.0;  RotX[2] = 0.0;
    RotY[0] = 0.0;  RotY[1] = 1.0;  RotY[2] = 0.0;
    RotZ[0] = 0.0;  RotZ[1] = 0.0;  RotZ[2] = 1.0;
    LastRX = LastRY = LastRZ = 0.0;
    CenX = CenY = CenZ = 0;
}


void InitialiseTransform( void )
{
    ResetColourMap();
    ResetTransform();

    ZoneBoth = True;
    HetaGroups = True;
    Hydrogens = True;
}

