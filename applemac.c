/* applemac.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include <QuickDraw.h>
#include <Controls.h>
#include <Palettes.h>
#include <Windows.h>
#include <Errors.h>
#include <Menus.h>
#include <Fonts.h>

#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#include <ToolUtils.h>
#include <Memory.h>
#include <Types.h>
#include <Scrap.h>

#include <stdlib.h>
#include <stdio.h>

#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"



#ifdef __CONDITIONALMACROS__
#define ArrowCursor   SetCursor(&qd.arrow)
#else
#define ArrowCursor   SetCursor(&arrow)
#endif



static PixPatHandle BackHand;
static PaletteHandle CMap;
static PixMap *PixelMap;



int CreateImage( void )
{
    register Long size;
    
    if( FBuffer ) _ffree(FBuffer);
    size = (Long)XRange*YRange*sizeof(Pixel);
    FBuffer = (Pixel*)_fmalloc( size+32 );

#ifdef EIGHTBIT
    PixelMap->rowBytes = XRange | 0x8000;
#else
    PixelMap->rowBytes = (XRange<<2) | 0x8000;
#endif

    PixelMap->baseAddr = (Ptr)FBuffer;
    PixelMap->bounds.right = XRange;
    PixelMap->bounds.bottom = YRange;
    return (int)FBuffer;
}


void TransferImage( void )
{
    register PixMapHandle pmHand;
    register GDHandle gdHand;
    register int mode;
    GrafPtr savePort;
    Rect rect;
   
    GetPort(&savePort);
    SetPort(CanvWin);
   
    ForeColor(blackColor);
    BackColor(whiteColor);
   
    rect.top = 0;   rect.bottom = YRange;
    rect.left = 0;  rect.right = XRange;

    gdHand = GetMaxDevice(&CanvWin->portRect);
    pmHand = (**gdHand).gdPMap;
    if( (**pmHand).pixelSize < (sizeof(Pixel)<<3) )
    {   mode = srcCopy+ditherCopy;
    } else mode = srcCopy;
    
    /* ClipRect(&rect); */
    CopyBits((BitMap*)PixelMap,
             (BitMap*)&CanvWin->portBits,
             &rect,&rect,mode,(RgnHandle)0);

    BackColor(blackColor);
    SetPort(savePort);
}


void ClearImage( void )
{
    GrafPtr savePort;
    RGBColor col;
    Rect rect;

    GetPort(&savePort);
    SetPort(CanvWin);
    
    rect.bottom = CanvWin->portRect.bottom-15;
    rect.right = CanvWin->portRect.right-15;
    rect.left = 0;    
    rect.top = 0;
    
    col.red   = BackR<<8 | BackR;
    col.green = BackG<<8 | BackG;
    col.blue  = BackB<<8 | BackB;
    
    MakeRGBPat(BackHand,&col);
    FillCRect(&rect,BackHand);
    SetPort(savePort);
}


static PicHandle CreateMacPicture( void )
{
    GrafPtr savePort;
    RgnHandle saveRgn;
    PicHandle pict;
    Rect rect;

    saveRgn = NewRgn();
    GetPort(&savePort);
    SetPort(CanvWin);
    GetClip( saveRgn );
    
    rect.top = 0;   rect.bottom = YRange;
    rect.left = 0;  rect.right = XRange;
    pict = OpenPicture(&rect);
  
    ForeColor(blackColor);
    BackColor(whiteColor);
   
    ClipRect( &rect );
    CopyBits((BitMap*)PixelMap,
             (BitMap*)&CanvWin->portBits,
             &rect,&rect,srcCopy,(RgnHandle)0);

    ClosePicture();
    
    SetClip(saveRgn);
    BackColor(blackColor);
    DisposeRgn(saveRgn);
    SetPort(savePort);
    return pict;
}


int PrintImage( void )
{
    register int xsize,ysize;
    register int high,wide;
    
    TPrStatus prStatus;
    TPPrPort printPort;
    GrafPtr savePort;
    PicHandle pict;
    short prErr;
    Rect *page;
    Rect rect;

    PrOpen();
    if( !PrintHand )
    {   PrintHand = (THPrint)NewHandle(sizeof(TPrint));
        PrintDefault(PrintHand);
    }

    if( PrJobDialog(PrintHand) )
    {   pict = CreateMacPicture();
        printPort = PrOpenDoc(PrintHand,0,0);
        GetPort(&savePort);    
        SetPort(&printPort->gPort);
        
        PrOpenPage( printPort, 0 );
        page = &(**PrintHand).prInfo.rPage; 
        wide = page->right - page->left;
        high = page->bottom - page->top;
        
        xsize = XRange;
        ysize = YRange;
        if( xsize > wide )
        {   ysize = (int)(((Long)ysize*wide)/xsize);
            xsize = wide;
        }
        if( ysize > high )
        {   xsize = (int)(((Long)xsize*high)/ysize);
            ysize = high;
        }

        rect.top  = page->top  + (high-ysize)>>1;
        rect.left = page->left + (wide-xsize)>>1;
        rect.bottom = rect.top + ysize;
        rect.right = rect.left + xsize;

        DrawPicture( pict, &rect );
        PrClosePage( printPort );

        SetPort(savePort);
        PrCloseDoc(printPort);
        KillPicture(pict);

        if( !(prErr = PrError()) )
            if( (*PrintHand)->prJob.bJDocLoop == bSpoolLoop )
            {   PrPicFile(PrintHand,0,0,0,&prStatus);
                prErr = PrError();
            }

    } else prErr = False;  /* Cancel Print */

    PrClose();
    return !prErr;
}


int ClipboardImage( void )
{
    register long clipErr;
    register long length;
    PicHandle pict;

    ZeroScrap();
    pict = CreateMacPicture();

    HLock((Handle)pict);
    length = (long)GetHandleSize((Handle)pict);
    clipErr = PutScrap(length,'PICT',(Ptr)*pict);
    HUnlock((Handle)pict);

    KillPicture(pict);   
    return !clipErr;
}


void AllocateColourMap( void )
{ 
#ifdef EIGHTBIT
    register PixMapHandle pmHand;
    register GDHandle gdHand;
    register CTabHandle cmap;
    register int i;
    RGBColor col;

    ULut[0] = True;
    ULut[255] = True;

    HLock((Handle)CMap);
    cmap = PixelMap->pmTable;
    (**cmap).ctSeed = GetCTSeed();
    for( i=0; i<256; i++ )
    {   if( ULut[i] )
        {   col.red   = RLut[i]<<8 | RLut[i];
            col.green = GLut[i]<<8 | GLut[i];
            col.blue  = BLut[i]<<8 | BLut[i];
            Lut[i] = i;
        } else col.red = col.green = col.blue = 0;
        (**cmap).ctTable[i].rgb = col;
        SetEntryColor(CMap,i,&col);
    }
    HUnlock((Handle)CMap);

    gdHand = GetMaxDevice(&CanvWin->portRect);
    pmHand = (**gdHand).gdPMap;
    if( (**pmHand).pixelSize >= 8 )
        ActivatePalette(CanvWin);
    FBClear = False;
#endif
}


void UpdateScrollBars( void )
{
    register int pos;

    pos = 50+(int)(50.0*DialValue[1]);
#ifdef USEOLDROUTINENAMES
    SetCtlValue(HScroll,pos);
#else
    SetControlValue(HScroll,pos);
#endif

    pos = 50+(int)(50.0*DialValue[0]);
#ifdef USEOLDROUTINENAMES
    SetCtlValue(VScroll,pos);
#else
    SetControlValue(VScroll,pos);
#endif
}


void SetMouseUpdateStatus( int bool )
{
    MouseUpdateStatus = bool;
}
                         
                         
void SetMouseCaptureStatus( int bool )
{
    MouseCaptureStatus = bool;
}
                         

int LookUpColour( char *name, int *r, int *g, int *b )
{
    UnusedArgument(name);
    UnusedArgument(r);
    UnusedArgument(g);
    UnusedArgument(b);

    return False;
}


void EnableMenus( int flag )
{
    MenuHandle hand;
    long offset;
    
    /* File Menu */
#ifdef USEOLDROUTINENAMES
    hand = GetMHandle(141);
#else
    hand = GetMenuHandle(141);
#endif
    if( flag && !Database )
    {   EnableItem(hand,1);
    } else DisableItem(hand,1);
    
    if( Database && flag )
    {   EnableItem(hand,2);
        EnableItem(hand,3);
        EnableItem(hand,5);
        EnableItem(hand,6);
    } else
    {   DisableItem(hand,2);
        DisableItem(hand,3);
        DisableItem(hand,5);
        DisableItem(hand,6);
    }
    
    /* Edit Menu */
#ifdef USEOLDROUTINENAMES
    hand = GetMHandle(142);
#else
    hand = GetMenuHandle(142);
#endif
    if( Database )
    {   EnableItem(hand,3);
    } else DisableItem(hand,3);
    
    if( flag && (GetScrap(0,'TEXT',&offset)>0) )
    {   EnableItem(hand,4);
    } else DisableItem(hand,4);
    
    if( flag && Database )
    {   EnableItem(hand,7);
    } else DisableItem(hand,7);
    
    /* Middle Menus */
#ifdef USEOLDROUTINENAMES
    if( Database && flag )
    {   EnableItem(GetMHandle(143),0);
        EnableItem(GetMHandle(144),0);
        EnableItem(GetMHandle(145),0);
    } else
    {   DisableItem(GetMHandle(143),0);
        DisableItem(GetMHandle(144),0);
        DisableItem(GetMHandle(145),0);
    }
#else 
    if( Database && flag )
    {   EnableItem(GetMenuHandle(143),0);
        EnableItem(GetMenuHandle(144),0);
        EnableItem(GetMenuHandle(145),0);
    } else
    {   DisableItem(GetMenuHandle(143),0);
        DisableItem(GetMenuHandle(144),0);
        DisableItem(GetMenuHandle(145),0);
    }
#endif
    
    /* Export Menu */
#ifdef USEOLDROUTINENAMES
    hand = GetMHandle(146);
#else 
    hand = GetMenuHandle(146);
#endif
    if( Database )
    {   EnableItem(hand,0);
    } else DisableItem(hand,0);
    DisableMenu = !flag;
    DrawMenuBar();
}


static void DefineMenus( void )
{
    MenuHandle hand;
    register int i;

    /* Apple Menu */
    hand = GetMenu(140);
#ifdef USEOLDROUTINENAMES
    AddResMenu(hand,'DRVR');
#else
    AppendResMenu(hand,'DRVR');
#endif
    InsertMenu(hand,0);

    for( i=141; i<148; i++ )
        InsertMenu( GetMenu(i), 0 );

    /* if( GetEnvirons(smRegionCode) == verUS ) */
    /* SetItem(GetMenuHandle(144),0,'\pColors");   */
    EnableMenus( True );
}


int OpenDisplay( int x, int y )
{
#ifndef EIGHTBIT
    register CTabHandle cmap;
#endif
    register int i;
    Rect rect;
    
    MouseCaptureStatus = False;
    MouseUpdateStatus = False;
    UseHourGlass = True;
    DisableMenu = False;

    for( i=0; i<8; i++ )
        DialValue[i] = 0.0;
    
    ULut[0] = ULut[255] = True;
    RLut[255] = GLut[255] = BLut[255] = 0;
    RLut[0] = GLut[0] = BLut[0] = 255;
    
    XRange = x;   WRange = XRange>>1;
    YRange = y;   HRange = YRange>>1;
    Range = MinFun(XRange,YRange);
    
    CanvWin = GetNewCWindow(150,0,(WindowPtr)-1);
    SetPort(CanvWin);
    
    SizeWindow(CanvWin,x+15,y+15,true);
    ShowWindow(CanvWin);
    
    /* Load Cursors */
    CanvCursor = GetCursor(160);
    CmndCursor = GetCursor(1);
    WaitCursor = GetCursor(4);   
    
    /* Create Scroll Bars */
    rect.left = -1;  rect.right = x+1;  
    rect.top = y;    rect.bottom = y+16;
    HScroll = NewControl(CanvWin,&rect,"\p",true,
                         50,0,100,scrollBarProc,0L);
    
    rect.left = x;   rect.right = x+16;
    rect.top = -1;   rect.bottom = y+1;
    VScroll = NewControl(CanvWin,&rect,"\p",true,
                         50,0,100,scrollBarProc,0L);
    
    DrawGrowIcon(CanvWin);

    /* PixMap! */
    PixelMap = (PixMap*)NewPtr(sizeof(PixMap));
    if( !PixelMap ) return True;

    /* 72.0 DPI Resolution! */
    PixelMap->hRes = 72;
    PixelMap->vRes = 72;
    PixelMap->bounds.left = 0;
    PixelMap->bounds.top = 0;
    PixelMap->cmpSize = 8;
    
    PixelMap->planeBytes = 0;
    PixelMap->pmReserved = 0;
    PixelMap->pmVersion = 0;
    PixelMap->packType = 0;
    PixelMap->packSize = 0;
    
#ifdef EIGHTBIT
    /* Indexed PixMap */
    PixelMap->pixelSize = 8;
    PixelMap->pixelType = 0;
    PixelMap->cmpCount = 1;
    
    PixelMap->pmTable = GetCTable(8);
    (**PixelMap->pmTable).ctSeed = GetCTSeed();
    
    CMap = NewPalette(256,(CTabHandle)0,pmTolerant,0);
    SetPalette(CanvWin,CMap,True);
#else
    /* Direct PixMap */
    PixelMap->pixelSize = 32;
    PixelMap->cmpSize = 8;
    PixelMap->pixelType = RGBDirect;
    PixelMap->cmpCount = 3;
    
    i = sizeof(ColorTable) - sizeof(CSpecArray);
    cmap = (CTabHandle)NewHandle(i);
    (**cmap).ctSeed = 32;
    (**cmap).ctFlags = 0;
    (**cmap).ctSize = 0;
    PixelMap->pmTable = cmap;
#endif
    
    /* Initialise Palette! */
    for( i=1; i<255; i++ )
        ULut[i] = False;
    AllocateColourMap();
    
    PrintHand = (THPrint)NULL;
    BackHand = NewPixPat();
    DefineMenus();
    return False;
}


void CloseDisplay( void )
{
}


void BeginWait( void )
{
    register WindowPtr win;
    
    if( UseHourGlass )
    {   win = FrontWindow();
        if( win==CanvWin || win==CmndWin )
            SetCursor(*WaitCursor);
    }
}



void EndWait( void )
{
    register WindowPtr win;
    GrafPtr savePort;
    Point pos;
    
    /* if( UseHourGlass )? */
    
    win = FrontWindow();
    if( win==CanvWin )
    {   GetPort(&savePort);
        SetPort(CanvWin);
        GetMouse(&pos);
        
        if( (pos.h>0) && (pos.v>0) &&
            (pos.v<CanvWin->portRect.bottom-15) &&
            (pos.h<CanvWin->portRect.right-15) )
        {   SetCursor(*CanvCursor);
        } else ArrowCursor;
        SetPort(savePort);
    } else if( win==CmndWin )
    {   GetPort(&savePort);
        SetPort(CmndWin);
        GetMouse(&pos);
        
        if( (pos.h>0) && (pos.v>0) &&
            (pos.v<CmndWin->portRect.bottom) &&
            (pos.h<CmndWin->portRect.right-15) )
        {   SetCursor(*CmndCursor);
        } else ArrowCursor;
        SetPort(savePort);
    }
}


void SetCanvasTitle( char *ptr )
{
    register char *dst;
    char buffer[255];

    if( ptr )
    {   dst = buffer+1;
        while( *ptr && (dst<buffer+254) )
            *dst++ = *ptr++;
        buffer[0] = (char)((dst-buffer)-1);
        *dst = '\0';
    } else  /* No Title! */
    {   buffer[0] = '\0';
        buffer[1] = '\0';
    }
    SetWTitle(CanvWin,(Byte*)buffer);
}

