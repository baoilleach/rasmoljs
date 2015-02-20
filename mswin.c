/* mswin.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include <windows.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>


#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"

#ifndef SET_BOUNDS
#define SET_BOUNDS   4109
#endif


static BITMAPINFO __far *BitInfo;
static HCURSOR WaitCursor;
static HCURSOR OldCursor;
static HMENU hMenu;


void AllocateColourMap( void )
{
#ifdef EIGHTBIT
    register COLORREF ref;
    register int count;
    register int i;
    
    if( ColourMap )
        DeleteObject(ColourMap);

    count = 0;
    for( i=0; i<256; i++ )
        if( ULut[i] ) 
        {  Palette->palPalEntry[count].peFlags = 0;
           Palette->palPalEntry[count].peRed   = RLut[i];
           Palette->palPalEntry[count].peGreen = GLut[i];
           Palette->palPalEntry[count].peBlue  = BLut[i];

           BitInfo->bmiColors[count].rgbBlue     = BLut[i];
           BitInfo->bmiColors[count].rgbGreen    = GLut[i];
           BitInfo->bmiColors[count].rgbRed      = RLut[i];
           BitInfo->bmiColors[count].rgbReserved = 0;
           count++;
        }   
    Palette->palNumEntries = count;
    BitInfo->bmiHeader.biClrUsed = count;   
    ColourMap = CreatePalette(Palette);

    for( i=0; i<256; i++ )
       if( ULut[i] )
       {   ref = RGB(RLut[i],GLut[i],BLut[i]);
           Lut[i] = GetNearestPaletteIndex(ColourMap,ref);
       }
#else /* EIGHTBIT */
#ifdef THIRTYTWOBIT
    register int i;

    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
            Lut[i] = (RLut[i]<<16) | (GLut[i]<<8) | BLut[i];
#else /* THIRTYTWOBIT */
    register int i;

    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
            Lut[i] = ((RLut[i]&0xf8)<<7)
                   | ((GLut[i]&0xf8)<<2)
                   | ((BLut[i]&0xf8)>>3);
#endif /* THIRTYTWOBIT */
#endif /* EIGHTBIT */
}


int CreateImage( void )
{
    register Long size;

    if( FBufHandle ) GlobalFree(FBufHandle);
    size = (Long)XRange*YRange*sizeof(Pixel)+16;
    FBufHandle = GlobalAlloc(GMEM_MOVEABLE,size);
    return (int)FBufHandle;
}


void TransferImage( void )
{
#ifdef EIGHTBIT
    register HPALETTE OldCMap;
#endif
    register HDC hDC;
        
    if( PixMap )
        DeleteObject(PixMap);

    BitInfo->bmiHeader.biWidth = XRange;
    BitInfo->bmiHeader.biHeight = YRange;

    hDC = GetDC(NULL);
    FBuffer = (Pixel  __huge*)GlobalLock(FBufHandle);
    /* CreateBitMap(XRange,YRange,1,8,FBuffer); */

#ifdef EIGHTBIT
    if( ColourMap )
    {   OldCMap = SelectPalette(hDC,ColourMap,FALSE);
        RealizePalette(hDC);  /* GDI Bug?? */
        PixMap = CreateDIBitmap( hDC, (BITMAPINFOHEADER __far *)BitInfo, 
                                 CBM_INIT, FBuffer, BitInfo, DIB_RGB_COLORS);
        if( OldCMap )                         
            SelectPalette(hDC,OldCMap,False);
    } else
#endif
        PixMap = CreateDIBitmap( hDC, (BITMAPINFOHEADER __far *)BitInfo, 
                                 CBM_INIT, FBuffer, BitInfo, DIB_RGB_COLORS);

    GlobalUnlock(FBufHandle);
    ReleaseDC(NULL,hDC);
    
    InvalidateRect(CanvWin,NULL,FALSE);
    UpdateWindow(CanvWin);
}


void ClearImage( void )
{
    auto RECT rect;
    register HBRUSH hand;
    register HDC hDC;
    
    hDC = GetDC(CanvWin);
    hand = CreateSolidBrush(RGB(RLut[0],GLut[0],BLut[0]));
    GetClientRect(CanvWin,&rect);
    FillRect(hDC,&rect,hand);
    ReleaseDC(CanvWin,hDC);
    DeleteObject(hand);

    if( PixMap )
    {   DeleteObject(PixMap);
        PixMap = NULL;
    }
}


int PrintImage( void )
{
    auto char printer[80];
    register char *device, *driver, *output;
    register int xsize, xres, yres;
    register int dx, dy, caps;

    register HDC hDC;
    auto DOCINFO info;
    auto RECT rect;

    GetProfileString("windows","device", "", printer, 80 );
    if( !(device = strtok(printer,",")) ) return( False );
    if( !(driver = strtok((char*)NULL,", ")) ) return( False );
    if( !(output = strtok((char*)NULL,", ")) ) return( False );

    hDC = CreateDC(driver,device,output,NULL);
    if( !hDC ) return( False );

    caps = GetDeviceCaps( hDC, RASTERCAPS );
    if( !(caps & RC_STRETCHDIB) ) return( False );
    
    xres = GetDeviceCaps( hDC, LOGPIXELSX );
    yres = GetDeviceCaps( hDC, LOGPIXELSY );
    xsize = GetDeviceCaps( hDC, HORZRES );

    dx = xsize - xres;
    dy = (int)(((long)dx*YRange)/XRange);

    /* Should set printer abort procedure */
    /* Position Image on Printed Page */
    rect.top = yres;        rect.bottom = rect.top + dy;
    rect.left = xres>>1;    rect.right = rect.left + dx;
    Escape( hDC, SET_BOUNDS, sizeof(RECT), (char __far*)&rect, NULL );

    /* Start RasWin Document */
    info.cbSize = sizeof(DOCINFO);
    info.lpszDocName = "RasWin";
    info.lpszOutput = NULL;
    StartDoc( hDC, &info );
    StartPage( hDC );

    BitInfo->bmiHeader.biWidth = XRange;
    BitInfo->bmiHeader.biHeight = YRange;
    FBuffer = (Pixel  __huge*)GlobalLock(FBufHandle);

    StretchDIBits( hDC, xres>>1, yres, dx, dy, 
                        0, 0, XRange, YRange, 
                        FBuffer, BitInfo, DIB_RGB_COLORS, SRCCOPY );

    GlobalUnlock(FBufHandle);

    EndPage( hDC );
    EndDoc( hDC );

    DeleteDC( hDC );
    return True;
}


int ClipboardImage( void )
{
#ifdef EIGHTBIT
    register int i;
#endif
    register BITMAPINFO __far *bitmap;
    register char __huge *src;
    register char __huge *dst;
    register long size,len;
    register HANDLE hand;

    if( OpenClipboard(CanvWin) )
    {   EmptyClipboard();

        /* SetClipboardData(CF_DIB,NULL);     */
        /* SetClipboardData(CF_PALETTE,NULL); */

        if( PixMap )
        {   len = (long)XRange*YRange*sizeof(Pixel);
#ifdef EIGHTBIT
            size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
#else
            size = sizeof(BITMAPINFOHEADER);
#endif
            if( (hand=GlobalAlloc(GHND,size+len)) )
            {   bitmap = (BITMAPINFO __far *)GlobalLock(hand);
                bitmap->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
                bitmap->bmiHeader.biWidth = XRange;
                bitmap->bmiHeader.biHeight = YRange;
                bitmap->bmiHeader.biPlanes = 1;
                bitmap->bmiHeader.biBitCount = sizeof(Pixel)<<3;
                bitmap->bmiHeader.biCompression = BI_RGB;
                bitmap->bmiHeader.biSizeImage = len;
                bitmap->bmiHeader.biXPelsPerMeter = 0;
                bitmap->bmiHeader.biYPelsPerMeter = 0;
                bitmap->bmiHeader.biClrImportant = 0;
                bitmap->bmiHeader.biClrUsed = 0;

#ifdef EIGHTBIT
                for( i=0; i<256; i++ )
                    if( ULut[i] )
                    {   bitmap->bmiColors[Lut[i]].rgbBlue  = BLut[i];
                        bitmap->bmiColors[Lut[i]].rgbGreen = GLut[i];
                        bitmap->bmiColors[Lut[i]].rgbRed   = RLut[i];
                    }
#endif

                src = (char __huge*)GlobalLock(FBufHandle);
                dst = ((char __huge*)bitmap)+size;

                /* Transfer the frame buffer */
                while( len-- ) *dst++ = *src++;

                GlobalUnlock(FBufHandle);
                GlobalUnlock(hand);
                SetClipboardData(CF_DIB,hand);
            }
        }

#ifdef EIGHTBIT
        if( ColourMap )
        {   if( (hand = CreatePalette(Palette)) )
                SetClipboardData(CF_PALETTE,hand);
        }
#endif
        CloseClipboard();
        return True;
    } else return False;
}


void SetCanvasTitle( char *ptr )
{
    SetWindowText(CanvWin,ptr);
}


void UpdateScrollBars( void )
{
    register int pos;

    pos = 50-(int)(50.0*DialValue[0]);
    SetScrollPos(CanvWin,SB_VERT,pos,TRUE);
    
    pos = (int)(50.0*DialValue[1])+50;
    SetScrollPos(CanvWin,SB_HORZ,pos,TRUE);
}


void SetMouseUpdateStatus( int stat )
{
    MouseUpdateStatus = stat;
}


void SetMouseCaptureStatus( int stat )
{
    if( stat )
    {   if( !MouseCaptureStatus )
            SetCapture(CanvWin);
    } else
        if( MouseCaptureStatus )
            ReleaseCapture();
    MouseCaptureStatus = stat;
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
    if( flag )
    {   SetMenu(CanvWin,hMenu);
    } else SetMenu(CanvWin,0);
    DisableMenu = !flag;
}


int OpenDisplay( HINSTANCE instance, int mode )
{
    register int i,size;
    auto long style;
    auto RECT rect;

    PixMap = NULL;

    MouseCaptureStatus = False;
    MouseUpdateStatus = False;
    UseHourGlass = True;
    DisableMenu = False;

    for( i=0; i<8; i++ )
         DialValue[i] = 0.0;

    XRange = DefaultWide;   WRange = XRange>>1;
    YRange = DefaultHigh;   HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    rect.top  = 0;   rect.bottom = YRange;
    rect.left = 0;   rect.right  = XRange;

    style = WS_OVERLAPPEDWINDOW | WS_HSCROLL | WS_VSCROLL;

    hMenu = LoadMenu(instance,"RasWinMenu");
    AdjustWindowRect(&rect,style,TRUE);
    CanvWin = CreateWindow("RasWinClass","RasWin Version 2.6.4", style,
                            CW_USEDEFAULT, CW_USEDEFAULT,
                            rect.right-rect.left, 
                            rect.bottom-rect.top,
                            NULL,hMenu,instance,NULL);

    if( !CanvWin ) return False;

#ifdef EIGHTBIT
    size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
    BitInfo = (BITMAPINFO __far*)_fmalloc( size );
    if( !BitInfo ) return False;

    size = sizeof(LOGPALETTE) + 256*sizeof(PALETTEENTRY);
    Palette = (LOGPALETTE __far*)_fmalloc( size );
    if( !Palette ) return False;
    Palette->palVersion = 0x300;

    ColourMap = NULL;
#else
    size = sizeof(BITMAPINFOHEADER);
    BitInfo = (BITMAPINFO __far*)_fmalloc( size );
    if( !BitInfo ) return False;
#endif
    
    WaitCursor = LoadCursor(NULL,IDC_WAIT);

    BitInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    BitInfo->bmiHeader.biBitCount = sizeof(Pixel)<<3;
    BitInfo->bmiHeader.biCompression = BI_RGB;
    BitInfo->bmiHeader.biXPelsPerMeter = 0;
    BitInfo->bmiHeader.biYPelsPerMeter = 0;
    BitInfo->bmiHeader.biClrImportant = 0;
    BitInfo->bmiHeader.biSizeImage = 0;
    BitInfo->bmiHeader.biClrUsed = 0;
    BitInfo->bmiHeader.biPlanes = 1;

    /* Initialise Palette! */
    ULut[0] = True;
    RLut[0] = GLut[0] = BLut[0] = 0;
    for( i=1; i<LutSize; i++ )
        ULut[i] = False;
    AllocateColourMap();

    ShowWindow(CanvWin,mode);
    UpdateScrollBars();
    UpdateWindow(CanvWin);
    return True;
}


void BeginWait( void )
{
    if( UseHourGlass )
        OldCursor = SetCursor(WaitCursor);
}


void EndWait( void )
{
    if( UseHourGlass )
        SetCursor(OldCursor);
}


int FetchEvent( int wait )
{
    /* Avoid Compiler Warning! */
    UnusedArgument(wait);
    return 0;
}


void CloseDisplay( void )
{
#ifdef EIGHTBIT
    if( ColourMap )
        DeleteObject(ColourMap);
#endif
    if( PixMap )
        DeleteObject(PixMap);
}

