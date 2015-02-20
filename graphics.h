/* graphics.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#ifdef APPLEMAC      
#define DefaultWide  400
#define DefaultHigh  400
#endif

#ifdef IBMPC
#define DefaultWide  480
#define DefaultHigh  480
#endif

#ifndef DefaultWide
#define DefaultWide  576
#define DefaultHigh  576
#endif


#ifdef EIGHTBIT
#define LutSize  256
#else
#define LutSize  1024
#endif


#define RFRotateX  0x0001
#define RFRotateY  0x0002
#define RFRotateZ  0x0004
#define RFZoom     0x0008
#define RFTransX   0x0010
#define RFTransY   0x0020
#define RFTransZ   0x0040
#define RFSlab     0x0080
#define RFReSize   0x0100
#define RFColour   0x0200
#define RFRefresh  0x0400

#define RFTrans    0x0070
#define RFRotate   0x0007
#define RFApply    0x017F
#define RFDials    0x00FF
#define RFMagnify  0x0108
#define RFInitial  0x01FF


#ifdef GRAPHICS
double DialValue[8];
int XRange, WRange;
int YRange, HRange;
int UseHourGlass;
int DisableMenu;
int ReDrawFlag;
int Range;

int MouseCaptureStatus;
int MouseUpdateStatus;

Pixel __huge *FBuffer;
short __huge *DBuffer;

unsigned char RLut[LutSize];
unsigned char GLut[LutSize];
unsigned char BLut[LutSize];
unsigned char ULut[LutSize];
Pixel Lut[LutSize];


#ifdef MSWIN
#ifdef EIGHTBIT
LOGPALETTE __far *Palette;
HPALETTE ColourMap;
#endif
HGLOBAL FBufHandle;
HGLOBAL DBufHandle;
HBITMAP PixMap;
HWND CanvWin;
#endif /* MSWIN */

#ifdef APPLEMAC
ControlHandle HScroll;
ControlHandle VScroll;
CursHandle CanvCursor;
CursHandle CmndCursor;
CursHandle WaitCursor;
WindowPtr CanvWin;
WindowPtr CmndWin;
THPrint PrintHand;
Handle FBufHandle;
Handle DBufHandle;
#endif /* APPLEMAC */

#else /* GRAPHICS */
extern double DialValue[8];
extern int XRange, WRange;
extern int YRange, HRange;
extern int UseHourGlass;
extern int DisableMenu;
extern int ReDrawFlag;
extern int Range;

extern int MouseCaptureStatus;
extern int MouseUpdateStatus;

extern Pixel __huge *FBuffer;
extern short __huge *DBuffer;

extern unsigned char RLut[LutSize];
extern unsigned char GLut[LutSize];
extern unsigned char BLut[LutSize];
extern unsigned char ULut[LutSize];
extern Pixel Lut[LutSize];


#ifdef MSWIN
#ifdef EIGHTBIT
extern LOGPALETTE __far *Palette;
extern HPALETTE ColourMap;
#endif
extern HGLOBAL FBufHandle;
extern HGLOBAL DBufHandle;
extern HBITMAP PixMap;
extern HWND CanvWin;
#endif /* MSWIN */

#ifdef APPLEMAC
extern ControlHandle HScroll;
extern ControlHandle VScroll;
extern CursHandle CanvCursor;
extern CursHandle CmndCursor;
extern CursHandle WaitCursor;
extern WindowPtr CanvWin;
extern WindowPtr CmndWin;
extern THPrint PrintHand;
extern Handle FBufHandle;
extern Handle DBufHandle;
#endif /* APPLEMAC */
#endif


int CreateImage( void );
void TransferImage( void );
int ClipboardImage( void );
void ClearImage( void );
int PrintImage( void );

void AllocateColourMap( void );
void UpdateScrollBars( void );
int LookUpColour( char*, int*, int*, int* );
void SetMouseUpdateStatus( int );
void SetMouseCaptureStatus( int );
void SetCanvasTitle( char* );
void EnableMenus( int );
void CloseDisplay( void );
int FetchEvent( int );
void BeginWait( void );
void EndWait( void );

#ifdef MSWIN
int OpenDisplay( HINSTANCE, int );
#else
int OpenDisplay( int, int );
#endif

