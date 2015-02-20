/* raswin.c
 * RasWin Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#include <windows.h>
#include <commdlg.h>
#include <direct.h>

#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#define RASMOL
#include "rasmol.h"
#include "raswin.idm"
#include "molecule.h"
#include "abstree.h"
#include "graphics.h"
#include "pixutils.h"
#include "transfor.h"
#include "cmndline.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "outfile.h"
#include "tmesh.h"

#ifdef SOCKETS
#include <winsock.h>
#define WM_WINSOCK  WM_USER+1
#endif

#ifdef DDEIPC
#include <dde.h>
#endif


/* Microsoft C vs Borland Turbo C */
#ifndef __TURBOC__
#ifndef stricmp
#define stricmp _stricmp
#endif
#define getcwd  _getcwd
char* _getcwd(char*,int);
#endif


#define TwoPi           6.28318531
#define IsIdentChar(x)  ((isalnum(x))||((x)=='_')||((x)=='$'))


/*===========================*/
/* Microsoft Windows DDE IPC */
/*===========================*/

#ifdef DDEIPC
#define DefaultDDETimeOut 10000
#define MaxDDEAdviseNum   32
#define MaxDDEConvNum     8

#define DDEColdLink       0x01
#define DDEWarmLink       0x02
#define DDEHotLink        0x03
#define DDEAckLink        0x04


typedef struct {
        HWND  server;
        HWND  client;
        char  closed;
    } DDEConv;

typedef struct {
        HANDLE data;
        HWND  server;
        HWND  client;
        ATOM  atom;
        char  mode;
        char  item;
        char  wait;
   } DDEAdvise;

static DDEAdvise DDEAdviseData[MaxDDEAdviseNum];
static DDEConv DDEConvData[MaxDDEConvNum];
static int RasWinDDEReady;
static int DDEAdviseCount;
static int DDEConvCount;
static int DDETimeOut;
#endif


#ifdef SOCKETS
/*===============================*/
/* Microsoft Windows Socket IPC  */
/*===============================*/

/* Supported Protocols */
#define ProtoRasMol   0x01

typedef struct {
        SOCKET socket;
        int protocol;
        int advise;
    } IPCConv;

#define MaxIPCConvNum     8
static IPCConv IPCConvData[MaxIPCConvNum];

static SOCKET SocketNo;
static int ServerPort;
static int UseSockets;
#endif


/*=============================*/
/* Generic Advise Information  */
/*=============================*/

#define AMNone         0x00
#define AMPickAtom     0x01
#define AMPickNumber   0x02
#define AMSelectCount  0x04
#define AMMolName      0x08
#define AMPickCoord    0x10

typedef struct {
        int bitmask;
        char *name;
    } AdviseType;

static AdviseType AdviseMap[AdvItemCount] = {
    { AMPickAtom,    "Pick"    },  /* AdvPickAtom    */
    { AMPickNumber,  "PickNo"  },  /* AdvPickNumber  */
    { AMSelectCount, "Count"   },  /* AdvSelectCount */
    { AMMolName,     "Name"    },  /* AdvName        */
    { AMNone,        "Ident"   },  /* AdvIdent       */
    { AMNone,        "Class"   },  /* AdvClass       */
    { AMNone,        "Image"   },  /* AdvImage       */
    { AMPickCoord,   "PickXYZ" }   /* AdvPickCoord   */
        };
 
static char AdviseBuffer[256];
static int AdviseLen;


/*================================*/
/* Command Line Terminal Emulator */
/*================================*/

#define CmndSize   (CmndRows*CmndCols)
#define ScrlMax    80
#define CmndRows   160
#define CmndCols   80


static int CmndStart;
static int ScrlStart;
static int TermCursor;
static int CharWide,CharHigh;
static int TermXPos,TermYPos;
static int TermRows,TermCols;
static char __far *TermScreen;
static HFONT TermFont;
static HWND CmndWin;

static int LabelOptFlag;
static int FileFormat;

static char snamebuf[128];
static char fnamebuf[128];
static char ifilters[512];
static char ofilters[512];
static OPENFILENAME ofn1;
static OPENFILENAME ofn2;
static HINSTANCE hInstance;
static char Text[256];


#ifdef _WIN32
LRESULT CALLBACK MainCallB(HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK CmndCallB(HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK DDECallB( HWND, UINT, WPARAM, LPARAM);
BOOL CALLBACK AboutCallB(  HWND, UINT, WPARAM, LPARAM);
BOOL CALLBACK InfoCallB(   HWND, UINT, WPARAM, LPARAM);
#else
LONG FAR PASCAL MainCallB( HWND, UINT, WPARAM, LPARAM);
LONG FAR PASCAL CmndCallB( HWND, UINT, WPARAM, LPARAM);
LONG FAR PASCAL DDECallB(  HWND, UINT, WPARAM, LPARAM); 
BOOL FAR PASCAL AboutCallB(HWND, UINT, WPARAM, LPARAM);
BOOL FAR PASCAL InfoCallB( HWND, UINT, WPARAM, LPARAM);
#endif


#ifdef DDEIPC
static void CloseDDELinks( void )
{
    register long alarm;
    register int i;
    auto MSG message;
    
    for( i=0; i<MaxDDEConvNum; i++ )
        if( DDEConvData[i].server )
        {   DDEConvData[i].closed = True;
            PostMessage( DDEConvData[i].client, WM_DDE_TERMINATE,
                         (WPARAM)DDEConvData[i].server, (LPARAM)0 );
        } 
   
    alarm = GetTickCount() + DDETimeOut;
    while( PeekMessage(&message,NULL,WM_DDE_FIRST,WM_DDE_LAST,PM_REMOVE) )
    {   DispatchMessage( &message );
        if( message.message == WM_DDE_TERMINATE )
            if( !DDEConvCount ) break;
            
        /* Terminate Time Out */
        if( (long)GetTickCount() > alarm )
        {   for( i=0; i<MaxDDEConvNum; i++ )
                if( DDEConvData[i].server )
                    DestroyWindow( DDEConvData[i].server );
            break;
        }
    }

    for( i=0; i<MaxDDEAdviseNum; i++ )
        if( DDEAdviseData[i].server && DDEAdviseData[i].wait )
        {   GlobalDeleteAtom(DDEAdviseData[i].atom);
            GlobalFree(DDEAdviseData[i].data);
            DDEAdviseData[i].wait = False;
            break;
        }
}
#endif


#ifdef SOCKETS
static void CloseSockets( void )
{
    register int i;

    if( UseSockets )
    {   if( SocketNo != INVALID_SOCKET )
            closesocket(SocketNo);

        for( i=0; i<MaxIPCConvNum; i++ )
            if( IPCConvData[i].protocol )
            {   closesocket(IPCConvData[i].socket);
                IPCConvData[i].protocol = 0;
            }

        SocketNo = INVALID_SOCKET;
        UseSockets = False;
        WSACleanup();
    }
}
#endif


void RasMolExit( void )
{
    DeleteObject(TermFont);
#ifdef DDEIPC
    CloseDDELinks();
#endif
#ifdef SOCKETS
    CloseSockets();
#endif
    CloseDisplay();
    exit(0);
}


void RasMolFatalExit( char *msg )
{
    MessageBox(NULL,msg,"RasMol Fatal Error!",
        MB_OK | MB_ICONEXCLAMATION | MB_APPLMODAL );
 
    /* PostQuitMessage(0); */
    DeleteObject(TermFont);
#ifdef DDEIPC
    CloseDDELinks();
#endif
#ifdef SOCKETS
    CloseSockets();
#endif
    CloseDisplay();
    exit(1);    
}


static void LoadInitFile( void )
{
    register char *src,*dst;
    register FILE *initrc;
    register char *fname;

    fname = "RASMOL.INI";
    initrc = fopen(fname,"rb");
    if( !initrc && (src=(char*)getenv("HOME")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '\\';

        src = fname; fname = fnamebuf;
        while( *src )
            *dst++ = *src++;
        *dst = '\0';
        initrc = fopen(fname,"rb");
    }

    if( initrc )
        LoadScriptFile(initrc,fname);
}


static void SetTermScroll( int pos )
{
    SetScrollPos(CmndWin,SB_VERT,pos,True);
    InvalidateRect(CmndWin,NULL,True);
    ScrlStart = ScrlMax - pos;
}


void WriteChar( int ch )
{
    register int i;
    auto RECT rect;

    /* Scroll to bottom! */
    if( ScrlStart )
        SetTermScroll( ScrlMax );
        
    switch( ch )
    {    case(0x07):  MessageBeep(0);
                      break;
                      
         case(0x08):  if( TermXPos>0 )
                      {   TermXPos--;
                          if( TermCursor )
                              SetCaretPos(TermXPos*CharWide,
                                          TermYPos*CharHigh);
                      }
                      break;
                      
         case(0x0D):  if( TermXPos )
                      {   if( TermCursor )
                              SetCaretPos(0,TermYPos*CharHigh);
                          TermXPos=0;
                      }
                      break;
                      
         case(0x0A):  if( TermYPos==TermRows-1 )
                      {   CmndStart++;
                          if( CmndStart == CmndRows )
                              CmndStart = 0;
                              
                          i = TermYPos + CmndStart;
                          if( i >= CmndRows ) i -= CmndRows;
                          _fmemset(TermScreen+i*CmndCols,' ',CmndCols);
                          InvalidateRect(CmndWin,NULL,FALSE);
                      } else TermYPos++;
                      TermXPos = 0;

                      if( TermCursor )
                          SetCaretPos(0,TermYPos*CharHigh);
                      UpdateWindow(CmndWin);
                      break;
                      
         
         default:     i = TermYPos + CmndStart;
                      if( i >= CmndRows ) i -= CmndRows;
                      TermScreen[i*CmndCols+TermXPos]=ch;
                      if( TermXPos < TermCols )
                      {   rect.top = TermYPos*CharHigh; 
                          rect.left = TermXPos*CharWide;
                          rect.bottom = rect.top+CharHigh;
                          rect.right = rect.left+CharWide;
                          InvalidateRect(CmndWin,&rect,FALSE);
                      }
                      
                      if( TermXPos==CmndCols-1 )
                      {   if( TermYPos==TermRows-1 )
                          {   CmndStart++;
                              if( CmndStart == CmndRows )
                                  CmndStart = 0;
                                  
                              i = TermYPos + CmndStart;
                              if( i >= CmndRows ) i -= CmndRows;
                              _fmemset(TermScreen+i*CmndCols,' ',CmndCols);
                              InvalidateRect(CmndWin,NULL,FALSE);
                          } else TermYPos++;
                          TermXPos=0;
                      } else TermXPos++;

                      if( TermCursor )
                          SetCaretPos(TermXPos*CharWide,
                                      TermYPos*CharHigh);
                      break;
                      
    }
}


void WriteString( char *ptr )
{
    while( *ptr )
        WriteChar(*ptr++);
}


static int InitTerminal( HINSTANCE instance )
{
    register HDC hDC;
    register DWORD style;
    auto TEXTMETRIC Text;
    auto LOGFONT LogFont;
    auto RECT rect;

    LogFont.lfHeight     = 12;
    LogFont.lfWidth      = 8;
    LogFont.lfEscapement = 0;
    LogFont.lfWeight     = 0;
    LogFont.lfItalic     = 0;
    LogFont.lfUnderline  = 0;
    LogFont.lfStrikeOut  = 0;
    
    LogFont.lfCharSet        = OEM_CHARSET;
    LogFont.lfOutPrecision   = OUT_DEFAULT_PRECIS;
    LogFont.lfClipPrecision  = CLIP_DEFAULT_PRECIS;
    LogFont.lfQuality        = DEFAULT_QUALITY;
    LogFont.lfPitchAndFamily = FIXED_PITCH | FF_MODERN;
    LogFont.lfFaceName[0]    = '\0';
    TermFont = CreateFontIndirect(&LogFont);

    /* TermFont = GetStockObject(ANSI_FIXED_FONT); */

    /* Default Window Size */
    TermCols = 80;  TermRows = 24;
    ScrlStart = CmndStart = 0;
    TermXPos = TermYPos = 0;

    TermScreen = (char __far*)_fmalloc(CmndSize*sizeof(char));
    if( !TermScreen ) return( False );
    _fmemset(TermScreen,' ',CmndSize);
    TermCursor = False;

    hDC = GetDC(NULL);
    SelectObject(hDC,TermFont);
    GetTextMetrics(hDC,&Text);  
    ReleaseDC(NULL,hDC);
    
    CharWide = Text.tmAveCharWidth;
    CharHigh = Text.tmHeight + Text.tmExternalLeading;

    rect.top  = 0;   rect.bottom = TermRows*CharHigh;
    rect.left = 0;   rect.right  = TermCols*CharWide;
    style = WS_OVERLAPPEDWINDOW | WS_VSCROLL;
    
    AdjustWindowRect(&rect,style,False);
#ifdef WS_EX_CLIENTEDGE
    CmndWin = CreateWindowEx(WS_EX_CLIENTEDGE,
                           "RasCliClass", "RasMol Command Line",
                           style, CW_USEDEFAULT, CW_USEDEFAULT,
                           rect.right-rect.left, rect.bottom-rect.top,
                           (HWND)0, (HMENU)0, instance, (LPVOID)0 );
#else
    CmndWin = CreateWindow("RasCliClass", "RasMol Command Line",
                           style, CW_USEDEFAULT, CW_USEDEFAULT,
                           rect.right-rect.left, rect.bottom-rect.top,
                           (HWND)0, (HMENU)0, instance, (LPVOID)0 );
#endif

    if( !CmndWin ) return( False );

    SetScrollRange(CmndWin,SB_VERT,0,ScrlMax,FALSE); 
    SetScrollPos(CmndWin,SB_VERT,ScrlMax,FALSE);
    ShowWindow(CmndWin,SW_SHOWMINNOACTIVE);
    return True;
}


static void PaintScreen( void )
{
    register int SRow,ERow,SCol,ECol;
    register char __far *ptr;
    register int row,len;
    register int x,y;
    
    auto RECT rect;
    auto PAINTSTRUCT ps;
    register HFONT font;
    register HDC hDC;
    
    hDC = BeginPaint(CmndWin,&ps);
    font = (HFONT)SelectObject(hDC,TermFont);
    SetBkColor(hDC,GetSysColor(COLOR_WINDOW));
    SetTextColor(hDC,RGB(0,0,0));
    SetBkMode(hDC,OPAQUE);
    
    SRow = ps.rcPaint.top/CharHigh;
    if( SRow >= TermRows )
    {   SRow = TermRows-1;
    } else if( SRow < 0 )
        SRow = 0;

    ERow = ps.rcPaint.bottom/CharHigh;
    if( ERow >= TermRows )
    {   ERow = TermRows-1;
    } else if( ERow < 0 ) 
        ERow = 0;
    
    SCol = ps.rcPaint.left/CharWide;
    if( SCol >= TermCols )
    {   SCol = TermCols-1;
    } else if( SCol < 0 ) 
        SCol = 0;
    
    ECol = ps.rcPaint.right/CharWide;
    if( ECol >= TermCols )
    {   ECol = TermCols-1;
    } else if( ECol < 0 ) 
        ECol = 0;

    len = ECol-SCol+1;
    x = SCol*CharWide;
    y = SRow*CharHigh;

    rect.right = x+len*CharWide;
    rect.left = x;   

    SRow += CmndStart - ScrlStart;
    if( SRow >= CmndRows )
    {   SRow -= CmndRows;
    } else if( SRow < 0 )
        SRow += CmndRows;
        
    ERow += CmndStart - ScrlStart;
    if( ERow >= CmndRows )
    {   ERow -= CmndRows;
    } else if( ERow < 0 )
        ERow += CmndRows;
        
    row = SRow;
    ptr = TermScreen + CmndCols*row + SCol;
    while( True )
    {   rect.top = y;    
        rect.bottom = y+CharHigh;
        ExtTextOut(hDC,x,y,ETO_OPAQUE|ETO_CLIPPED,
                   &rect,ptr,len,NULL);
                   
        if( row != ERow )
        {   ptr += CmndCols;
            row++;
            if( row == CmndRows )
            {   ptr = TermScreen + SCol;
                row = 0;
            }
        } else break;
        y += CharHigh;
    }
    
    SelectObject(hDC,font);
    EndPaint(CmndWin,&ps);
    
    if( TermCursor )
    {   row = TermYPos + ScrlStart;
        if( row < TermRows )
        {   SetCaretPos(TermXPos*CharWide,row*CharHigh);
            ShowCaret(CmndWin);
        } else HideCaret(CmndWin);
    }    
}


static void DetermineHostInfo( void )
{
#ifdef _WIN32
    auto OSVERSIONINFO osinfo;
    auto SYSTEM_INFO sysinfo;
    register char *winver;
    register char *cpu;
    register int count;

    osinfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    GetVersionEx(&osinfo);
    switch( osinfo.dwPlatformId )
    {   case 0:  /* VER_PLATFORM_WIN32s */
                 winver = "Win32s";
                 break;

        case 1:  /* VER_PLATFORM_WIN32_WINDOWS */
                 if( (osinfo.dwMajorVersion==4) &&
                     (osinfo.dwMinorVersion==0) ) 
                 {   winver = "Windows95";     /* 4.00 */
                 } else if( (osinfo.dwMajorVersion==4) &&
                            (osinfo.dwMinorVersion<90) )
                 {   winver = "Windows98";     /* 4.10 */
                 } else winver = "WindowsME";  /* 4.90 */
                 break;

        case 2:  /* VER_PLATFORM_WIN32_NT */
                 if( osinfo.dwMajorVersion < 5 )
                 {   winver = "WindowsNT";
                 } else winver = "Windows2000";
                 break;

        case 3:  /* VER_PLATFORM_WIN32_CE */
                 winver = "WindowsCE";
                 break;

        default: winver = "Win32";
    }

    if( osinfo.dwPlatformId == VER_PLATFORM_WIN32_NT )
    {   GetSystemInfo(&sysinfo);
        /* sysinfo.wProcessorArchitecture */
        switch( *((WORD*)&sysinfo) )
        {   case PROCESSOR_ARCHITECTURE_INTEL:
                switch( sysinfo.wProcessorLevel )
                {   case 3:  cpu = "Intel 386";         break;
                    case 4:  cpu = "Intel 486";         break;
                    case 5:  cpu = "Intel Pentium";     break;
                    case 6:  cpu = "Intel Pentium III"; break;
                    default: cpu = "Unknown Intel";     break;
                }
                break;

            case PROCESSOR_ARCHITECTURE_MIPS:
                switch( sysinfo.wProcessorLevel )
                {   case  1: cpu = "MIPS R2000";    break;
                    case  2: cpu = "MIPS R3000";    break;
                    case  3: cpu = "MIPS R6000";    break;
                    case  4: cpu = "MIPS R4000";    break;
                    case  6: cpu = "MIPS R6000A";   break;
                    case  9: cpu = "MIPS R10000";   break;
                    case 10: cpu = "MIPS R4200";    break;
                    case 11: cpu = "MIPS R4300";    break;
                    case 16: cpu = "MIPS R8000";    break;
                    case 32: cpu = "MIPS R4600";    break;
                    case 33: cpu = "MIPS R4700";    break;
                    case 34: cpu = "MIPS R4650";    break;
                    case 35: cpu = "MIPS R5000";    break;
                    default: cpu = "Unknown MIPS";  break;
                }
                break;

            case PROCESSOR_ARCHITECTURE_ALPHA:
                switch( sysinfo.wProcessorLevel )
                {   case 21064:  cpu = "Alpha 21064";  break;
                    case 21066:  cpu = "Alpha 21066";  break;
                    case 21164:  cpu = "Alpha 21164";  break;
                    default: cpu = "Unknown Alpha";    break;
                }
                break;

            case PROCESSOR_ARCHITECTURE_PPC:
                switch( sysinfo.wProcessorLevel )
                {   case 1:  cpu = "PPC 601";     break;
                    case 3:  cpu = "PPC 603";     break;
                    case 4:  cpu = "PPC 604";     break;
                    case 6:  cpu = "PPC 603+";    break;
                    case 9:  cpu = "PPC 604+";    break;
                    case 20: cpu = "PPC 620";     break;
                    default: cpu = "Unknown PPC"; break;
                }
                break;

            default:
                cpu = "unrecognised";
                break;
        }
    } else /* Windows 95 or Windows98 */
    {   GetSystemInfo(&sysinfo);
        switch( sysinfo.dwProcessorType )
        {   case(386):   cpu = "Intel 386";          break;
            case(486):   cpu = "Intel 486";          break;
            case(586):   cpu = "Intel Pentium";      break;
            case(686):   cpu = "Intel Pentium III";  break;
            case(860):   cpu = "Intel i860";         break;
            case(2000):  cpu = "MIPS R2000";         break;
            case(3000):  cpu = "MIPS R3000";         break;
            case(4000):  cpu = "MIPS R4000";         break;
            case(4400):  cpu = "MIPS R4400";         break;
            case(4600):  cpu = "MIPS R4600";         break;
            case(5000):  cpu = "MIPS R5000";         break;
            case(8000):  cpu = "MIPS R8000";         break;
            case(10000): cpu = "MIPS R10000";        break;
            case(21064): cpu = "DEC Alpha 21064";    break;
            case(21066): cpu = "DEC Alpha 21066";    break;
            case(21164): cpu = "DEC Alpha 21164";    break;
            default:     cpu = "unrecognised";
        }
    }

    count = sysinfo.dwNumberOfProcessors;
    if( count > 1 )
    {   sprintf(Text,"%s on %d %s CPUs\n",winver,count,cpu);
    } else sprintf(Text,"%s on a single %s CPU",winver,cpu);
#else
    register DWORD flags;
    register char *cpu;

    flags = GetWinFlags();
    if( flags & WF_CPU286 )
    {      cpu = "286";
    } else if( flags & WF_CPU386 )
    {      cpu = "386";
    } else cpu = "486";
                              
    if( !(flags&WF_80x87) )
    {   sprintf(Text,"%s without maths coprocessor",cpu);
    } else sprintf(Text,"%s with maths coprocessor",cpu); 
#endif
}


#ifdef _WIN32
BOOL CALLBACK AboutCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#else
BOOL FAR PASCAL AboutCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#endif
{
    UnusedArgument(lArg);

    switch(uMsg)
    {   case(WM_INITDIALOG):  
                DetermineHostInfo();
                SetDlgItemText(hWin,IDD_HARDWARE,Text);         
                return TRUE;
    
        case(WM_COMMAND):     
#ifdef _WIN32
                if( LOWORD(wArg) == IDOK )
#else
                if( wArg == IDOK )
#endif
                {   EndDialog(hWin,TRUE);
                    return TRUE;
                }
                break;
    }
    return 0;
}


static void DisplayMoleculeInfo( HWND hWin )
{
    register int line;
    register int len;

    line = IDD_INFOTEXT;
    
    if( *Info.moleculename )
    {   sprintf(Text," Molecule name ...... %s",Info.moleculename);
        SetDlgItemText(hWin,line++,Text);
    }
    
    if( *Info.classification )
    {   sprintf(Text," Classification ......... %s",Info.classification);
        SetDlgItemText(hWin,line++,Text);
    }
    
    if( *Info.identcode )
    {   sprintf(Text," Brookhaven code .. %s",Info.identcode);
        SetDlgItemText(hWin,line++,Text);
    }
    
    if( Info.chaincount>1 )
    {   sprintf(Text," Number of chains .. %d",Info.chaincount);
        SetDlgItemText(hWin,line++,Text);
    }
    
    len = sprintf(Text," Number of groups .. %d",MainGroupCount);
    if( HetaGroupCount ) sprintf(Text+len," (%d)",HetaGroupCount);
    SetDlgItemText(hWin,line++,Text);

    len = sprintf(Text," Number of atoms ... %ld",MainAtomCount);
    if( HetaAtomCount ) sprintf(Text+len," (%ld)",HetaAtomCount);
    SetDlgItemText(hWin,line++,Text);

    sprintf(Text," Number of bonds ... %ld",Info.bondcount);
    SetDlgItemText(hWin,line++,Text);
}


#ifdef _WIN32
BOOL CALLBACK InfoCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#else
BOOL FAR PASCAL InfoCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#endif
{
    UnusedArgument(lArg);

    switch(uMsg)
    {   case(WM_INITDIALOG):  
                DisplayMoleculeInfo(hWin);
                return TRUE;
                              
        case(WM_COMMAND):
#ifdef _WIN32
                if( LOWORD(wArg) == IDOK )
#else
                if( wArg == IDOK )
#endif
                {   EndDialog(hWin,TRUE);
                    return TRUE;
                }
                break;
    }
    return 0;
}


#ifdef DDEIPC
static char *GetItemName( int item )
{
    switch( item )
    {   case  0:  return (char*)0;
        case -1:  return "Topics";
        case -2:  return "SysItems";
        case -3:  return "Formats";
        case -4:  return "Status";
        case -5:  return "Items";
    }

    if( item <= AdvItemCount )
        return AdviseMap[item-1].name;
    return "";
}
#endif


static HANDLE RenderClipboard( WPARAM format )
{
#ifdef EIGHTBIT
    register int i;
#endif
    register BITMAPINFO __far *bitmap;
    register char __huge *src;
    register char __huge *dst;
    register HANDLE result;
    register long size,len;
   
    if( format==CF_PALETTE )
    {
#ifdef EIGHTBIT
        if( ColourMap )
            return CreatePalette(Palette);
#endif
        return NULL;
    }
    
    if( !PixMap || (format!=CF_DIB) )
        return NULL;

    len = (long)XRange*YRange*sizeof(Pixel);
#ifdef EIGHTBIT
    size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
#else
    size = sizeof(BITMAPINFOHEADER);
#endif
    if( !(result=GlobalAlloc(GHND,size+len)) ) return NULL;
    
    bitmap = (BITMAPINFO __far *)GlobalLock(result);
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
    GlobalUnlock(result);
    return result;    
}


#ifdef DDEIPC
static void SendItemData( HWND hSrc, HWND hDst,
                          int mode, int item, int advise )
{
    register DDEDATA FAR *data;
    register HANDLE FAR *hImage;
    register HANDLE hData;
    register LPARAM lArg;
    register ATOM atom;

    register char __far *dest;
    register char *src, *dst;
    register char *name;
    register Long len;
    register int i;

    name = GetItemName(item);

    if( mode == DDEWarmLink )
    {   atom = GlobalAddAtom(name);
#ifdef _WIN32
        lArg = PackDDElParam(WM_DDE_DATA,0,atom);
#else
        lArg = MAKELONG(0,atom);
#endif
        if( !PostMessage(hDst,WM_DDE_DATA,(WPARAM)hSrc,lArg) )
            GlobalDeleteAtom(atom);
        return;
    }

    dst = AdviseBuffer;
    switch( (item<0)? item : item-1 )
    {   case(-1): /* Topics */
                  src="System\tRemoteControl"; 
                  while( *dst++ = *src++ ); break;

        case(-2): /* SysItems */
                  src = "Topics\tSysItems\tFormats\tStatus\tItems";
                  while( *dst++ = *src++ ); break;

        case(-3): /* Formats */
#ifdef EIGHTBIT
                  src = "DIB\tTEXT\tPalette\tLink";
#else
                  src = "DIB\tTEXT\tLink";
#endif
                  while( *dst++ = *src++ ); break;

        case(-4): /* Status */
                  src = RasWinDDEReady? "Ready" : "Busy";
                  while( *dst++ = *src++ ); break;

        case(-5): /* Items */
                  for( i=0; i<AdvItemCount; i++ )
                  {   if( i ) *dst++ = '\t';
                      src = AdviseMap[i].name; 
                      while( *src )
                          *dst++ = *src++;
                  }
                  *dst = '\0';
                  break;

        case(AdvPickAtom):
                  if( QAtom )
                  {   src = Residue[QGroup->refno];
                      if( src[0]!=' ' ) *dst++ = src[0];
                      *dst++  = src[1]; *dst++ = src[2];

                      sprintf(dst,"%d",QGroup->serno);
                      for( dst=Text; *dst; dst++ );
                      if( QChain->ident!=' ' )
                      {   *dst++ = ':';
                          *dst++ = QChain->ident;
                      }
                      *dst++ = '.';
                      
                      src = ElemDesc[QAtom->refno];
                      if( src[0]!=' ' ) *dst++ = src[0];
                      *dst++  = src[1]; *dst++ = src[2];
                      if( src[3]!=' ' ) *dst++ = src[3];
                  } 
                  *dst = '\0';
                  break;

        case(AdvPickNumber):
                  if( QAtom )
                  { sprintf(dst,"%d",QAtom->serno);
                  } else *dst = '\0';
                  break;

        case(AdvSelectCount):
                  sprintf(dst,"%ld",SelectCount);
                  break;
                  
        case(AdvName):
                  src = Info.moleculename;
                  while( *dst++ = *src++ );
                  break;

        case(AdvPickCoord):
                  if( QAtom )
                  { sprintf( dst, "%ld\t%ld\t%ld",
                             QAtom->xorg, QAtom->yorg, QAtom->zorg);
                  } else *dst = '\0';
                  break;

        default:  *dst = '\0';
                  break;
    }

    len = sizeof(DDEDATA);
    if( item == AdvImage )
    {   len += sizeof(HANDLE);
        AdviseLen = 0;
    } else
    {   AdviseLen = strlen(AdviseBuffer)+1;
        len += AdviseLen;
    }

    if( hData = GlobalAlloc(GHND|GMEM_DDESHARE,len) )
    {   if( data = (DDEDATA FAR*)GlobalLock(hData) )
        {   data->fResponse = (mode!=DDEAckLink);
            data->fAckReq = (mode==DDEColdLink);
            data->fRelease = True;

            if( item == AdvImage )
            {   data->cfFormat = CF_DIB;
                hImage = (HANDLE __far*)&data->Value[0];
                *hImage = RenderClipboard(CF_DIB);

            } else 
            {   data->cfFormat = CF_TEXT;
                dest = (char __far*)&data->Value[0];
                memcpy(dest,AdviseBuffer,len+1);
                /* Correctly terminate the data string */
                /* dest[len]='\r'; dest[len+1]='\n';   */
            }
            
            GlobalUnlock(hData);
            atom = GlobalAddAtom(name);
#ifdef _WIN32
            lArg = PackDDElParam(WM_DDE_DATA,(UINT)hData,(UINT)atom);
#else
            lArg = MAKELONG((UINT)hData,(UINT)atom);
#endif
            if( PostMessage(hDst,WM_DDE_DATA,(WPARAM)hSrc,lArg) )
            {   if( mode == DDEAckLink )
                {   SetTimer( hSrc, (UINT)hDst, DDETimeOut, NULL );
                    DDEAdviseData[advise].data = hData;
                    DDEAdviseData[advise].atom = atom;
                    DDEAdviseData[advise].wait = True;
                }
                return;
            }
            GlobalDeleteAtom(atom);
        }
        GlobalFree(hData);
    }
}
#endif


#ifdef DDEIPC
static int GetItemNumber( ATOM atom )
{
    register int i;

    GlobalGetAtomName(atom,Text,240);

    for( i=1; i<6; i++ )
        if( !stricmp(Text,GetItemName(-i)) )
            return -i;

    for( i=0; i<AdvItemCount; i++ )
        if( !stricmp(Text,AdviseMap[i].name) )
            return i+1;
    return 0;
}
#endif


#ifdef SOCKETS
static void PrepareIPCAdviseItem( int item )
{
    register char *src, *dst;
    register int i,flag;

    dst = AdviseBuffer;
    src = AdviseMap[item].name;
    while( *src ) *dst++ = *src++;
    *dst++ = ':';  *dst++ = ' ';

    switch( item )
    {   case(AdvPickAtom):
                  if( QAtom )
                  {   src = Residue[QGroup->refno];
                      flag = False;
                      for( i=0; i<3; i++ )
                          if( (src[i]!=' ') && !isalpha(src[i]) )
                              flag = True;

                      if( flag ) *dst++ = '[';
                      for( i=0; i<3; i++ )
                          if( src[i]!=' ' )
                              *dst++ = src[i];
                      if( flag ) *dst++ = ']';
                      sprintf(dst,"%d",QGroup->serno);
                      for( dst=AdviseBuffer; *dst; dst++ );
                      if( QChain->ident!=' ' )
                      {   if( isdigit(QChain->ident) ) *dst++ = ':';
                          *dst++ = QChain->ident;
                      }
                      *dst++ = '.';

                      src = ElemDesc[QAtom->refno];
                      for( i=0; i<4; i++ )
                          if( src[i]!=' ' )
                              *dst++ = src[i];
                  }
                  *dst++ = '\n';
                  *dst = '\0';
                  break;

        case(AdvPickNumber):
                  if( !QAtom )
                  {   *dst++ = '\n'; *dst = '\0';
                  } else sprintf(dst,"%ld\n",(long)QAtom->serno);
                  break;

        case(AdvSelectCount):
                  sprintf(dst,"%ld\n",(long)SelectCount);
                  break;

        case(AdvName):
                  src = Info.moleculename;
                  while( *src ) *dst++ = *src++;
                  *dst++ = '\n'; *dst = '\0';
                  break;

        case(AdvPickCoord):
                  if( !QAtom )
                  {   *dst++ = '\n'; *dst = '\0';
                  } else sprintf( dst, "%ld\t%ld\t%ld\n",
                             (long)QAtom->xorg,
                             (long)QAtom->yorg,
                             (long)QAtom->zorg);
                  break;

        default:  *dst++ = '\n';
                  *dst = '\0';
                  break;
    }
    AdviseLen = strlen(AdviseBuffer)+1;
}
#endif


#ifdef SOCKETS
static void AdviseSocketUpdate( int item )
{
    register int mask;
    register int i;

    if( (item>=0) && UseSockets )
    {   mask = AdviseMap[item].bitmask;
        if( mask )
        {   AdviseLen = 0;
            for( i=0; i<MaxIPCConvNum; i++ )
                if( IPCConvData[i].protocol && (IPCConvData[i].advise&mask) )
                {   if( !AdviseLen ) PrepareIPCAdviseItem(item);
                    send(IPCConvData[i].socket,AdviseBuffer,AdviseLen,0);
                }
        }
    }
}
#endif


#ifdef DDEIPC
static void AdviseDDEUpdate( int item )
{
    register DDEAdvise *ptr;
    register int i;

    if( DDEAdviseCount )
    {   if( item >= 0 ) item++;
        ptr = DDEAdviseData;
        for( i=0; i<MaxDDEAdviseNum; i++ )
        {   if( ptr->server && (ptr->item==(char)item) && !ptr->wait )
                SendItemData(ptr->server,ptr->client,ptr->mode,item,i);
            ptr++;
        }
    }
}
#endif


void AdviseUpdate( int item )
{
    /* Avoid compiler warning! */
    UnusedArgument(item);

#ifdef SOCKETS
    AdviseSocketUpdate(item);
#endif
#ifdef DDEIPC
    AdviseDDEUpdate(item);
#endif
}


void RefreshScreen( void )
{
    ReDrawFlag &= ~RFTransZ;

    if( ReDrawFlag )
    {
#ifdef DDEIPC
        if( RasWinDDEReady )
        {   RasWinDDEReady = False;
            AdviseUpdate( -4 );
        }
#endif
        
        if( ReDrawFlag & RFReSize )
            ReSizeScreen();

        if( ReDrawFlag & RFColour )
        {   ClearImage();
            DefineColourMap();
        }

        if( Database )
        {   BeginWait();
            if( ReDrawFlag & RFApply ) 
                ApplyTransform();
            DrawFrame();
            TransferImage();
            EndWait();
        } else
        {   ClearBuffers();
            TransferImage();
        }
        ReDrawFlag = 0;
    }
}


#ifdef DDEIPC
static void GenerateDDEReply( HWND hDst, HWND hSrc, int flag, UINT data )
{
    register LPARAM lArg;

#ifdef _WIN32
    lArg = PackDDElParam(WM_DDE_ACK,(flag?0x8000:0),data);
#else
    lArg = MAKELPARAM((flag?0x8000:0),data);
#endif
    PostMessage( (HWND)hDst, WM_DDE_ACK, (WPARAM)hSrc, lArg );
}
#endif


#ifdef DDEIPC
#ifdef _WIN32
LRESULT CALLBACK DDECallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#else
LONG FAR PASCAL DDECallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#endif
{
    register DDEADVISE FAR *options;
    register HWND hDest;
    auto UINT loword;
    auto UINT hiword;

    register DDEConv *ptr;
    register char __huge *cmnd;
    register int item, stat;
    register int format,i;
    register int flag;
    
    switch( uMsg )
    {   case( WM_TIMER ):    
                    KillTimer( hWin, wArg );
                    for( i=0; i<MaxDDEAdviseNum; i++ )
                        if( (DDEAdviseData[i].server==hWin) &&
                             DDEAdviseData[i].wait )
                        {   GlobalDeleteAtom(DDEAdviseData[i].atom);
                            GlobalFree(DDEAdviseData[i].data);
                            DDEAdviseData[i].wait = False;
                            break;
                        }
                    return 0L;

        case( WM_DDE_ACK ):
                    KillTimer( hWin, wArg );
#ifdef _WIN32
                    if( UnpackDDElParam(uMsg,lArg,&loword,&hiword) )
                    {   FreeDDElParam(uMsg,lArg);
                    } else return 0L;
#else
                    loword = LOWORD(lArg);
                    hiword = HIWORD(lArg);
#endif

                    if( hiword )
                    {   item = GetItemNumber( (ATOM)hiword );
                        GlobalDeleteAtom( (ATOM)hiword );
                    } else item = 0;

                    if( !item )
                        return 0L;

                    if( !(loword&0x8000) )
                    {   for( i=0; i<MaxDDEAdviseNum; i++ )
                            if( (DDEAdviseData[i].server==hWin) &&
                                (DDEAdviseData[i].item==(char)item) &&
                                 DDEAdviseData[i].wait )
                            {   if( DDEAdviseData[i].atom != (ATOM)hiword )
                                    GlobalDeleteAtom(DDEAdviseData[i].atom);
                                GlobalFree(DDEAdviseData[i].data);
                                DDEAdviseData[i].wait = False;
                                break;
                            }
                    }
                    return 0L;

        case( WM_DDE_REQUEST ):
                    hiword = HIWORD(lArg);
                    if( hiword )
                    {   item = GetItemNumber( (ATOM)hiword );
                    } else item = 0;

                    if( !item )
                    {   GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                        return 0L;
                    }

                    format = (item==AdvImage+1)? CF_DIB : CF_TEXT; 
                    if( format == (int)LOWORD(lArg) )
                    {   SendItemData(hWin,(HWND)wArg,DDEColdLink,item,0);
                        GlobalDeleteAtom( (ATOM)hiword );
                    } else GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                    return 0L;

        case( WM_DDE_UNADVISE ):
                    hiword = HIWORD(lArg);
                    if( hiword )
                    {   item = GetItemNumber( (ATOM)hiword );
                        if( !item )
                        {   GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                            return 0L;
                        }
                    } else item = 0;

                    flag = False;
                    for( i=0; i<MaxDDEAdviseNum; i++ )
                        if( (DDEAdviseData[i].server==hWin) &&
                            ( !hiword || DDEAdviseData[i].item==(char)item ) )
                        {   if( DDEAdviseData[i].wait )
                            {   GlobalDeleteAtom(DDEAdviseData[i].atom);
                                GlobalFree(DDEAdviseData[i].data);
                            }
                            DDEAdviseData[i].server = NULL;
                            DDEAdviseCount--;
                            flag = True;
                        }
                    GenerateDDEReply((HWND)wArg,hWin,flag,hiword);
                    return 0L;

        case( WM_DDE_ADVISE ):
#ifdef _WIN32
                    if( UnpackDDElParam(uMsg,lArg,&loword,&hiword) )
                    {   FreeDDElParam(uMsg,lArg);
                    } else return 0L;
#else
                    loword = LOWORD(lArg);
                    hiword = HIWORD(lArg);
#endif
                    if( hiword )
                    {   item = GetItemNumber( (ATOM)hiword );
                    } else item = 0;

                    if( !item || (DDEAdviseCount==MaxDDEAdviseNum ) )
                    {   GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                        return 0L;
                    }

                    /* Check for established link! */
                    for( i=0; i<MaxDDEAdviseNum; i++ )
                        if( (DDEAdviseData[i].server==hWin) &&
                            (DDEAdviseData[i].item==(char)item) )
                            break;

                    if( i < MaxDDEAdviseNum )
                    {   /* Should we reuse the existing advise? */
                        GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                        return 0L;
                    }

                    options = (DDEADVISE FAR*)GlobalLock((HGLOBAL)loword);
                    if( !options )
                    {   GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                        return 0L;
                    }

                    format = (item==AdvImage+1)? CF_DIB : CF_TEXT;
                    if( options->cfFormat == format ) 
                    {   for( i=0; i<MaxDDEConvNum; i++ )
                           if( DDEConvData[i].server == hWin )
                           {   hDest = DDEConvData[i].client;
                               break;
                           }

                       for( i=0; i<MaxDDEAdviseNum; i++ )
                           if( !DDEAdviseData[i].server )
                               break;

                       DDEAdviseData[i].server = hWin;
                       DDEAdviseData[i].client = hDest;
                       DDEAdviseData[i].atom = hiword;
                       DDEAdviseData[i].wait = False;
                       DDEAdviseData[i].item = item;
                       DDEAdviseCount++;

                       if( options->fDeferUpd )
                       {      DDEAdviseData[i].mode = DDEWarmLink;
                       } else if( options->fAckReq )
                       {      DDEAdviseData[i].mode = DDEAckLink;
                       } else DDEAdviseData[i].mode = DDEHotLink;

                       GenerateDDEReply((HWND)wArg,hWin,True,hiword);

                /* We could advise client of the current item value!      */
                /* SendItemData(hWin,hDest,DDEAdviseData[i].mode,item,i); */

                    } else GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                    GlobalUnlock((HGLOBAL)loword);
                    return 0L;

        case( WM_DDE_EXECUTE ):
#ifdef _WIN32
                    hiword = lArg;
#else
                    hiword = HIWORD(lArg);
#endif
                    cmnd = (char __huge*)GlobalLock((HANDLE)hiword);
                    if( !cmnd )
                    {   GenerateDDEReply((HWND)wArg,hWin,False,hiword);
                        return 0L;
                    }

                    stat = ExecuteIPCCommand( cmnd );
                    GlobalUnlock((HANDLE)hiword);
                    GenerateDDEReply((HWND)wArg,hWin,stat,hiword);

                    if( (stat==IPC_Quit) || (stat==IPC_Exit) )
                        RasMolExit();

                    if( ReDrawFlag )
                        RefreshScreen();
                    if( !CommandActive )
                        ResetCommandLine(0);
                    return 0L;

        case( WM_DDE_TERMINATE ):
                    /* Destroy all Hot/Warm Links */
                    for( i=0; i<MaxDDEAdviseNum; i++ )
                        if( DDEAdviseData[i].server == hWin )
                        {   DDEAdviseData[i].server = NULL;
                            if( DDEAdviseData[i].wait )
                            {   GlobalDeleteAtom(DDEAdviseData[i].atom);
                                GlobalFree(DDEAdviseData[i].data);
                            }
                            DDEAdviseCount--;
                        }

                    /* Remove the Conversation */
                    for( i=0; i<MaxDDEConvNum; i++ )
                        if( DDEConvData[i].server == hWin )
                        {   ptr = &DDEConvData[i];
                            if( !ptr->closed )
                                PostMessage( ptr->client, WM_DDE_TERMINATE,
                                             (WPARAM)ptr->server, 0L );
                            DestroyWindow( ptr->server );
                            ptr->server = NULL;
                            DDEConvCount--;
                            break;
                        }
                    return 0L;

        default:    return DefWindowProc(hWin,uMsg,wArg,lArg);
    }
}
#endif


static void ResizeTerminal( int x, int y )
{
    register int rows, cols;
    register int sr, er;

    register HBRUSH hBr;
    register HDC hDC;
    auto RECT rect;
    
    if( x > CharWide )
    {   cols = x/CharWide;
    } else cols = 1;
    
    if( y > CharHigh )
    {   rows = y/CharHigh;
    } else rows = 1;

    /* Scroll to bottom! */
    if( ScrlStart )
        SetTermScroll( ScrlMax );

    if( rows < TermRows )
    {   if( TermYPos >= rows )
        {   CmndStart += (TermYPos - rows) + 1;
            if( CmndStart >= CmndRows )
                CmndStart -= CmndRows;
            TermYPos = rows - 1;
        
            hDC = GetDC(CmndWin);
            GetClientRect(CmndWin,&rect);
            hBr = CreateSolidBrush(GetSysColor(COLOR_WINDOW));
            FillRect(hDC,&rect,hBr);
            ReleaseDC(CmndWin,hDC);
        } 

    } else if( rows > TermRows )
    {   sr = TermRows + CmndStart;
        if( sr >= CmndRows )
            sr -= CmndRows;
            
        er = CmndStart + rows;
        if( er >= CmndRows )
            er -= CmndRows;
            
        do {
            _fmemset(TermScreen+sr*CmndCols,' ',CmndCols);
            sr++; if( sr == CmndRows ) sr = 0;
        } while( sr != er );
    }
    
    InvalidateRect(CmndWin,NULL,False);
    if( cols > CmndCols )
    {   TermCols = CmndCols;
    } else TermCols = cols;
    TermRows = rows;
}


#ifdef _WIN32
LRESULT CALLBACK CmndCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#else
LONG FAR PASCAL CmndCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#endif
{
    register int row;
    
    switch(uMsg)
    {   case(WM_CLOSE):       DestroyWindow(CanvWin);
                              DestroyWindow(CmndWin);
                              CommandActive = True;
                              ReDrawFlag = False;
                              break;

        case(WM_DESTROY):     /* Destroy RasWin */
                              PostQuitMessage(0);
                              break;

         case(WM_SYSCHAR):    if( lArg & (1L<<29) )  /* ALT-key pressed? */
                                  return(DefWindowProc(hWin,uMsg,wArg,lArg));

         case(WM_CHAR):       if( ProcessCharacter(LOBYTE(wArg)) )
                                  if( ExecuteCommand() )
                                      RasMolExit();
                              break;

         case(WM_PAINT):      PaintScreen();
                              return(0L);
                              
        case(WM_SYSKEYDOWN):
        case(WM_KEYDOWN):     switch(LOBYTE(wArg))
                              {   case(0x23): ProcessCharacter(0x05); break;
                                  case(0x24): ProcessCharacter(0x01); break;
                                  case(0x25): ProcessCharacter(0x02); break;
                                  case(0x26): ProcessCharacter(0x10); break;
                                  case(0x27): ProcessCharacter(0x06); break;
                                  case(0x28): ProcessCharacter(0x0e); break;
                                  case(0x2e): ProcessCharacter(0x04); break;
                                  
                                  default:
                                  return(DefWindowProc(hWin,uMsg,wArg,lArg));
                              }
                              break;
        
         case(WM_SETFOCUS):   if( !TermCursor )
                              {   CreateCaret(hWin,NULL,CharWide,CharHigh);
                                  TermCursor = True;
                              }
                              
                              row = TermYPos + ScrlStart;
                              if( row < TermRows )
                              {   SetCaretPos(TermXPos*CharWide,row*CharHigh);
                                  ShowCaret(hWin);
                              } else HideCaret(hWin);
                              return(0L);
                              
         case(WM_SIZE):       if( wArg != SIZE_MINIMIZED )
                                  ResizeTerminal(LOWORD(lArg),HIWORD(lArg));
                              return(0L);
                              
         case(WM_KILLFOCUS):  if( TermCursor )
                              {   TermCursor=False;
                                  HideCaret(hWin);
                                  DestroyCaret();
                              }
                              return(0L);

         case(WM_VSCROLL):
#ifdef _WIN32
                              switch( LOWORD(wArg) )
#else
                              switch( wArg )
#endif
                              {  case(SB_TOP):    SetTermScroll(0);  break;
                                 case(SB_BOTTOM): SetTermScroll(ScrlMax);  
                                                  break;
                                 
                                 case(SB_LINEUP):   
                                     if( ScrlStart < ScrlMax )
                                         SetTermScroll((ScrlMax-ScrlStart)-1);
                                     break;
                                     
                                 case(SB_LINEDOWN):
                                     if( ScrlStart > 0 )
                                         SetTermScroll((ScrlMax-ScrlStart)+1);
                                     break;
                                     
                                 case(SB_PAGEUP):
                                     if( ScrlStart < (ScrlMax-10) )
                                     {   SetTermScroll((ScrlMax-ScrlStart)-10);
                                     } else SetTermScroll(0);
                                     break;
                                     
                                 case(SB_PAGEDOWN):
                                     if( ScrlStart > 10 )
                                     {   SetTermScroll((ScrlMax-ScrlStart)+10);
                                     } else SetTermScroll(ScrlMax);
                                     break;
                                     
                                 case(SB_THUMBTRACK):
                                 case(SB_THUMBPOSITION):
#ifdef _WIN32
                                     SetTermScroll(HIWORD(wArg));
#else
                                     SetTermScroll(LOWORD(lArg));
#endif
                                     break;
                              }
                              break;
                                                            
         default:  return DefWindowProc(hWin,uMsg,wArg,lArg);
    }

    if( ReDrawFlag )
        RefreshScreen();
    if( !CommandActive )
        ResetCommandLine(0);
    return 0L;
}


static void LoadInputFile( void )
{
    ofn1.lpstrDefExt = "PDB";
    ofn1.nFilterIndex = 1;
    *fnamebuf = '\0';

    if( GetOpenFileName(&ofn1) )
    {   switch( ofn1.nFilterIndex )
        {   case(1): FetchFile(FormatPDB,False,fnamebuf);     break;
            case(2): FetchFile(FormatAlchemy,False,fnamebuf); break;
            case(3): FetchFile(FormatMol2,False,fnamebuf);    break;
            case(4): FetchFile(FormatMDL,False,fnamebuf);     break;
            case(5): FetchFile(FormatXYZ,False,fnamebuf);     break;
            case(6): FetchFile(FormatCharmm,False,fnamebuf);  break;
            case(7): FetchFile(FormatMOPAC,False,fnamebuf);   break;
        }
        DefaultRepresentation();
    }
}


static void SaveOutputFile( int format )
{
    register char *ext;
    register int num;

    switch( format )
    {   case(IDM_BMP):   ext="BMP";  num=1;   break;
        case(IDM_GIF):   ext="GIF";  num=2;   break;
        case(IDM_EPSF):  ext="PS";   num=3;   break;
        case(IDM_PPM):   ext="PPM";  num=6;   break;
        case(IDM_RAST):  ext="RAS";  num=10;  break;
        default: return;
    }

    ofn2.nFilterIndex = num;
    ofn2.lpstrDefExt = ext;
    *fnamebuf = '\0';
    
/*  Default Filename   
 *  dst = fnamebuf;
 *  for( src="RASWIN."; *src; src++ ) *dst++ = *src;
 *  for( src=ext; *src; src++ ) *dst++ = *src;
 *  *dst++ = '\0';
 */
    
    if( GetSaveFileName(&ofn2) )    
        switch( ofn2.nFilterIndex )
        {   case(1):  WriteBMPFile(fnamebuf);             break;
            case(2):  WriteGIFFile(fnamebuf);             break;
            case(3):  WriteEPSFFile(fnamebuf,True,True);  break;
            case(4):  WriteEPSFFile(fnamebuf,False,True); break;
            case(5):  WriteVectPSFile(fnamebuf);          break;
            case(6):  WritePPMFile(fnamebuf,True);        break;
            case(7):  WritePPMFile(fnamebuf,False);       break;
            case(8):  WritePICTFile(fnamebuf);            break;
            case(9):  WriteIRISFile(fnamebuf);            break;
            case(10): WriteRastFile(fnamebuf,True);       break;
            case(11): WriteRastFile(fnamebuf,False);      break;
        }
}


static void HandlePrintSetUp( void )
{
    auto PRINTDLG pd;

    memset(&pd,0,sizeof(PRINTDLG));
    pd.lStructSize = sizeof(PRINTDLG);
    pd.hwndOwner = CanvWin;
    pd.Flags = PD_PRINTSETUP;

    PrintDlg(&pd);

    if( pd.hDevNames ) GlobalFree(pd.hDevNames);
    if( pd.hDevMode )  GlobalFree(pd.hDevMode);
}


static BOOL HandleMenu( WPARAM option )
{
#ifndef _WIN32
    register FARPROC lpProc;
#endif
    register char *ptr;
    register int mask;
   
    switch(option)
    {   /* File Menu */
        case(IDM_OPEN):   if( !Database )
                              LoadInputFile();
                          break;
                          
        case(IDM_INFO):
#ifdef _WIN32
                          DialogBox(hInstance,"InfoBox",CanvWin,InfoCallB);
#else   
                          lpProc = MakeProcInstance(InfoCallB,hInstance);
                          DialogBox(hInstance,"InfoBox",CanvWin,lpProc);
                          FreeProcInstance(lpProc);
#endif
                          break;
                          
        case(IDM_CLOSE):  ZapDatabase();
                          break;

        case(IDM_PRINT):  if( !PrintImage() )
                          {   if( CommandActive )
                                  WriteChar('\n');
                              WriteString("Warning: No suitable printer!\n");
                              CommandActive = False;
                          }
                          break;
        
        case(IDM_SETUP):  HandlePrintSetUp();
                          break;

        case(IDM_EXIT):   PostMessage(CanvWin,WM_CLOSE,0,0L);
                          break;

        /* Edit Menu */
        case(IDM_SELECT): mask = NormAtomFlag;
                          if( HetaGroups ) mask |= HeteroFlag;
                          if( Hydrogens )  mask |= HydrogenFlag;
                          SelectZone(mask);
                          break;

        case(IDM_COPY):   if( !ClipboardImage() )
                          {   if( CommandActive )
                                  WriteChar('\n');
                              WriteString("Unable to copy to clipboard!\n");
                              CommandActive = False;
                          }
                          break;
        
        /* Help Menu */
        case(IDM_ABOUT):
#ifdef _WIN32
                          DialogBox(hInstance,"AboutBox",CanvWin,AboutCallB);
#else
                          lpProc = MakeProcInstance(AboutCallB,hInstance);
                          DialogBox(hInstance,"AboutBox",CanvWin,lpProc);
                          FreeProcInstance(lpProc);
#endif
                          break;

        case(IDM_HELP):   if( getcwd(fnamebuf,100) )
                          {   ptr = fnamebuf;
                              while( *ptr ) ptr++;
                              if( (ptr!=fnamebuf) && (ptr[-1]!='\\') ) 
                                  *ptr++ = '\\';

                              strcpy(ptr,"RASWIN.HLP");
                              WinHelp(CanvWin,fnamebuf,HELP_INDEX,0L);
                          }
                          break;
       
        /* Display Menu */
        case(IDM_WIREFRAME):  DisableSpacefill();
                              EnableWireframe(WireFlag,0);
                              SetRibbonStatus(False,0,0);
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_BACKBONE):   DisableSpacefill();
                              DisableWireframe();
                              SetRibbonStatus(False,0,0);
                              EnableBackbone(CylinderFlag,80);
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_STICKS):     DisableSpacefill();
                              if( MainAtomCount<256 )
                              {   EnableWireframe(CylinderFlag,40);
                              } else EnableWireframe(CylinderFlag,80);
                              SetRibbonStatus(False,0,0);
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_SPHERES):    SetVanWaalRadius();
                              DisableWireframe();
                              SetRibbonStatus(False,0,0);
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_BALLSTICK):  SetRadiusValue(120);
                              EnableWireframe(CylinderFlag,40);
                              SetRibbonStatus(False,0,0);
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_RIBBONS):    DisableSpacefill();
                              DisableWireframe();
                              SetRibbonStatus(True,RibbonFlag,0);
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_STRANDS):    DisableSpacefill();
                              DisableWireframe();
                              SetRibbonStatus(True,StrandFlag,0);
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        case(IDM_CARTOONS):   /* Cartoons */
                              DisableSpacefill();
                              DisableWireframe();
                              SetRibbonCartoons();
                              DisableBackbone();
                              ReDrawFlag |= RFRefresh;
                              break;

        /* Colours Menu */
        case(IDM_MONO):     MonoColourAttrib(255,255,255);
                            ReDrawFlag |= RFColour;  break;
        case(IDM_CPK):      CPKColourAttrib();
                            ReDrawFlag |= RFColour;  break;
        case(IDM_SHAPELY):  ShapelyColourAttrib();
                            ReDrawFlag |= RFColour;  break;
        case(IDM_STRUCT):   StructColourAttrib();
                            ReDrawFlag |= RFColour;  break;
        case(IDM_GROUP):    ScaleColourAttrib( GroupAttr );
                            ReDrawFlag |= RFColour;  break;
        case(IDM_CHAIN):    ScaleColourAttrib( ChainAttr );
                            ReDrawFlag |= RFColour;  break;
        case(IDM_TEMPER):   ScaleColourAttrib( TempAttr );
                            ReDrawFlag |= RFColour;  break;
        case(IDM_USER):     UserMaskAttrib(MaskColourFlag);
                            ReDrawFlag |= RFColour;  break;
        
        /* Options Menu */
        case(IDM_SLAB):      ReDrawFlag |= RFRefresh;
                             UseSlabPlane = !UseSlabPlane;
                             if( UseSlabPlane )
                                 UseShadow = False;
                             break;

        case(IDM_HYDROGEN):  mask = NormAtomFlag;
                             if( HetaGroups )
                                 mask |= HeteroFlag;
                             Hydrogens = !Hydrogens;
                             ReDrawFlag |= RFRefresh;

                             if( Hydrogens )
                             {      SelectZone(mask|HydrogenFlag);
                             } else RestrictZone(mask);
                             break;
        
        case(IDM_HETERO):    mask = NormAtomFlag;
                             if( Hydrogens )
                                 mask |= HydrogenFlag;
                             HetaGroups = !HetaGroups;
                             ReDrawFlag |= RFRefresh;

                             if( HetaGroups )
                             {      SelectZone(mask|HeteroFlag);
                             } else RestrictZone(mask);
                             break;
        
        case(IDM_SPECULAR):  FakeSpecular = !FakeSpecular;
                             ReDrawFlag |= RFColour;
                             break;
        
        case(IDM_SHADOW):    ReDrawFlag |= RFRefresh;
                             UseShadow = !UseShadow;
                             if( UseShadow )
                             {   ReviseInvMatrix();
                                 VoxelsClean = False;
                                 UseSlabPlane = False;
                                 ReAllocBuffers();
                             }
                             break;

        case(IDM_STEREO):    /* Stereo */
                             if( UseStereo )
                             {   SetStereoMode(False);
                             } else SetStereoMode(True);
                             ReDrawFlag |= RFRefresh;
                             break;

        case(IDM_LABELS):    /* Labels */
                             LabelOptFlag = !LabelOptFlag;
                             DefaultLabels(LabelOptFlag);
                             ReDrawFlag |= RFRefresh;
                             break;

        /* Save Menu */
        case(IDM_BMP):   case(IDM_GIF):
        case(IDM_EPSF):  case(IDM_PPM):
        case(IDM_RAST):    SaveOutputFile( option ); 
                           break;
        
        default:  return FALSE;
    }
    return TRUE;
}    


#ifdef DDEIPC
static void InitiateServer( HWND hWinCli, LPARAM lParam )
{
    register ATOM aTopicIn, aTopicOut;
    register ATOM aApplIn, aApplOut;
    register HWND hWinServ;
    
    auto char TopicName[16];
    auto char ApplName[16];
    register int i;

    if( DDEConvCount == MaxDDEConvNum )
        return;
            
    if( aApplIn = LOWORD(lParam) )
    {   GlobalGetAtomName(aApplIn,ApplName,14);
        if( stricmp(ApplName,"RasWin") &&
            stricmp(ApplName,"RasWin32") )
            return;
    } else return;
    
    if( aTopicIn = HIWORD(lParam) )
    {   GlobalGetAtomName(aTopicIn,TopicName,14);
        /* Test for Valid Topic */
        /* if( _stricmp(Topic,"System") &&
         *     _stricmp(Topic,"RemoteControl") )
         * return;
         */
    } else *TopicName = '\0';
   
    hWinServ = CreateWindow("RasDDEClass","RasWinDDE",
                            WS_CHILD, 0, 0, 0, 0,
                            CanvWin, NULL, hInstance, NULL );
    if( !hWinServ ) return;
         
    for( i=0; i<MaxDDEConvNum; i++ )
        if( !DDEConvData[i].server )
            break;
            
    DDEConvData[i].server = hWinServ;
    DDEConvData[i].client = hWinCli;
    DDEConvData[i].closed = False;
    DDEConvCount++;       
          
    /* Main DDE Server */       
    aTopicOut = GlobalAddAtom("System");
#ifdef _WIN32
    aApplOut = GlobalAddAtom("RasWin32");
    SendMessage( hWinCli, WM_DDE_ACK, (WPARAM)hWinServ,
                 MAKELPARAM(aApplOut,aTopicOut) );
#else
    aApplOut = GlobalAddAtom("RasWin");
    SendMessage( hWinCli, WM_DDE_ACK, (WPARAM)hWinServ,
                 MAKELONG(aApplOut,aTopicOut) ); 
#endif
}
#endif


#ifdef SOCKETS
static int IsIPCAdviseRequest( char *ptr, int conv )
{
    auto char item[34];
    register char *dst;
    register char *src;
    register int i,ch;
    register int flag;

    if( !strncmp(ptr,"Advise:",7) )
    {   src = ptr+7;
        flag = True;
    } else if( !strncmp(ptr,"Unadvise:",9) )
    {   src = ptr+9;
        flag = False;
    } else return False;

    while( True )
    {   ch = *src++;
        if( isspace(ch) )
            continue;

        if( isalpha(ch) )
        {   dst = item;
            *dst++ = ch;
            while( IsIdentChar(*src) )
            {   if( dst < item+32 )
                {   *dst++ = *src++;
                } else src++;
            }
            *dst = '\0';

            for( i=0; i<AdvItemCount; i++ )
                if( !stricmp(item,AdviseMap[i].name) )
                {   if( flag )
                    {      IPCConvData[conv].advise |=  AdviseMap[i].bitmask;
                    } else IPCConvData[conv].advise &= ~AdviseMap[i].bitmask;
                    break;
                }

           /* Warning: Unknown Advise Item! */
        } else if( ch != ',' )
            break;
    }
    return True;
}


static void OpenSocket( HWND hWin )
{
    auto int length;
    auto struct sockaddr_in addr;
    auto WSADATA wsadata;
    register WORD vers;
    register int i;

    UseSockets = False;
    SocketNo = INVALID_SOCKET;

    vers = MAKEWORD(1,1);
    if( WSAStartup(vers,&wsadata) )
        return;

    if( LOBYTE(wsadata.wVersion) != 1 ||
        HIBYTE(wsadata.wVersion) != 1 )
    {   WSACleanup();
        return;
    }

    SocketNo = socket(AF_INET,SOCK_STREAM,0);
    if( SocketNo == INVALID_SOCKET )
    {   WSACleanup();
        return;
    }

    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);
    addr.sin_port = htons((short)ServerPort);

    if( bind(SocketNo,(struct sockaddr __far*)&addr,sizeof(addr)) )
    {   closesocket(SocketNo);
        WSACleanup();
        return;
    }

    if( !ServerPort )
    {   length = sizeof(addr);
        if( !getsockname(SocketNo,(struct sockaddr __far*)&addr,&length) )
        {   ServerPort = ntohs(addr.sin_port);
            sprintf(Text,"RasMol Server TCP/IP Port: %d\n",ServerPort);
            WriteString(Text);
        }
    }

    UseSockets = True;
    for( i=0; i<MaxIPCConvNum; i++ )
        IPCConvData[i].protocol = 0;

    listen(SocketNo,5);
    WSAAsyncSelect(SocketNo,hWin,WM_WINSOCK,FD_ACCEPT|FD_READ|FD_CLOSE);
}


static int OpenIPCConnection( SOCKET sock )
{
    register int i;

    if( sock == INVALID_SOCKET )
        return False;

    for( i=0; i<MaxIPCConvNum; i++ )
        if( !IPCConvData[i].protocol )
        {   IPCConvData[i].protocol = ProtoRasMol;
            IPCConvData[i].socket = sock;
            IPCConvData[i].advise = AMNone;
            return True;
        }
    closesocket(sock);
    return False;
}


static void CloseIPCConnection( SOCKET sock )
{
    register int i;

    if( sock != INVALID_SOCKET )
    {   for( i=0; i<MaxIPCConvNum; i++ )
            if( IPCConvData[i].protocol && (IPCConvData[i].socket==sock) )
                IPCConvData[i].protocol = 0;
        closesocket(sock);
    }
}


static void HandleSocketData( SOCKET sock )
{
    auto char buffer[4097];
    register char *src;
    register char *dst;
    register int result;
    register int len;
    register int ch;
    register int i;

    if( sock == INVALID_SOCKET )
        return;

    for( i=0; i<MaxIPCConvNum; i++ )
        if( IPCConvData[i].protocol && (IPCConvData[i].socket==sock) )
        {   len = recv( sock, buffer, 4096, 0 );
            if( len > 0 )
            {   buffer[len] = '\0';
                src = dst = buffer;
                while( (ch = *src++) )
                    if( (ch>=' ') && (ch<='~') )
                        *dst++ = ch;
                *dst = '\0';

                if( !IsIPCAdviseRequest(buffer,i) )
                {   result = ExecuteIPCCommand(buffer);
                    if( result == IPC_Exit )
                    {   CloseIPCConnection(sock);
                    } else if( result == IPC_Quit )
                        RasMolExit();
                }
            } else CloseIPCConnection(sock);
            return;
        }

    closesocket(sock);
}
#endif


static int GetStatus( int mask )
{
    register int status;
    
    status = 0;                             
    if( mask & MK_LBUTTON ) status |= MMLft;
    if( mask & MK_MBUTTON ) status |= MMMid;
    if( mask & MK_RBUTTON ) status |= MMRgt;
    if( mask & MK_CONTROL ) status |= MMCtl;          
    if( mask & MK_SHIFT )   status |= MMSft;
    return status;
}


static void AdjustMenus( WPARAM wArg, LPARAM lArg )
{
    register HMENU hMenu;
    register int status;
    
    if( lArg == 4 )
    {   /* Options Menu */
        hMenu = (HMENU)wArg;

        status = UseSlabPlane ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_SLAB,status);

        status = Hydrogens ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_HYDROGEN,status);

        status = HetaGroups ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_HETERO,status);

        status = FakeSpecular ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_SPECULAR,status);

        status = UseShadow ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_SHADOW,status);

        status = UseStereo ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_STEREO,status);

        status = LabelOptFlag ? MF_CHECKED : MF_UNCHECKED;
        CheckMenuItem(hMenu,IDM_LABELS,status);
    }
}


#ifdef _WIN32
LRESULT CALLBACK MainCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#else
LONG FAR PASCAL MainCallB( HWND hWin, UINT uMsg, WPARAM wArg, LPARAM lArg )
#endif
{
#ifdef EIGHTBIT
    register HPALETTE hCMap;
#endif
    register int pos,status;
    register int x,y;

    register COLORREF col;    
    register HANDLE hand;
    register HDC hMemDC;
    register HDC hDC;
    auto PAINTSTRUCT ps;
    auto RECT rc;
    
    CanvWin = hWin;

    switch(uMsg)
    {   case(WM_DROPFILES):   /* Disable Drag & Drop */
                              if( IsPaused ) break;

                              ZapDatabase();
                              *fnamebuf = '\0';
                              DragQueryFile((HDROP)wArg,0,fnamebuf,127);
                              FetchFile(FormatPDB,False,fnamebuf);
                              DefaultRepresentation();
                              DragFinish((HDROP)wArg);
                              break;

        case(WM_DESTROY):     /* Destroy RasWin */
                              DragAcceptFiles(CanvWin,FALSE);
                              PostQuitMessage(0);
                              break;

        case(WM_CLOSE):       DestroyWindow(CanvWin);
                              DestroyWindow(CmndWin);
                              break;

        case(WM_ACTIVATE):    
#ifdef _WIN32
                              if( !LOWORD(wArg) || HIWORD(lArg) ) break;
#else
                              if( !wArg || HIWORD(lArg) ) break;
#endif

#ifdef EIGHTBIT
        case(WM_QUERYNEWPALETTE):
                              if( ColourMap )
                              {   hDC = GetDC(hWin);
                                  hCMap = SelectPalette(hDC,ColourMap,False);
                                  status = RealizePalette(hDC);
                                  if( hCMap ) SelectPalette(hDC,hCMap,False);
                                  ReleaseDC(hWin,hDC);
                                  
                                  if( status )
                                  {   InvalidateRect(hWin,NULL,True);
                                      return True;
                                  }
                              }
                              return 0L;
                              
        case(WM_PALETTECHANGED):
                              if( ColourMap && ((HWND)wArg != hWin) )
                              {   hDC = GetDC(hWin);
                                  hCMap = SelectPalette(hDC,ColourMap,False);
                                  if( RealizePalette(hDC) )
                                      InvalidateRect(hWin,NULL,True);
                                  if( hCMap ) SelectPalette(hDC,hCMap,False);
                                  ReleaseDC(hWin,hDC);
                              }
                              return 0L;
#endif

        case(WM_INITMENUPOPUP):  /* Initialise Checks */
                              AdjustMenus( wArg, lArg );
                              return 0L;                      

        case(WM_SIZE):        if( wArg != SIZE_MINIMIZED )
                              {   GetClientRect(hWin,&rc);
                                  YRange = rc.bottom;
                                  XRange = rc.right;
                              
                                  /* Ensure Long Aligned */
                                  XRange = (XRange+3)&~3;
                              
                                  Range = MinFun(XRange,YRange);
                                  ReDrawFlag |= RFReSize;
                                  HRange = YRange>>1;
                                  WRange = XRange>>1;
                                  ClearImage();
                              }
                              break;

        case(WM_HSCROLL):     /* Horizontal Scroll */
                              pos = GetScrollPos(hWin,SB_HORZ);
#ifdef _WIN32
                              switch( LOWORD(wArg) )
#else
                              switch( wArg )
#endif
                              {   case(SB_LINEDOWN):  pos += 5;   break;
                                  case(SB_PAGEDOWN):  pos += 10;  break;
                                  case(SB_PAGEUP):    pos -= 10;  break;
                                  case(SB_LINEUP):    pos -= 5;   break;
                                  default:            return(0L);

                                  case(SB_THUMBTRACK):
                                  case(SB_THUMBPOSITION):
#ifdef _WIN32
                                             pos = HIWORD(wArg);
#else
                                             pos = LOWORD(lArg);
                                             break;
#endif
                              }
                              
                              if( pos>100 ) 
                              {   pos -= 100;
                              } else if( pos<0 ) 
                                  pos += 100; 
                             
                              SetScrollPos(hWin,SB_HORZ,pos,TRUE);
                              DialValue[1] = (pos/50.0)-1.0;
                              ReDrawFlag |= RFRotateY;
                              break;

        case(WM_VSCROLL):     /* Vertical Scroll */
                              pos = GetScrollPos(hWin,SB_VERT);
#ifdef _WIN32
                              switch( LOWORD(wArg) )
#else
                              switch( wArg )
#endif
                              {   case(SB_LINEDOWN):  pos += 5;   break;
                                  case(SB_PAGEDOWN):  pos += 10;  break;
                                  case(SB_PAGEUP):    pos -= 10;  break;
                                  case(SB_LINEUP):    pos -= 5;   break;
                                  default:            return(0L);

                                  case(SB_THUMBTRACK):
                                  case(SB_THUMBPOSITION):
#ifdef _WIN32
                                             pos = HIWORD(wArg);
#else
                                             pos = LOWORD(lArg);
#endif
                                             break;
                              }
                              
                              if( pos>100 ) 
                              {   pos -= 100;
                              } else if( pos<0 ) 
                                  pos += 100; 
                             
                              SetScrollPos(hWin,SB_VERT,pos,TRUE);
                              DialValue[0] = 1.0-(pos/50.0);
                              ReDrawFlag |= RFRotateX;
                              break;

        case(WM_LBUTTONDOWN): 
        case(WM_MBUTTONDOWN):
        case(WM_RBUTTONDOWN): x = (short)LOWORD(lArg);
                              y = YRange-(short)HIWORD(lArg);
                              status = GetStatus(wArg);
                              ProcessMouseDown(x,y,status);
                              break;
        
        case(WM_LBUTTONUP):
        case(WM_MBUTTONUP):
        case(WM_RBUTTONUP):   /* Mouse Buttons */
                              x = (short)LOWORD(lArg);
                              y = YRange-(short)HIWORD(lArg);
                              status = GetStatus(wArg);
                              ProcessMouseUp(x,y,status);
                              break;


        case(WM_MOUSEMOVE):   /* Mouse Movement */
                              x = (short)LOWORD(lArg);
                              y = YRange-(short)HIWORD(lArg);
                              status = GetStatus(wArg);
                              ProcessMouseMove(x,y,status);
                              break;
                                      
        case(WM_SETFOCUS):    /* Obtain Window Focus */ 
        case(WM_KILLFOCUS):   /* Release Window Focus */
                              SendMessage(CmndWin,uMsg,wArg,lArg);     
                              return( 0L );
        
        case(WM_PAINT):       hDC = BeginPaint(hWin,&ps);
                              SetBkMode(hDC,TRANSPARENT);
                              if( PixMap )
                              {   
#ifdef EIGHTBIT
                                  hCMap = SelectPalette(hDC,ColourMap,False);
                                  RealizePalette(hDC);
#endif
#ifdef _WIN32
                                  SetWindowOrgEx(hDC,0,0,NULL);
#else
                                  SetWindowOrg(hDC,0,0);
#endif
                                  hMemDC = CreateCompatibleDC(hDC);
                                  SelectObject(hMemDC,PixMap);
                                  BitBlt(hDC,0,0,XRange,YRange,
                                         hMemDC,0,0,SRCCOPY);
#ifdef EIGHTBIT
                                  SelectPalette(hDC,hCMap,False);
#endif
                                  DeleteDC(hMemDC);
                              } else /* Erase Update Region */
                              {
#ifdef EIGHTBIT
                                   if( ColourMap )
                                   {   hCMap=SelectPalette(hDC,ColourMap,0);
                                       RealizePalette(hDC);
                                       col = RGB(BackR,BackG,BackB);
                                       hand = CreateSolidBrush(col);
                                       GetUpdateRect(hWin,&rc,False);
                                       FillRect(hDC,&rc,(HBRUSH)hand);
                                       if( hCMap )
                                           SelectPalette(hDC,hCMap,False);
                                   } else
                                   {   col = RGB(BackR,BackG,BackB);
                                       hand = CreateSolidBrush(col);
                                       GetUpdateRect(hWin,&rc,False);
                                       FillRect(hDC,&rc,(HBRUSH)hand);
                                   }
#else
                                   col = RGB(BackR,BackG,BackB);
                                   hand = CreateSolidBrush(col);
                                   GetUpdateRect(hWin,&rc,False);
                                   FillRect(hDC,&rc,hand);
#endif
                                   DeleteObject(hand);
                              }
                              EndPaint(hWin,&ps);
#ifdef DDEIPC
                              if( !RasWinDDEReady )
                              {   RasWinDDEReady = True;
                                  AdviseUpdate(-4);
                              }
#endif
                              return 0L;
        
        case(WM_SYSCHAR):     if( lArg & (1L<<29) )  /* ALT-key pressed? */
                                  return(DefWindowProc(hWin,uMsg,wArg,lArg));
        case(WM_CHAR):        if( ProcessCharacter(LOBYTE(wArg)) )
                                  if( ExecuteCommand() )
                                      RasMolExit();
                              break;

        case(WM_SYSKEYDOWN):
        case(WM_KEYDOWN):     switch(LOBYTE(wArg))
                              {   case(0x23): ProcessCharacter(0x05); break;
                                  case(0x24): ProcessCharacter(0x01); break;
                                  case(0x25): ProcessCharacter(0x02); break;
                                  case(0x26): ProcessCharacter(0x10); break;
                                  case(0x27): ProcessCharacter(0x06); break;
                                  case(0x28): ProcessCharacter(0x0e); break;
                                  case(0x2e): ProcessCharacter(0x04); break;
                                  
                                  default:
                                  return(DefWindowProc(hWin,uMsg,wArg,lArg));
                              }
                              break;

        case(WM_RENDERALLFORMATS):
                              OpenClipboard(hWin);
                              SendMessage(hWin,WM_RENDERFORMAT,CF_DIB,0L);
                              SendMessage(hWin,WM_RENDERFORMAT,CF_PALETTE,0L);
                              CloseClipboard();
                              return 0L;
                              
        case(WM_RENDERFORMAT):
                              hand = RenderClipboard(wArg);
                              if( hand )
                                  SetClipboardData(wArg,hand);
                              return 0L;
                              

#ifdef DDEIPC
        case(WM_DDE_INITIATE): /* DDE Server Connection */
                              InitiateServer((HWND)wArg,lArg);
                              return 0L;
#endif

#ifdef SOCKETS
        case(WM_WINSOCK):     /* IPC Server Connection */
                              switch( WSAGETSELECTEVENT(lArg) )
                              {   case FD_ACCEPT:
                                      OpenIPCConnection(accept(wArg,0,0));
                                      return 0L;

                                  case FD_READ:
                                      HandleSocketData((SOCKET)wArg);
                                      break;

                                  case FD_CLOSE:
                                      CloseIPCConnection((SOCKET)wArg);
                                      return 0L;

                                  default:
                                      return 0L;
                              }
                              break;
#endif

        case(WM_COMMAND):     
#ifdef _WIN32
                              if( !IsPaused && HandleMenu(LOWORD(wArg)) ) 
                                  break;
#else
                              if( !IsPaused && HandleMenu(wArg) ) break;
#endif                              
        default:              return( DefWindowProc(hWin,uMsg,wArg,lArg) );

    }

    if( ReDrawFlag )
        RefreshScreen();
    if( !CommandActive )
        ResetCommandLine(0);
    return 0L;
}


static int InitialiseApplication( void )
{
    WNDCLASS wc;
    
    wc.hIcon = LoadIcon(hInstance,"RasWinIcon");
    wc.hInstance = hInstance;
    wc.cbWndExtra = 0;
    wc.cbClsExtra= 0;

    /* Canvas Window Class */
    wc.style = 0;
    wc.lpfnWndProc = MainCallB;
    wc.hbrBackground = CreateSolidBrush(RGB(0,0,0));
    wc.hCursor = LoadCursor(hInstance,"RasWinCursor");
    wc.lpszClassName = "RasWinClass";
    wc.lpszMenuName = NULL;

    if( !RegisterClass(&wc) )
        return False;

    /* Terminal Window Class */
    wc.style = CS_NOCLOSE;
    wc.lpfnWndProc = CmndCallB;
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
    wc.hCursor = LoadCursor(NULL,IDC_ARROW);
    wc.lpszClassName = "RasCliClass";
    wc.lpszMenuName = NULL;

    if( !RegisterClass(&wc) )
        return False;

#ifdef DDEIPC
    /* DDE Server Window Class */
    wc.lpfnWndProc = DDECallB;
    wc.lpszClassName = "RasDDEClass";
    wc.hbrBackground = NULL;
    wc.hCursor = NULL;
    wc.hIcon = NULL;

    if( !RegisterClass(&wc) )
        return False;
#endif
    return True;
}


static char *RegisterFormat( char *buffer, char *desc, char *ext )
{
    while( *desc )
        *buffer++ = *desc++;
    *buffer++ = '\0';

    while( *ext )
        *buffer++ = *ext++;
    *buffer++ = '\0';
    return buffer;
}


static void InitFileDialogBoxes( void )
{
    register char *dst;

    dst = ifilters;
    dst = RegisterFormat(dst,"Brookhaven Databank","*.PDB;*.ENT");
    dst = RegisterFormat(dst,"Alchemy File Format","*.ALC;*.MOL");
    dst = RegisterFormat(dst,"Sybyl MOL2 Format","*.SYB;*.ML2;*.SY2;*.MOL");
    dst = RegisterFormat(dst,"MDL Mol File Format","*.MDL;*.MOL");
    dst = RegisterFormat(dst,"MSC (XMol) XYZ Format","*.XYZ");
    dst = RegisterFormat(dst,"CHARMm File Format","*.CHM");
    dst = RegisterFormat(dst,"MOPAC File Format","*.MOP");
    *dst = '\0';

    /* Load File Common Dialog Box */
    ofn1.lStructSize=sizeof(OPENFILENAME);
    ofn1.Flags = OFN_NOCHANGEDIR | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
    ofn1.lpstrFilter = ifilters;
    ofn1.lpstrTitle = "Select Molecular Coordinate File";
    ofn1.lpstrFile = fnamebuf;
    ofn1.nMaxFile = 128;

    ofn1.lpstrCustomFilter = NULL;
    ofn1.lpstrInitialDir = NULL;
    ofn1.lpstrFileTitle = NULL;
    ofn1.hwndOwner = NULL;

    dst = ofilters;
    dst = RegisterFormat(dst,"Microsoft Bitmap","*.BMP");
    dst = RegisterFormat(dst,"CompuServe GIF","*.GIF");
    dst = RegisterFormat(dst,"Colour PostScript","*.PS;*.EPS");
    dst = RegisterFormat(dst,"Mono PostScript","*.PS;*.EPS");
    dst = RegisterFormat(dst,"Vector PostScript","*.PS;*.EPS");
    dst = RegisterFormat(dst,"Raw Portable Pixmap","*.PPM");
    dst = RegisterFormat(dst,"ASCII Portable Pixmap","*.PPM");
    dst = RegisterFormat(dst,"Apple Macintosh PICT","*.PIC");
    dst = RegisterFormat(dst,"Silicon Graphics RGB","*.RGB");
    dst = RegisterFormat(dst,"RLE Sun Rasterfile","*.RAS;*.IM8");
    dst = RegisterFormat(dst,"Sun Rasterfile","*.RAS");
    *dst = '\0';

    /* Save File Common Dialog Box */
    ofn2.lStructSize=sizeof(OPENFILENAME);
    ofn2.Flags = OFN_NOCHANGEDIR | OFN_HIDEREADONLY | OFN_NOREADONLYRETURN
               | OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
    ofn2.lpstrFilter = ofilters;
    ofn2.lpstrTitle = "Select Graphics Ouptut File";
    ofn2.lpstrFile = fnamebuf;
    ofn2.nMaxFile = 128;

    ofn2.lpstrCustomFilter = NULL;
    ofn2.lpstrInitialDir = NULL;
    ofn2.lpstrFileTitle = NULL;
    ofn2.hwndOwner = NULL;
}


static void InitWindowsProfiles( void )
{
#ifdef _WIN32
    if( !GetProfileString("extensions","pdb","",Text,128) )
        WriteProfileString("extensions","pdb","raswin32.exe ^.pdb");
    if( !GetProfileString("extensions","ent","",Text,128) )
        WriteProfileString("extensions","ent","raswin32.exe ^.ent");
#else
    if( !GetProfileString("extensions","pdb","",Text,128) )
        WriteProfileString("extensions","pdb","raswin.exe ^.pdb");
    if( !GetProfileString("extensions","ent","",Text,128) )
        WriteProfileString("extensions","ent","raswin.exe ^.ent");
#endif
#ifdef DDEIPC
    DDETimeOut = GetPrivateProfileInt("RasWin","DDETimeOut",
                                      DefaultDDETimeOut,"RASWIN.INI");
#endif
}


static void InitDefaultValues( void )
{
#ifdef DDEIPC
    register int i;

    DDEConvCount = 0;
    DDEAdviseCount = 0;
    for( i=0; i<MaxDDEConvNum; i++ )
        DDEConvData[i].server = NULL;
    for( i=0; i<MaxDDEAdviseNum; i++ )
        DDEAdviseData[i].server = NULL;
    RasWinDDEReady = True;
#endif

    Interactive = True;

    fnamebuf[0] = '\0';
    snamebuf[0] = '\0';

    CalcBondsFlag = True;
    LabelOptFlag = False;
    FileFormat = FormatPDB;
    AllowWrite = False;

#ifdef SOCKETS
    ServerPort = 21069;
#endif
}


#define FORMATOPTMAX   15
static struct {
        char *ident;
        int format;
        int len;
    } FormatOpt[FORMATOPTMAX] = {
            { "alchemy",    FormatAlchemy,   7 },
            { "biosym",     FormatBiosym,    6 },
            { "cif",        FormatCIF,       3 },
            { "charmm",     FormatCharmm,    6 },
            { "fdat",       FormatFDAT,      4 },
            { "gaussian",   FormatGaussian,  8 },
            { "macromodel", FormatMacroMod, 10 },
            { "mdl",        FormatMDL,       3 },
            { "mmdb",       FormatMMDB,      4 },
            { "mol2",       FormatMol2,      4 },
            { "mopac",      FormatMOPAC,     5 },
            { "nmrpdb",     FormatNMRPDB,    6 },
            { "pdb",        FormatPDB,       3 },
            { "shelx",      FormatSHELX,     5 },
            { "xyz",        FormatXYZ,       3 }
                                };


static int ProcessOptions( char __far *ptr )
{
    register char *dst;
    register int i;

    while( *ptr )
    {   if( (*ptr==' ') || (*ptr=='=') )
        {   ptr++;
        } else if( (*ptr=='/') || (*ptr=='-') )
        {   ptr++;
            for( i=0; i<FORMATOPTMAX; i++ )
                if( !_fstrnicmp(ptr,FormatOpt[i].ident,FormatOpt[i].len) )
                    break;

            if( i < FORMATOPTMAX )
            {   FileFormat = FormatOpt[i].format;
                ptr += FormatOpt[i].len;
            } else if( !_fstrnicmp(ptr,"sybyl",5) )
            {   FileFormat = FormatMol2;
                ptr += 5;
            } else if( !_fstrnicmp(ptr,"pdbnmr",6) )
            {   FileFormat = FormatNMRPDB;
                ptr += 6;
#ifdef CEXIOLIB
            } else if( !_fstrnicmp(ptr,"cex",3) )
            {   FileFormat = FormatCEX;
                ptr += 3;
#endif

            } else if( !_fstrnicmp(ptr,"script",6) )
            {   ptr += 6;
                while( *ptr && (*ptr==' ') )
                    ptr++;

                if( *ptr )
                {   dst = snamebuf;
                    while( *ptr && (*ptr!=' ') )
                        *dst++ = *ptr++;
                    *dst = '\0';
                } else return False;

#ifdef SOCKETS
            } else if( !_fstrnicmp(ptr,"port",4) )
            {   ptr += 4;
                while( *ptr && (*ptr==' ') )
                    ptr++;

                if( isdigit(*ptr) )
                {   ServerPort = (*ptr++)-'0';
                    while( isdigit(*ptr) )
                        ServerPort = (10*ServerPort)+(*ptr++)-'0';
                } else return False;
#endif
            } else if( !_fstrnicmp(ptr,"connect",7) )
            {   CalcBondsFlag = True;
                ptr += 7;
            } else if( !_fstrnicmp(ptr,"noconnect",9) )
            {   CalcBondsFlag = False;     
                ptr += 8;
            } else if( !_fstrnicmp(ptr,"insecure",8) )
            {   AllowWrite = True;
                ptr += 8;
            } else if( !_fstrnicmp(ptr,"secure",6) )
            {   AllowWrite = False;
                ptr += 6;
            } else return False;

        } else if( !*fnamebuf )
        {   dst = fnamebuf;
            if( *ptr == '"' )
            {   ptr++;
                while( *ptr && (*ptr!='"') )
                    *dst++ = *ptr++;
                if( *ptr=='"' ) 
                    ptr++;
            } else /* Unquoted! */
            {   while( *ptr && (*ptr!=' ') )
                    *dst++ = *ptr++;
            }
            *dst = '\0';
        } else return False;
    }
    return True;
}


int PASCAL WinMain( HINSTANCE hCurrent, HINSTANCE hPrevious,
                    LPSTR lpCmdLine, int nCmdShow )
{
    register FILE *fp;
    auto MSG event;

    hInstance = hCurrent;
    if( !hPrevious && !InitialiseApplication() )
        return False;

    InitDefaultValues();
    InitFileDialogBoxes();
    InitWindowsProfiles();
    ProcessOptions(lpCmdLine);

    /* Avoid Windows NT problems! */
    ReDrawFlag = RFInitial;
    CommandActive = True;

    if( !InitTerminal(hInstance) ||
        !OpenDisplay(hInstance,nCmdShow) )
       return False;

    WriteString("RasMol Molecular Renderer\n");
    WriteString("Roger Sayle, December 1998\n");
    WriteString("Version 2.6.4\n");

#ifdef EIGHTBIT
    WriteString("[8bit version]\n\n");
#else
#ifdef SIXTEENBIT
    WriteString("[16bit version]\n\n");
#else
    WriteString("[32bit version]\n\n");
#endif
#endif

    InitialiseCmndLine();
    InitialiseCommand(); 
    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();
    InitialiseRepres();
    InitialiseTMesh();
   
#ifdef SOCKETS
    OpenSocket(CanvWin);
#endif

    if( *fnamebuf && FetchFile(FileFormat,True,fnamebuf) )
        DefaultRepresentation();
    ResetCommandLine(1);

    LoadInitFile();
    if( *snamebuf )
    {   if( !(fp=fopen(snamebuf,"rb")) )
        {   if( CommandActive )
                WriteChar('\n');
            WriteString("Error: File '");
            WriteString(snamebuf);
            WriteString("' not found!\n");
            CommandActive = False;
        } else LoadScriptFile(fp,snamebuf);
    }

    RefreshScreen();
    if( !CommandActive )
        ResetCommandLine(0);

    DragAcceptFiles(CanvWin,TRUE);
    while( GetMessage(&event,NULL,0,0) )
    {   TranslateMessage(&event);
        DispatchMessage(&event);
    }
    DeleteObject(TermFont);
#ifdef DDEIPC
    CloseDDELinks();
#endif
#ifdef SOCKETS
    CloseSockets();
#endif
    CloseDisplay();

    return event.wParam;
}

