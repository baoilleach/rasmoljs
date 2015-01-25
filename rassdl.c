/* rastxt.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#include <SDL/SDL.h>
static SDL_Surface* Screen;



#ifndef sun386
#include <stdlib.h>
#endif
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <math.h>


#define RASMOL
#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"
#include "molecule.h"
#include "infile.h"
#include "abstree.h"
#include "transfor.h"
#include "cmndline.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "pixutils.h"
#include "outfile.h"


#ifdef IBMPC
#include <conio.h>
#endif
#ifdef UNIX
#include <signal.h>
#endif

#ifdef TERMIOS
#include <sys/types.h>
#include <sys/time.h>
#include <termios.h>

#ifdef esv
#include <sysv/unistd.h>
#else
#include <unistd.h>
#endif
#endif /* TERMIOS */



#define TwoPi           6.28318531
#define IsIdentChar(x)  ((isalnum(x))||((x)=='_')||((x)=='$'))

/* Either stdout or stderr */
#define OutFp stdout

#ifdef TERMIOS
static struct termios OrigTerm;
static struct termios IntrTerm;
#endif


static int InitialWide;
static int InitialHigh;

static char *FileNamePtr;
static char *ScriptNamePtr;
static int FileFormat;
static int ProfCount;



void WriteChar( int ch )
{   putc(ch,OutFp);
}


void WriteString( char *ptr )
{   fputs(ptr,OutFp);
}


static int InitTerminal( void )
{
#ifdef TERMIOS
    register int fd;
#endif

    setbuf(stdin,(char*)NULL);

#ifdef SIGTTIN
    signal(SIGTTIN,SIG_IGN);
#endif
#ifdef SIGTTOU
    signal(SIGTTOU,SIG_IGN);
#endif

#ifdef TERMIOS
    fd = fileno(stdin);
    if( !isatty(fd) )
        return False;

    if( tcgetattr(fd,&OrigTerm) < 0 )
        return False;

    IntrTerm = OrigTerm;
    IntrTerm.c_iflag |= IGNBRK|IGNPAR;
    IntrTerm.c_iflag &= ~(BRKINT|PARMRK|INPCK|IXON|IXOFF);
    IntrTerm.c_lflag &= ~(ICANON|ISIG|ECHO|ECHOE|ECHOK|ECHONL|NOFLSH);
    /* IntrTerm.c_lflag |= ISIG; */

    IntrTerm.c_cc[VMIN] = 1;
    IntrTerm.c_cc[VTIME] = 0;

#ifdef VSUSP /* Disable ^Z */
    IntrTerm.c_cc[VSUSP] = 0;
#endif

    tcsetattr(fd,TCSANOW,&IntrTerm);
#endif /* TERMIOS */
    return True;
}


static void ResetTerminal( void )
{
#ifdef TERMIOS
    register int fd;

    fd = fileno(stdin);
    if( isatty(fd) )
        tcsetattr(fd,TCSANOW,&OrigTerm);

#endif
}


void RasMolExit( void )
{
    WriteChar('\n');
    if( CommandActive )
        WriteChar('\n');
    ResetTerminal();
    exit(0);
}


void RasMolFatalExit( char *msg )
{
    WriteChar('\n');
    WriteString(msg);
    WriteChar('\n');
    WriteChar(0x07);
    ResetTerminal();
    exit(1);
}


void RasMolSignalExit( int i )
{
    UnusedArgument(i);

    RasMolFatalExit("*** Quit ***");
}


static int FetchCharacter( void )
{
#ifdef IBMPC
    return _getch();
#else
    return getc(stdin);
#endif
}


static int ReadCharacter( void )
{
    register int ch;

    ch = FetchCharacter();
#ifdef IBMPC
    if( ch && (ch!=0xe0) ) 
        return ch;

    ch = FetchCharacter();
    switch( ch )
    {   case('H'): return( 0x10 ); /* Up    */
        case('P'): return( 0x0e ); /* Down  */
        case('K'): return( 0x02 ); /* Left  */
        case('M'): return( 0x06 ); /* Right */
    }
#else
    if( ch != 0x1b )
        return ch;
    ch = FetchCharacter();
    if( (ch!='[') && (ch!='O') )
        return ch;

    ch = FetchCharacter();
    switch( ch )
    {   case('A'): return( 0x10 ); /* Up    */
        case('B'): return( 0x0e ); /* Down  */
        case('C'): return( 0x06 ); /* Right */
        case('D'): return( 0x02 ); /* Left  */
    }
#endif
    return 0;
}


static void LoadInitFile( void )
{
    register char *src,*dst;
    register FILE *initrc;
    register char *fname;
    char fnamebuf[128];

    fname = "RASMOL.INI";
    initrc = fopen(fname,"rb");
    if( !initrc && (src=(char*)getenv("HOME")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '\\';

        src = fname; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"rb");
    }

    if( !initrc && (src=(char*)getenv("RASMOLPATH")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '\\';

        src = "rasmolrc"; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"rb");
    }

    if( initrc )
        LoadScriptFile(initrc,fname);
}


int CreateImage( void )
{
    register Long size;
    
    if( FBuffer ) _ffree(FBuffer);
    size = (Long)XRange*YRange*sizeof(Pixel);
    FBuffer = (Pixel*)_fmalloc( size+32 );
    return( (int)FBuffer );
}


void TransferImage( void )
{
  Pixel __huge *src;
  int x, y;
  Uint32 *dst;

  if (SDL_MUSTLOCK(Screen)) SDL_LockSurface(Screen);
  
  dst = (Uint32 *)Screen->pixels;
  src = FBuffer;
  for( y=0; y<YRange; y++ ) {
    for( x=0; x<XRange; x++ ) {
      *dst++ = *src++;
    }
  }
  
  if (SDL_MUSTLOCK(Screen)) SDL_UnlockSurface(Screen);
  
  SDL_Flip(Screen);
}


int ClipboardImage( void )
{
    return False;
}


void ClearImage( void )
{
}


int PrintImage( void )
{
    return False;
}



void AllocateColourMap( void )
{
}


void UpdateScrollBars( void )
{
}


int LookUpColour( char *name, int *r, int *g, int *b )
{
    UnusedArgument(name);
    UnusedArgument(r);
    UnusedArgument(g);
    UnusedArgument(b);

    return False;
}


void SetMouseUpdateStatus( int bool )
{
    MouseUpdateStatus = bool;
}
                         
                         
void SetMouseCaptureStatus( int bool )
{
    MouseCaptureStatus = bool;
}
                         

void SetCanvasTitle( char *ptr )
{
    UnusedArgument(ptr);
}


void EnableMenus( int flag )
{
    UnusedArgument(flag);
}


void CloseDisplay( void )
{
}


void BeginWait( void )
{
}


void EndWait( void )
{
}


int OpenDisplay( int x, int y )
{
    register int i;

    for( i=0; i<8; i++ )
        DialValue[i] = 0.0;
    
    XRange = x;   WRange = XRange>>1;
    YRange = y;   HRange = YRange>>1;
    Range = MinFun(XRange,YRange);
    
    /* Initialise Palette! */
    for( i=0; i<256; i++ )
        ULut[i] = False;
    AllocateColourMap();
    return False;
}


void AdviseUpdate( int item )
{
    UnusedArgument(item);
}


void RefreshScreen( void )
{
    if( !UseSlabPlane )
    {   ReDrawFlag &= ~RFTransZ|RFSlab;
    } else ReDrawFlag &= ~RFTransZ;

    if( ReDrawFlag )
    {   if( ReDrawFlag & RFReSize )
            ReSizeScreen();

        if( ReDrawFlag & RFColour )
            DefineColourMap();

        if( Database )
        {   if( ReDrawFlag & RFApply ) 
                ApplyTransform();
            DrawFrame();
        }
        ReDrawFlag = 0;
    }
}



static void ProfileExecution( void )
{
    register long start,stop;
    register Real delta;
    register int i;

    delta = TwoPi/ProfCount;

    printf("Profiling Execution!\n");

    start = time((time_t *)NULL);
    for( i=0; i<ProfCount; i++ )
    {   DrawFrame();
        ReDrawFlag |= RFRotateY;
        DialValue[1] += delta;
        /* UpdateScrollBars(); */
        ApplyTransform();
    }

    stop = time((time_t *)NULL);
    fprintf(stderr,"Execution of %d frames\n",ProfCount);
    fprintf(stderr,"Duration = %ld seconds\n",stop-start);
}


static void InitDefaultValues( void )
{
    Interactive = False;

    FileNamePtr = NULL;
    ScriptNamePtr = NULL;
    InitialWide = DefaultWide;
    InitialHigh = DefaultHigh;
    ProfCount = 0;

    FileFormat = FormatPDB;
    CalcBondsFlag = True;
}


static void DisplayUsage( void )
{
    fputs("usage: rasdos [-script scriptfile] ",OutFp);
    fputs("[[-format] file]\n    formats: -pdb -nmrpdb ",OutFp);
    fputs("-mopac -mdl -mol2 -xyz -alchemy -charmm\n\n",OutFp);
    exit(1);
}


#define FORMATOPTMAX   15
static struct {
        char *ident;
        int format;
    } FormatOpt[FORMATOPTMAX] = { 
            { "alchemy",    FormatAlchemy  },
            { "biosym",     FormatBiosym   },
            { "cif",        FormatCIF      },
            { "charmm",     FormatCharmm   },
            { "fdat",       FormatFDAT     },
            { "gaussian",   FormatGaussian },
            { "macromodel", FormatMacroMod },
            { "mdl",        FormatMDL      },
            { "mmdb",       FormatMMDB     },
            { "mol2",       FormatMol2     },
            { "mopac",      FormatMOPAC    },
            { "nmrpdb",     FormatNMRPDB   },
            { "pdb",        FormatPDB      },
            { "shelx",      FormatSHELX    },
            { "xyz",        FormatXYZ      }
                                };

int InitializeSDL()
{  
  if(SDL_Init(SDL_INIT_VIDEO)<0) {
    printf("Failed SDL_Init %s\n", SDL_GetError());
    return False;
  }

  Screen=SDL_SetVideoMode(InitialWide,InitialHigh,32,SDL_ANYFORMAT);
  if(Screen==NULL) {
    printf("Failed SDL_SetVideoMode: %s\n", SDL_GetError());
    SDL_Quit();
    return False;
  }
  
  TransferImage();
  
  return True;
}

void MainLoop()
{
  SDL_Event event;
  int quit = False;
  
  while(quit==False) {
    while(SDL_PollEvent(&event)) {
      if(event.type == SDL_QUIT)
        quit = True;
    }
  }
  
}

int main( int argc, char *argv[] )
{
    register int done;

    InitDefaultValues();

    ReDrawFlag = 0;

    done = OpenDisplay(InitialWide,InitialHigh);

    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();
    InitialiseRepres();
    
    FileFormat = FormatPDB;
    done = FetchFile(FileFormat,True,"pdb1crn.ent");
    if( !done ) {
      RasMolFatalExit("Profile Error: Unable to read data file!");
    }
    DefaultRepresentation();
 
    SetRadiusValue(120);
    EnableWireframe(CylinderFlag,40);
    RefreshScreen();    

    InitializeSDL();
    MainLoop();
    
    SDL_Quit();
    return 0;
}

