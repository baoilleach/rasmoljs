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
#include "tmesh.h"


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
static int LabelOptFlag;
static int FileFormat;
static int ProfCount;

static int MainLoopTicks;
static unsigned int ProfileStartTime;
static int ProfileMainLoop;


void WriteChar( int ch )
{   
#ifdef __EMSCRIPTEN__
  char buf[2];
  buf[0] = ch; buf[1] = '\0';
  WriteString(buf);
#else
  putc(ch,OutFp);
#endif
}


void WriteString( char *ptr )
{
#ifdef __EMSCRIPTEN__
  char buf[MAXBUFFLEN];
  char *start = "handleEcho(\"";
  char *end = "\");";
  char *dst, *src;
  int i;
  
  dst = buf;
  src = start;
  while(*src != '\0')
  {
    *dst++ = *src++;
  }
  src = ptr;
  while(*src != '\0')
  {
    if (*src == '\n') { // Replace with literal \n
      *dst++ = '\\';
      *dst++ = 'n';
      src++;
    }
    else if (*src == '"') {
      *dst++ = '\\';
      *dst++ = '"';
      src++;
    }
    else
      *dst++ = *src++;
  }
  src = end;
  while(*src != '\0')
  {
    *dst++ = *src++;
  }
  *dst = '\0';
  //for (i=0; i<strlen(buf); i++)
  //  printf("%c (%d) ", buf[i], buf[i]);
  emscripten_run_script(buf);
#else
  fputs(ptr,OutFp);
#endif
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

unsigned int GetTickCount()
{
    struct timeval tv;
    if(gettimeofday(&tv, NULL) != 0)
        return 0;

    return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}

void TransferImage( void )
{
  Pixel __huge *src;
  int x, y;
  Uint32 *dst;

  if (SDL_MUSTLOCK(Screen)) SDL_LockSurface(Screen);
  
  dst = (Uint32 *)Screen->pixels;
  src = FBuffer;  
  memcpy(dst, src, YRange*XRange*sizeof(Uint32));

  if (SDL_MUSTLOCK(Screen)) SDL_UnlockSurface(Screen);
  
  SDL_Flip(Screen);
}


int ClipboardImage( void )
{
    return False;
}


void ClearImage( void )
{
  memset(FBuffer, 0xff000000, XRange*YRange*sizeof(Uint32));
  TransferImage();
}


int PrintImage( void )
{
    return False;
}



void AllocateColourMap( void )
{
  register int i;

  for( i=0; i<LutSize; i++ )
    if( ULut[i] )
      Lut[i] = RLut[i] | (GLut[i]<<8) | BLut[i]<<16 | 0xff000000;
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

int FetchEvent( int wait )
{
    /* Avoid compiler warning! */
    UnusedArgument(wait);
    return 0;
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

void HandleMenu( int hand )
{
    register int menu;
    register int item;
    register int mask;

    menu = hand>>8;
    item = hand&0xff;
    switch( menu )
    {   case(0):  /* File Menu */
                  switch( item )
                  {   case(1):  /* Open */
                                if( !Database )
                                    ResetCommandLine(2);
                                break;

                      case(2):  /* Save As */
                                if( Database )
                                    ResetCommandLine(4);
                                break;

                      case(3):  /* Close */
                                ZapDatabase();
                                break;

                      case(5):  /* Exit */
                                RasMolExit();
                                break;
                  } 
                  break;

        case(1):  /* Display Menu */
                  switch( item )
                  {   case(1):  /* Wireframe */
                                DisableSpacefill();
                                EnableWireframe(WireFlag,0);
                                SetRibbonStatus(False,0,0);
                                DisableBackbone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(2):  /* Backbone */
                                DisableSpacefill();
                                DisableWireframe();
                                SetRibbonStatus(False,0,0);
                                EnableBackbone(CylinderFlag,80);
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(3):  /* Sticks */
                                DisableSpacefill();
                                if( MainAtomCount<256 )
                                {   EnableWireframe(CylinderFlag,40);
                                } else EnableWireframe(CylinderFlag,80);
                                SetRibbonStatus(False,0,0);
                                ReDrawFlag |= RFRefresh;
                                DisableBackbone();
                                break;

                      case(4):  /* Spheres */
                                SetVanWaalRadius();
                                DisableWireframe();
                                SetRibbonStatus(False,0,0);
                                DisableBackbone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(5):  /* Ball & Stick */
                                SetRadiusValue(120);
                                EnableWireframe(CylinderFlag,40);
                                SetRibbonStatus(False,0,0);
                                DisableBackbone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(6):  /* Ribbons */
                                DisableSpacefill();
                                DisableWireframe();
                                SetRibbonStatus(True,RibbonFlag,0);
                                DisableBackbone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(7):  /* Strands */
                                DisableSpacefill();
                                DisableWireframe();
                                SetRibbonStatus(True,StrandFlag,0);
                                DisableBackbone();
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(8):  /* Cartoons */
                                DisableSpacefill();
                                DisableWireframe();
                                SetRibbonCartoons();
                                DisableBackbone();
                                ReDrawFlag |= RFRefresh;
                  }
                  break;

        case(2):  /* Colours Menu */
                  switch( item )
                  {   case(1):  /* Monochrome */
                                MonoColourAttrib(255,255,255);
                                ReDrawFlag |= RFColour;  break;
                      case(2):  /* CPK */
                                CPKColourAttrib();
                                ReDrawFlag |= RFColour;  break;
                      case(3):  /* Shapely */
                                ShapelyColourAttrib();
                                ReDrawFlag |= RFColour;  break;
                      case(4):  /* Group */
                                ScaleColourAttrib( GroupAttr );
                                ReDrawFlag |= RFColour;  break;
                      case(5):  /* Chain */
                                ScaleColourAttrib( ChainAttr );
                                ReDrawFlag |= RFColour;  break;
                      case(6):  /* Temperature */
                                ScaleColourAttrib( TempAttr );
                                ReDrawFlag |= RFColour;  break;
                      case(7):  /* Structure */
                                StructColourAttrib();
                                ReDrawFlag |= RFColour;  break;
                      case(8):  /* User */
                                UserMaskAttrib(MaskColourFlag);
                                ReDrawFlag |= RFColour;  break;
                  }
                  break;

        case(3):  /* Option Menu */
                  switch( item )
                  {   case(1):  /* Slabbing */
                                ReDrawFlag |= RFRefresh;
                                UseSlabPlane = !UseSlabPlane;
                                if( UseSlabPlane )
                                    UseShadow = False;
                                break;

                      case(2):  /* Hydrogens */
                                mask = NormAtomFlag;
                                if( HetaGroups )
                                    mask |= HeteroFlag;
                                Hydrogens = !Hydrogens;
                                ReDrawFlag |= RFRefresh;
                                      
                                if( Hydrogens )
                                {   SelectZone(mask|HydrogenFlag);
                                } else RestrictZone(mask);
                                break;

                      case(3):  /* Hetero Atoms */
                                mask = NormAtomFlag;
                                if( Hydrogens )
                                    mask |= HydrogenFlag;
                                HetaGroups = !HetaGroups;
                                ReDrawFlag |= RFRefresh;
                                
                                if( HetaGroups )
                                {   SelectZone(mask|HeteroFlag);
                                } else RestrictZone(mask);
                                break;

                      case(4):  /* Specular */
                                FakeSpecular = !FakeSpecular;
                                ReDrawFlag |= RFColour;
                                break;

                      case(5):  /* Shadows */
                                ReDrawFlag |= RFRefresh;
                                UseShadow = !UseShadow;
                                if( UseShadow )
                                {   ReviseInvMatrix();
                                    VoxelsClean = False;
                                    UseSlabPlane = False;
                                    ReAllocBuffers();
                                }
                                break;

                      case(6):  /* Stereo */
                                if( UseStereo )
                                {   SetStereoMode(False);
                                } else SetStereoMode(True);
                                ReDrawFlag |= RFRefresh;
                                break;

                      case(7):  /* Labels */
                                LabelOptFlag = !LabelOptFlag;
                                DefaultLabels(LabelOptFlag);
                                ReDrawFlag |= RFRefresh;
                                break;
                  }
                  break;

        case(4):  /* Export Menu */
                  ResetCommandLine(3);
                  StateOption = item;
                  break;

        case(5):  /* Help Menu */
                  break;
    }
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
            TransferImage();
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
    Interactive = True;

    FileNamePtr = NULL;
    ScriptNamePtr = NULL;
    InitialWide = DefaultWide;
    InitialHigh = DefaultHigh;
    ProfCount = 0;

    FileFormat = FormatPDB;
    CalcBondsFlag = True;
    
    ProfileMainLoop = False;
    MainLoopTicks = 0;
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

  Screen=SDL_SetVideoMode(InitialWide,InitialHigh,32,SDL_ANYFORMAT|SDL_RESIZABLE);
  if(Screen==NULL) {
    printf("Failed SDL_SetVideoMode: %s\n", SDL_GetError());
    SDL_Quit();
    return False;
  }
  SDL_WM_SetCaption("Rasmol 2.6.4", NULL);
  return True;
}

static void UpdateStatusKeys(int *status)
{
  Uint8 *keys;
  
#ifdef __EMSCRIPTEN__
  keys = SDL_GetKeyboardState(NULL);
#else
  keys = SDL_GetKeyState(NULL);
#endif
  if ( keys[SDLK_LSHIFT] == SDL_PRESSED )
    *status |= MMSft;
  if ( keys[SDLK_LCTRL] == SDL_PRESSED )
    *status |= MMCtl;
}

static int GetStatusFromMotion(SDL_MouseMotionEvent bevent)
{
  int status = 0;
  
  if (bevent.state & SDL_BUTTON_LMASK)
    status |= MMLft;
  if (bevent.state & SDL_BUTTON_MMASK)
    status |= MMMid;
  if (bevent.state & SDL_BUTTON_RMASK)
    status |= MMRgt;

  UpdateStatusKeys(&status);

  return status;    
}

static int GetStatusFromButton(SDL_MouseButtonEvent bevent)
{
  int status = 0;
  
  switch(bevent.button)
  {
    case SDL_BUTTON_LEFT:
      status |= MMLft;
      break;
    case SDL_BUTTON_RIGHT:
      status |= MMRgt;
      break;
    case SDL_BUTTON_MIDDLE:
      status |= MMMid;
      break;
    default:
      break;
  }

  UpdateStatusKeys(&status);

  return status;
}

void SetScreenSize(int width, int height)
{
  Screen = SDL_SetVideoMode( width, height, 32, SDL_ANYFORMAT | SDL_RESIZABLE );
  XRange = width;   WRange = XRange>>1;
  YRange = height;   HRange = YRange>>1;
  Range = MinFun(XRange, YRange);
  ReDrawFlag |= RFReSize;
  ReSizeScreen();
}

void MainLoop()
{
  if (ProfileMainLoop) {
    Real delta = TwoPi / 100;
    if (MainLoopTicks == 0)
      ProfileStartTime = GetTickCount();
    MainLoopTicks++;
    ReDrawFlag |= RFRotateY;
    DialValue[1] += delta;
    ApplyTransform();
    if (MainLoopTicks >= 100) {
      char buf[25];
      int fps = MainLoopTicks * 1000 / (GetTickCount() - ProfileStartTime);
      sprintf(buf, "Profile FPS: %d\n", fps);
      WriteString(buf);
      ProfileMainLoop = False;
    }
  }
    
  SDL_Event event;
  while(SDL_PollEvent(&event)) {
    switch(event.type)
    {
    case SDL_QUIT:
      SDL_Quit();
      exit(0);
      break;
    case SDL_MOUSEBUTTONDOWN:
      ProcessMouseDown(event.button.x, event.button.y, GetStatusFromButton(event.button));
      break;
    case SDL_MOUSEBUTTONUP:
      ProcessMouseUp(event.button.x, event.button.y, GetStatusFromButton(event.button));
      break;
    case SDL_MOUSEMOTION:
      ProcessMouseMove(event.button.x, event.button.y, GetStatusFromMotion(event.motion));
      break;
    case SDL_VIDEORESIZE:
      SetScreenSize(event.resize.w, event.resize.h);
      break;
    default:
      break;
    }
  }
  if(ReDrawFlag)
    RefreshScreen();
}

void HandleCommand(const char* command)
{
  if (strlen(command) >= MAXBUFFLEN-1)
    return;
    
  strcpy(CurLine, command);
  ExecuteCommand();
}

void SetProfileOn()
{
  ProfileMainLoop = True;
  MainLoopTicks = 0;
}

int main( int argc, char *argv[] )
{
    register FILE *fp;
    register int done;
    register char ch;

    InitDefaultValues();

    ReDrawFlag = 0;
    
    setbuf(OutFp,(char *)NULL);
    OpenDisplay(InitialWide,InitialHigh);
    InitTerminal();

    WriteString("RasMol Molecular Renderer\n");
    WriteString("Roger Sayle, December 1998\n");
    WriteString("Version 2.6.4\n");
    WriteString("[32bit version]\n\n");

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
    
    HandleCommand("load pdb1crn.ent");
    
    DefaultRepresentation();
 
    SetRadiusValue(120);
    EnableWireframe(CylinderFlag,40);
    
    InitializeSDL();
    RefreshScreen();

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(MainLoop, 0, True);
#else    
    ResetCommandLine(1);
    while(True)
      MainLoop();
#endif    
    return 0;
}

