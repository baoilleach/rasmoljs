/* cmndline.c
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
#ifndef sun386
#include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdio.h>

#define CMNDLINE
#include "cmndline.h"
#include "molecule.h"
#include "command.h"
#include "render.h"
#include "graphics.h"


/* Determine Mouse Sensitivity! */
#define IsClose(u,v) (((u)>=(v)-1) && ((u)<=(v)+1))


#define MM_NONE   0x00
#define MM_ROTX   0x01
#define MM_ROTY   0x02
#define MM_ROTZ   0x03
#define MM_TRNX   0x04
#define MM_TRNY   0x05
#define MM_ZOOM   0x06
#define MM_CLIP   0x07

#define MM_PICK   0x01
#define MM_NEXT   0x02
#define MM_LABL   0x03
#define MM_SELE   0x04



typedef struct {
        unsigned char dxfunc, dxinvrt;
        unsigned char dyfunc, dyinvrt;
        unsigned char click;
    } MouseMapping;


/* RasMol Mouse Bindings:
 *  XY Rotation:   Left (and any other) Button
 *  XY Translaion: Middle or Right Buttons
 *  Scale:         Shift (and Control) Left (and any other) Button
 *  Z Rotation:    Shift (and Control) Middle or Right Buttons
 *  Clipping:      Control Left (and any other) Button
 */
 
static MouseMapping MouseRasMol[32] = {
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /*                     */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft                 */
        { MM_TRNX, 0, MM_TRNY, 0, MM_PICK },  /*     Mid             */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft+Mid             */
        { MM_TRNX, 0, MM_TRNY, 0, MM_PICK },  /*         Rgt         */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft    +Rgt         */
        { MM_TRNX, 0, MM_TRNY, 0, MM_PICK },  /*     Mid+Rgt         */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft+Mid+Rgt         */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*             Sft     */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft        +Sft     */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*     Mid    +Sft     */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid    +Sft     */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt+Sft     */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft    +Rgt+Sft     */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*     Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*                 Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft            +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*     Mid        +Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid        +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt    +Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft    +Rgt    +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*     Mid+Rgt    +Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid+Rgt    +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*             Sft+Ctl */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft        +Sft+Ctl */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*     Mid    +Sft+Ctl */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid    +Sft+Ctl */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt+Sft+Ctl */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft    +Rgt+Sft+Ctl */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*     Mid+Rgt+Sft+Ctl */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft+Ctl */
    };        
           
           
/* MSI (Biosym) Insight Bindings:
 *     Needs to be verified!
 */
static MouseMapping MouseInsight[32] = {
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /*                     */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft                 */
        { MM_TRNX, 0, MM_TRNY, 0, MM_PICK },  /*     Mid             */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_PICK },  /* Lft+Mid             */
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /*         Rgt         */
        { MM_ROTZ, 0, MM_ROTZ, 0, MM_PICK },  /* Lft    +Rgt         */
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /*     Mid+Rgt         */
        { MM_CLIP, 0, MM_CLIP, 0, MM_PICK },  /* Lft+Mid+Rgt         */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*             Sft     */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft        +Sft     */
        { MM_TRNY, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid    +Sft     */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid    +Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt+Sft     */
        { MM_ROTZ, 0, MM_ROTZ, 0, MM_NEXT },  /* Lft    +Rgt+Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*     Mid+Rgt+Sft     */
        { MM_CLIP, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*                 Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft            +Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid        +Ctl */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid        +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt    +Ctl */
        { MM_ROTZ, 0, MM_ROTZ, 0, MM_NEXT },  /* Lft    +Rgt    +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*     Mid+Rgt    +Ctl */
        { MM_CLIP, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid+Rgt    +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*             Sft+Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft        +Sft+Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid    +Sft+Ctl */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_NEXT },  /* Lft+Mid    +Sft+Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt+Sft+Ctl */
        { MM_ROTZ, 0, MM_ROTZ, 0, MM_NEXT },  /* Lft    +Rgt+Sft+Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*     Mid+Rgt+Sft+Ctl */
        { MM_CLIP, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft+Ctl */
    };        


/* MSI Quanta Mouse Bindings:
 *     Needs to be verified!
 */
static MouseMapping MouseQuanta[32] = {
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /*                     */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft                 */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /*     Mid             */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft+Mid             */
        { MM_ROTZ, 0, MM_NONE, 0, MM_PICK },  /*         Rgt         */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft    +Rgt         */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /*     Mid+Rgt         */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /* Lft+Mid+Rgt         */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /*             Sft     */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft        +Sft     */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid    +Sft     */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid    +Sft     */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*         Rgt+Sft     */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft    +Rgt+Sft     */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*                 Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft            +Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /*     Mid        +Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft+Mid        +Ctl */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /*         Rgt    +Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft    +Rgt    +Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /*     Mid+Rgt    +Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /* Lft+Mid+Rgt    +Ctl */
        { MM_NONE, 0, MM_ZOOM, 0, MM_NEXT },  /*             Sft+Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft        +Sft+Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid    +Sft+Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid    +Sft+Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*         Rgt+Sft+Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft    +Rgt+Sft+Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid+Rgt+Sft+Ctl */
        { MM_NONE, 0, MM_CLIP, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft+Ctl */
    };        
                 
/* Tripos Sybyl 6.2 Mouse Bindings:
 *    Tripos Graphics Manual, July 1995
 *  Selection:      Left Button
 *  Scale:          Middle & Right Buttons
 *  XY Rotation:    Right Button
 *  Z Rotation:     Left & Right Buttons
 *  XY Translation: Middle Button
 *  Z Translation:  Left & Middle Buttons
 */
 
static MouseMapping MouseSybyl[32] = {
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /*                     */
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /* Lft                 */
        { MM_TRNX, 0, MM_TRNY, 0, MM_PICK },  /*     Mid             */
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /* Lft+Mid             */
        { MM_ROTY, 0, MM_ROTX, 0, MM_PICK },  /*         Rgt         */
        { MM_ROTZ, 0, MM_NONE, 0, MM_PICK },  /* Lft    +Rgt         */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_PICK },  /*     Mid+Rgt         */
        { MM_NONE, 0, MM_NONE, 0, MM_PICK },  /* Lft+Mid+Rgt         */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*             Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft        +Sft     */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid    +Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft+Mid    +Sft     */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /*         Rgt+Sft     */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /* Lft    +Rgt+Sft     */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_NEXT },  /*     Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft     */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*                 Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft            +Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid        +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft+Mid        +Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /*         Rgt    +Ctl */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /* Lft    +Rgt    +Ctl */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_NEXT },  /*     Mid+Rgt    +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft+Mid+Rgt    +Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /*             Sft+Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft        +Sft+Ctl */
        { MM_TRNX, 0, MM_TRNY, 0, MM_NEXT },  /*     Mid    +Sft+Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft+Mid    +Sft+Ctl */
        { MM_ROTY, 0, MM_ROTX, 0, MM_NEXT },  /*         Rgt+Sft+Ctl */
        { MM_ROTZ, 0, MM_NONE, 0, MM_NEXT },  /* Lft    +Rgt+Sft+Ctl */
        { MM_ZOOM, 0, MM_ZOOM, 0, MM_NEXT },  /*     Mid+Rgt+Sft+Ctl */
        { MM_NONE, 0, MM_NONE, 0, MM_NEXT },  /* Lft+Mid+Rgt+Sft+Ctl */
    };        



#define HISTSIZE    4096
#define HISTMASK    4095
static char HistBuff[HISTSIZE];
static int MinHist,MaxHist;
static int CurHist;

static char *CurPrompt;
static int CurPos,MaxPos;

static MouseMapping *MouseBinding;
static int PointX,PointY;
static int InitX,InitY;
static int HeldButton;



static void UpdateLine( void )
{
    register int i;

    for( i=CurPos; i<MaxPos; i++ )
        WriteChar(CurLine[i]);
    WriteChar(' ');
    for( i=MaxPos+1; i>CurPos; i-- )
        WriteChar(0x08);
}


static void CopyHistory( void )
{
    register int i;

    for( i=CurPos; i>0; i-- )
        WriteChar(0x08);
    for( i=0; i<MaxPos; i++ )
        WriteChar(' ');
    WriteChar(0x0D);
    WriteString(CurPrompt);

    CurPos = 0;
    if( (i=CurHist) != MaxHist )
        while( HistBuff[i] )
        {   CurLine[CurPos++] = HistBuff[i];
            WriteChar(HistBuff[i]);
            i = (i+1) & HISTMASK;
        }
    CurLine[CurPos] = 0;
    MaxPos = CurPos;
}


int ProcessCharacter( int ch )
{
    register int i;

    if( !ch ) return( False );

    if( IsPaused )
    {   if( (ch==0x04) || (ch==0x1a) )
        {   InterruptPauseCommand();
        } else ResumePauseCommand();
        return( False );
    }

    if( (ch>=' ') && (ch<='~') )
    {   if( MaxPos<MAXLINELEN )
        {   for( i=MaxPos; i>CurPos; i-- )
                CurLine[i] = CurLine[i-1];
            CurLine[CurPos++] = ch;
            CurLine[++MaxPos] = 0;

            WriteChar(ch);
            if( CurPos<MaxPos )
                UpdateLine();
        } else
            WriteChar(0x07);

    } else
        switch( ch )
        {    case( 0x7f ):  /* DEL and ^H */
             case( 0x08 ):  if( CurPos>0 )
                            {   for( i=CurPos; i<=MaxPos; i++ )
                                    CurLine[i-1] = CurLine[i];
                                CurPos--; MaxPos--;
                                WriteChar(0x08);
                                UpdateLine();
                            }
                            break;

             case( 0x04 ):  if( CurPos<MaxPos ) /* ^D */
                            {   for( i=CurPos; i<MaxPos; i++ )
                                    CurLine[i] = CurLine[i+1];
                                MaxPos--; UpdateLine();
                            }
                            break;

             case( 0x0d ):  /* ^M and ^J */
             case( 0x0a ):  WriteChar('\n');
                            if( MaxPos )
                                for( i=0; i<=MaxPos; i++ )
                                {    HistBuff[MaxHist] = CurLine[i];
                                     MaxHist=(MaxHist+1)&HISTMASK;
                                     if( MaxHist==MinHist )
                                     {   while( HistBuff[MinHist] )
                                             MinHist=(MinHist+1)&HISTMASK;
                                         MinHist=(MinHist+1)&HISTMASK;
                                     }
                                }
                            CommandActive = False;
                            return( True );

             case( 0x02 ):  if( CurPos>0 )  /* ^B */
                            {    WriteChar(0x08);
                                 CurPos--;
                            }
                            break;

             case( 0x06 ):  if( CurPos<MaxPos )  /* ^F */
                                WriteChar(CurLine[CurPos++]);
                            break;

             case( 0x01 ):  while( CurPos>0 )   /* ^A */
                            {    WriteChar(0x08);
                                 CurPos--;
                            }
                            break;

             case( 0x05 ):  while( CurPos<MaxPos )  /* ^E */
                                WriteChar(CurLine[CurPos++]);
                            break;

             case( 0x0c ):  WriteChar('\n');    /* ^L */
                            WriteString(CurPrompt);
                            for( i=0; i<MaxPos; i++ )
                                WriteChar(CurLine[i]);
                            for( i=CurPos; i<MaxPos; i++ )
                                WriteChar(0x08);
                            break;

             case( 0x10 ):  if( CurHist != MinHist ) /* ^P */
                            {   CurHist -= 2;
                                if( CurHist<0 )
                                    CurHist += HISTSIZE;
                                while( HistBuff[CurHist] )
                                    CurHist=CurHist?CurHist-1:HISTMASK;
                                CurHist = (CurHist+1)&HISTMASK;
                                CopyHistory();
                            }
                            break;

             case( 0x0e ):  if( CurHist != MaxHist ) /* ^N */
                            {   while( HistBuff[CurHist] )
                                    CurHist = (CurHist+1)&HISTMASK;
                                CurHist = (CurHist+1)&HISTMASK;
                                CopyHistory();
                            }
                            break;
        }
    return False;
}


void ResetCommandLine( int state )
{
    if( state )
    {   EnableMenus(state==1);
        switch( CurState=state )
        {   case(1):   CurPrompt="RasMol> ";            break;
            case(2):   CurPrompt="PDB file name:";      break;
            case(3):   CurPrompt="Image file name:";    break;
            case(4):   CurPrompt="Molecule file name:"; break;
        }
    }

    if( CommandActive )
        WriteChar('\n');
    CommandActive = True;
    WriteString(CurPrompt);

    CurHist = MaxHist;
    CurPos = MaxPos = 0;
    CurLine[0] = 0;
}


void InvalidateCmndLine( void )
{
    if( CommandActive )
        WriteChar('\n');
    CommandActive=False;
}


static void ClampDial( int dial, Real value )
{
    register Real temp;

    temp = DialValue[dial] + value;

    if( temp > 1.0 )
    {   DialValue[dial] = 1.0;
    } else if( temp < -1.0 )
    {   DialValue[dial] = -1.0;
    } else DialValue[dial] = temp;
}


static void WrapDial( int dial, Real value )
{
    register Real temp;

    temp = DialValue[dial] + value;
    while( temp < -1.0 )  temp += 2.0;
    while( temp > 1.0 )   temp -= 2.0;
    DialValue[dial] = temp;
}


void ProcessMouseDown( int x, int y, int stat )
{   
    register MouseMapping *map;
    
    map = &MouseBinding[stat];
    if( map->dxfunc || map->dyfunc )
        SetMouseCaptureStatus(True);
    InitX = PointX = x;
    InitY = PointY = y;
    HeldButton = True;
}


void ProcessMouseUp( int x, int y, int stat )
{                      
    register MouseMapping *map;
    
    map = &MouseBinding[stat];
    if( map->dxfunc || map->dyfunc )
        SetMouseCaptureStatus(False);

    if( !HeldButton )
        return;

    HeldButton = False;
                                  
    if( Database && IsClose(x,InitX) && IsClose(y,InitY) )
    {   switch( map->click )
        {   case(MM_PICK):  if( PickAtom(False,x,y) )
                            {   AdviseUpdate(AdvPickNumber);
                                AdviseUpdate(AdvPickAtom);
                                AdviseUpdate(AdvPickCoord);
                            }
                            break;
                   
            case(MM_NEXT):  if( PickAtom(True,x,y) )
                            {   AdviseUpdate(AdvPickNumber);
                                AdviseUpdate(AdvPickAtom);
                                AdviseUpdate(AdvPickCoord);
                            }
                            break;

        }
    }
}


static int BindingActive( int stat )
{
    return( MouseBinding[stat].dxfunc || MouseBinding[stat].dyfunc );
}


void SetMouseMode( int mode )
{
    register int stat;
    
    switch( mode )
    {   case(MMRasMol):  MouseBinding = MouseRasMol;  break;
        case(MMInsight): MouseBinding = MouseInsight; break;
        case(MMQuanta):  MouseBinding = MouseQuanta;  break;
        case(MMSybyl):   MouseBinding = MouseSybyl;   break;
    }

    /* Should also test for BindingActive(MMSft|MMCtl)! */
    stat = BindingActive(MMSft) || BindingActive(MMCtl);
    SetMouseUpdateStatus(stat);            
}


static void PerformMouseFunc( int func, int delta, int max )
{
    switch( func )
    {   case(MM_ROTX):  WrapDial( 0, (Real)(2*delta)/max );
                        ReDrawFlag |= RFRotateX;
                        break;
                        
        case(MM_ROTY):  WrapDial( 1, (Real)(2*delta)/max );
                        ReDrawFlag |= RFRotateY;
                        break;
                        
        case(MM_ROTZ):  WrapDial( 2, (Real)(2*delta)/max );
                        ReDrawFlag |= RFRotateZ;
                        break;
                        
        case(MM_ZOOM):  ClampDial( 3, (Real)(2*delta)/max );
                        ReDrawFlag |= RFZoom;
                        break;
                        
        case(MM_TRNX):  ClampDial( 4, (Real)delta/max );
                        ReDrawFlag |= RFTransX;
                        break;

        case(MM_TRNY):  ClampDial( 5, (Real)delta/max );
                        ReDrawFlag |= RFTransY;
                        break;

        case(MM_CLIP):  ClampDial( 7, (Real)delta/max );
                        ReDrawFlag |= RFSlab;
                        break;
    }     
}


void ProcessMouseMove( int x, int y, int stat )
{   
    register MouseMapping *map;
    register int dx,dy;
    
    map = &MouseBinding[stat];
    if( map->dxfunc || map->dyfunc )
    {   SetMouseCaptureStatus(True);
        if( !HeldButton )
        {   InitX = PointX = x;
            InitY = PointY = y;
            HeldButton = True;
            return;
        }
    
        if( IsClose(x,InitX) && IsClose(y,InitY) )
            return;                           
        
        if( map->dxinvrt )
        {      dx = PointX - x;
        } else dx = x - PointX;  
    
        if( map->dyinvrt )
        {      dy = PointY - y;
        } else dy = y - PointY;
    
        if( dx ) PerformMouseFunc(map->dxfunc,dx,XRange);
        if( dy ) PerformMouseFunc(map->dyfunc,dy,YRange);
    
        if( ReDrawFlag & (RFRotateX|RFRotateY) )
            UpdateScrollBars();
    
        PointX = x;
        PointY = y; 
    } else              
    {   SetMouseCaptureStatus(False);
        HeldButton = False;
    }
}


void InitialiseCmndLine( void )
{
    MaxHist = MinHist = 1;
    HistBuff[0] = 0;

    CommandActive = False;

    SetMouseMode(MMRasMol);
    HeldButton = False;
}

