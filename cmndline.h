/* cmndline.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

#define MAXBUFFLEN   256
#define MAXLINELEN   256

/* Button State Bits */
#define MMLft     0x01
#define MMMid     0x02
#define MMRgt     0x04
#define MMSft     0x08
#define MMCtl     0x10


#define MMRasMol   0x00
#define MMInsight  0x01
#define MMQuanta   0x02
#define MMSybyl    0x03


#ifdef CMNDLINE
char CurLine[MAXBUFFLEN];
int CurState,StateOption;
int CommandActive;
int MouseMode;

#else
extern char CurLine[MAXBUFFLEN];
extern int CurState,StateOption;
extern int CommandActive;
extern int MouseMode;
#endif

int ProcessCharacter( int );
void ResetCommandLine( int );
void InvalidateCmndLine( void );
void InitialiseCmndLine( void );

void SetMouseMode( int );
void ProcessMouseMove( int, int, int );
void ProcessMouseDown( int, int, int );
void ProcessMouseUp( int, int, int );

